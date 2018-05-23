from shortuuid import uuid
import os
import re
from configparser import ConfigParser, ExtendedInterpolation
import subprocess

from ..annotate.constants import DEFAULTS as ANNOTATE_DEFAULTS
from ..annotate.main import ANNOTATION_PASS
from ..cluster import main as cluster_main
from ..cluster.constants import DEFAULTS as CLUSTER_DEFAULTS
from ..constants import SUBCOMMAND, PROTOCOL
from ..summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from ..tools import convert_tool_output
from ..util import stringify_args_to_command, mkdirp, output_tabbed_file, log
from ..validate.constants import DEFAULTS as VALIDATE_DEFAULTS
from ..validate.main import VALIDATION_PASS
from .job import Job, ArrayJob, SlurmScheduler, PbsScheduler, LocalScheduler, SgeScheduler, STD_OPTIONS

PROGNAME = 'mavis'
SHEBANG = '#!/bin/bash'
SCHEDULERS_BY_NAME = {sched.NAME: sched for sched in [SlurmScheduler, PbsScheduler, LocalScheduler, SgeScheduler]}


def run_conversion(config, libconf, conversion_dir, assume_no_untemplated=True):
    """
    Converts files if not already converted. Returns a list of filenames
    """
    inputs = []
    # run the conversions
    for input_file in libconf.inputs:
        output_filename = os.path.join(conversion_dir, input_file + '.tab')
        if input_file in config.convert:
            if not os.path.exists(output_filename):
                command = config.convert[input_file]
                if command[0] == 'convert_tool_output':  # convert_tool_output FILEPATH [FILEPATH...] TOOL stranded
                    log('converting input command:', command)
                    output_tabbed_file(convert_tool_output(
                        command[1:-2], command[-2], command[-1], log=log, assume_no_untemplated=assume_no_untemplated
                    ), output_filename)
                else:
                    command = ' '.join(command) + ' -o {}'.format(output_filename)
                    log('converting input command:')
                    log('>>>', command, time_stamp=False)
                    subprocess.check_output(command, shell=True)
            inputs.append(output_filename)
        else:
            inputs.append(input_file)
    return inputs


def validate_args(config, libconf):
    """
    returns the mavis command for running the validate step

    writes the bash script for running the validation job

    creates a job instance with input/output paths

    """
    allowed_args = [
        'masking',
        'reference_genome',
        'aligner_reference',
        'library',
        'bam_file',
        'protocol',
        'read_length',
        'stdev_fragment_size',
        'median_fragment_size',
        'strand_specific',
        'annotations'
    ] + list(VALIDATE_DEFAULTS.keys())

    # overwrite args in order of increasing specificity
    args = {}
    args.update(VALIDATE_DEFAULTS.items())
    args.update(config.reference.items())
    args.update(config.validate.items())
    args.update(libconf.items())
    args = {k:v for k, v in args.items() if k in allowed_args}
    return args

def annotate_args(config, libconf):
    allowed_args = [
        'reference_genome',
        'template_metadata',
        'masking',
        'annotations',
        'min_orf_size',
        'max_orf_cap',
        'library',
        'protocol',
        'min_domain_mapping_match',
        'domain_name_regex_filter',
        'max_proximity'
    ] + list(ANNOTATE_DEFAULTS.keys())
    args = {}
    args.update(ANNOTATE_DEFAULTS.items())
    args.update(config.reference.items())
    args.update(config.cluster.items())
    args.update(config.illustrate.items())
    args.update(config.annotate.items())
    args.update(libconf.items())
    args = {k:v for k, v in args.items() if k in allowed_args}
    return args

def summary_args(config):
    allowed_args = [
        'flanking_call_distance',
        'split_call_distance',
        'contig_call_distance',
        'spanning_call_distance',
        'dgv_annotation',
        'annotations'
    ] + list(SUMMARY_DEFAULTS.keys())
    args = {}
    args.update(config.reference.items())
    args.update(config.pairing.items())
    args.update(config.summary.items())
    args = {k:v for k, v in args.items() if k in allowed_args}
    return args


class Pipeline:

    def __init__(
        self,
        outputdir,
        validations=None,
        annotations=None,
        pairing=None,
        summary=None,
        batch_id='batch-{}'.format(uuid())
    ):
        """
        Args
            outputdir (str): path to main output directory for all mavis pipeline results
            validations (list of Job): list of validation jobs
            annotations (list of Job): list of annotation jobs
            pairing (Job): pairing job
            summary (Job): summary job
            batch_id (str): the batch id for this pipeline run. Used in avoinfing job name conflicts
        """
        self.outputdir = outputdir
        self.validations = [] if validations is None else validations
        self.annotations = [] if annotations is None else annotations
        self.pairing = pairing
        self.summary = summary
        self.batch_id = batch_id

    @classmethod
    def build(cls, config):
        """
        Args
            config (MavisConfig): the main configuration. Note this is the config after all reference inputs have been loaded
        Returns
            Pipeline: the pipeline instance with job dependencies information etc.
        """
        conversion_dir = mkdirp(os.path.join(config.output, 'converted_inputs'))
        config.output = os.path.abspath(config.output)
        if config.schedule.scheduler not in SCHEDULERS_BY_NAME:
            raise NotImplementedError('unsupported scheduler', config.schedule.scheduler, list(SCHEDULERS_BY_NAME.keys()))

        scheduler = SCHEDULERS_BY_NAME[config.schedule.scheduler]()
        pipeline = Pipeline(outputdir=config.output)

        annotation_output_files = []

        for libconf in config.libraries.values():
            base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
            log('setting up the directory structure for', libconf.library, 'as', base)
            libconf.inputs = run_conversion(config, libconf, conversion_dir)

            # run the cluster stage
            cluster_output = mkdirp(os.path.join(base, SUBCOMMAND.CLUSTER))  # creates the clustering output dir
            merge_args = {'batch_id': pipeline.batch_id, 'output': cluster_output}
            merge_args['split_only'] = SUBCOMMAND.CLUSTER in config.get('skip_stage', [])
            merge_args.update(config.reference.items())
            merge_args.update(config.cluster.items())
            merge_args.update(libconf.items())
            log('clustering', '(split only)' if merge_args['split_only'] else '')
            clustered_files = cluster_main.main(log_args=True, **merge_args)

            # make a validation job for each cluster file
            validate_job = None

            if SUBCOMMAND.VALIDATE not in config.skip_stage:
                mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE))
                args = validate_args(config, libconf)
                args['output'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
                args['inputs'] = os.path.join(cluster_output, '{}-${}.tab'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))


                script_name = os.path.join(base, SUBCOMMAND.VALIDATE, 'submit.sh')
                log('writing:', script_name)
                with open(script_name, 'w') as fh:
                    fh.write(SHEBANG + '\n\n')
                    commands = ['{} {}'.format(PROGNAME, SUBCOMMAND.VALIDATE)] + stringify_args_to_command(args)
                    fh.write('\\\n\t'.join(commands) + '\n')

                job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
                job_options['memory_limit'] = config.schedule.validation_memory
                if libconf.protocol == PROTOCOL.TRANS:
                    job_options['memory_limit'] = config.schedule.trans_validation_memory

                validate_job = ArrayJob(
                    stage=SUBCOMMAND.VALIDATE,
                    tasks=len(clustered_files),
                    script=script_name,
                    name='MV_{}_{}'.format(libconf.library, pipeline.batch_id),
                    **job_options
                )
                pipeline.validations.append(validate_job)

            # make an annotation job for each validation/cluster job/file
            mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE))
            args = annotate_args(config, libconf)
            args['output'] = os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
            args['inputs'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT), VALIDATION_PASS)
            script_name = os.path.join(base, SUBCOMMAND.ANNOTATE, 'submit.sh')
            log('writing:', script_name)
            with open(script_name, 'w') as fh:
                fh.write(SHEBANG + '\n\n')
                commands = ['{} {}'.format(PROGNAME, SUBCOMMAND.ANNOTATE)] + stringify_args_to_command(args)
                fh.write('\\\n\t'.join(commands) + '\n')
            job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            job_options['memory_limit'] = config.schedule.annotation_memory
            if validate_job:
                job_options['dependencies'] = [validate_job]
            annotate_job = ArrayJob(
                stage=SUBCOMMAND.ANNOTATE,
                tasks=len(clustered_files),
                script=script_name,
                name='MA_{}_{}'.format(libconf.library, pipeline.batch_id),
                **job_options
            )
            pipeline.annotations.append(annotate_job)

            # add the expected output file names for input to pairing
            for taskid in range(1, len(clustered_files) + 1):
                fname = os.path.join(args['output'], ANNOTATION_PASS)
                fname = re.sub(r'\${}'.format(scheduler.ENV_TASK_IDENT), str(taskid), fname)
                annotation_output_files.append(fname)

        # set up the pairing job
        args = {}
        args.update(config.pairing.items())
        args['output'] = os.path.join(config.output, SUBCOMMAND.PAIR)
        mkdirp(args['output'])
        args['inputs'] = '\\\n\t'.join(annotation_output_files)
        args = ['{} {}'.format(PROGNAME, SUBCOMMAND.PAIR)] + stringify_args_to_command(args)
        script_name = os.path.join(config.output, SUBCOMMAND.PAIR, 'submit.sh')
        log('writing:', script_name)
        with open(script_name, 'w') as fh:
            fh.write(SHEBANG + '\n\n')
            fh.write('\\\n\t'.join(args))
        pipeline.pairing = Job(
            SUBCOMMAND.PAIR,
            script=script_name,
            name='MP_{}'.format(pipeline.batch_id),
            dependencies=pipeline.annotations,
            **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
        )
        # set up the summary job
        args = summary_args(config)
        args['output'] = os.path.join(config.output, SUBCOMMAND.SUMMARY)
        mkdirp(args['output'])
        args['inputs'] = os.path.join(config.output, SUBCOMMAND.PAIR, 'mavis_paired*.tab')
        script_name = os.path.join(args['output'], 'submit.sh')
        log('writing:', script_name)
        with open(script_name, 'w') as fh:
            fh.write(SHEBANG + '\n\n')
            commands = ['{} {}'.format(PROGNAME, SUBCOMMAND.SUMMARY)] + stringify_args_to_command(args)
            fh.write('\\\n\t'.join(commands) +  '\n')
        pipeline.summary = Job(
            stage=SUBCOMMAND.SUMMARY,
            name='MS_{}'.format(pipeline.batch_id),
            script=script_name,
            dependencies=[pipeline.pairing],
            **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
        )
        return pipeline

    def check_status(self):
        pass

    @classmethod
    def read_build_file(cls, filepath):
        """
        read the configuration file which stored the build information concerning jobs and dependencies

        Args
            filepath (str): path to the input config file
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError('File does not exist: {}'.format(filepath))
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser.read(filepath)

        jobs = {}
        print(parser)
        for sec in parser.sections():
            print(parser[sec])
            if sec != 'general':
                section = {}
                for attr, value in parser[sec].items():
                    if attr in ['dependencies', 'inputs', 'outputs'] and value:
                        section[attr] = re.split(r'[;\s]+', value)
                    else:
                        section[attr] = value
                if 'tasks' in parser[sec]:
                    jobs[sec] = ArrayJob(**section)
                else:
                    jobs[sec] = Job(**section)

        for job in jobs.values():
            print(job)
            print(job.dependencies)
            for i, prior_job_name in enumerate(job.dependencies):
                job.dependencies[i] = jobs[prior_job_name]

        pipeline = cls(**dict(parser['general'].items()))

        for job in jobs.values():
            if job.stage == SUBCOMMAND.VALIDATE:
                pipeline.validations.append(job)
            elif job.stage == SUBCOMMAND.ANNOTATE:
                pipeline.annotations.append(job)
            elif job.stage == SUBCOMMAND.PAIR:
                if pipeline.pairing:
                    raise ValueError('mavis pipeline expects a single pairing job')
                pipeline.pairing = job
            elif job.stage == SUBCOMMAND.SUMMARY:
                if pipeline.summary:
                    raise ValueError('mavis pipeline expects a single summary job')
                pipeline.summary = job
            else:
                raise NotImplementedError('unexpected job stage for MAVIS pipeline: {}'.format(job.stage), job)

        return pipeline


    def write_build_file(self, filename):
        """
        write the build.cfg file for the current pipeline. This is the file used in re-loading the pipeline
        to check the status and report failures, etc. later.

        Args
            filename (str): path to the output config file
        """
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser['general'] = {
            'batch_id': self.batch_id,
            'outputdir': self.outputdir
        }

        for job in [self.summary, self.pairing] + self.validations + self.annotations:
            parser[job.name] = job.flatten()

        with open(filename, 'w') as configfile:
            parser.write(configfile)
