from configparser import ConfigParser, ExtendedInterpolation
from shortuuid import uuid
import os
import re
import shutil
import subprocess

from ..annotate.constants import DEFAULTS as ANNOTATE_DEFAULTS
from ..annotate.main import ANNOTATION_PASS
from ..cluster import main as cluster_main
from ..constants import SUBCOMMAND, PROTOCOL
from ..summary.constants import DEFAULTS as SUMMARY_DEFAULTS
from ..tools import convert_tool_output
from ..util import stringify_args_to_command, mkdirp, output_tabbed_file, LOG, DEVNULL
from ..validate.constants import DEFAULTS as VALIDATE_DEFAULTS
from ..validate.main import VALIDATION_PASS
from .job import Job, ArrayJob, LogFile
from .scheduler import SlurmScheduler, TorqueScheduler, LocalScheduler, SgeScheduler
from .constants import JOB_STATUS, STD_OPTIONS

PROGNAME = shutil.which('mavis')
SHEBANG = '#!/bin/bash'
SCHEDULERS_BY_NAME = {sched.NAME: sched for sched in [SlurmScheduler, TorqueScheduler, LocalScheduler, SgeScheduler]}


def parse_run_time(filename):
    """
    parses the run time listed at the end of a file following mavis conventions
    """
    with open(filename, 'r') as fh:
        for line in fh.readlines()[::-1]:
            match = re.match(r'^\s*run time \(s\): (\d+)\s*$', line)
            if match:
                return int(match.group(1))
    return -1


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
                    LOG('converting input command:', command)
                    output_tabbed_file(convert_tool_output(
                        command[1:-2], command[-2], command[-1], log=LOG, assume_no_untemplated=assume_no_untemplated
                    ), output_filename)
                else:
                    command = ' '.join(command) + ' -o {}'.format(output_filename)
                    LOG('converting input command:')
                    LOG('>>>', command, time_stamp=False)
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
    args.update({k.replace('_filename', ''): v for k, v in config.reference.items() if k.endswith('_filename')})
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
    args.update({k.replace('_filename', ''): v for k, v in config.reference.items() if k.endswith('_filename')})
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
    args.update({k.replace('_filename', ''): v for k, v in config.reference.items() if k.endswith('_filename')})
    args.update(config.pairing.items())
    args.update(config.summary.items())
    args = {k:v for k, v in args.items() if k in allowed_args}
    return args


class Pipeline:

    def __init__(
        self,
        output_dir,
        scheduler,
        validations=None,
        annotations=None,
        pairing=None,
        summary=None,
        checker=None,
        batch_id='batch-{}'.format(uuid())
    ):
        """
        Args
            output_dir (str): path to main output directory for all mavis pipeline results
            scheduler (Scheduler): the class for interacting with a job scheduler
            validations (list of Job): list of validation jobs
            annotations (list of Job): list of annotation jobs
            pairing (Job): pairing job
            summary (Job): summary job
            batch_id (str): the batch id for this pipeline run. Used in avoinfing job name conflicts
        """
        self.scheduler = scheduler
        self.output_dir = output_dir
        self.validations = [] if validations is None else validations
        self.annotations = [] if annotations is None else annotations
        self.pairing = pairing
        self.summary = summary
        self.checker = checker
        self.batch_id = batch_id
        self.args = {}  # for local runs only, store config to be passed to MAVIS stage

    def write_submission_script(self, subcommand, job, args):
        """
        Args
            subcommand (SUBCOMMAND): the pipeline step this script will run
            job (Job): the job the script is for
            args (dict): arguments for the subcommand
        """
        LOG('writing:', job.script)
        with open(job.script, 'w') as fh:
            fh.write(SHEBANG + '\n')
            commands = ['{} {}'.format(PROGNAME, subcommand)] + stringify_args_to_command(args)
            fh.write(' \\\n\t'.join(commands) + '\n')

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
        pipeline = Pipeline(output_dir=config.output, scheduler=scheduler)

        annotation_output_files = []

        for libconf in config.libraries.values():
            base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
            LOG('setting up the directory structure for', libconf.library, 'as', base)
            libconf.inputs = run_conversion(config, libconf, conversion_dir)

            # run the cluster stage
            cluster_output = mkdirp(os.path.join(base, SUBCOMMAND.CLUSTER))  # creates the clustering output dir
            merge_args = {'batch_id': pipeline.batch_id, 'output': cluster_output}
            merge_args['split_only'] = SUBCOMMAND.CLUSTER in config.get('skip_stage', [])
            merge_args.update(config.reference.items())
            merge_args.update(config.cluster.items())
            merge_args.update(libconf.items())
            LOG('clustering', '(split only)' if merge_args['split_only'] else '')
            clustered_files = cluster_main.main(log_args=True, **merge_args)

            # make a validation job for each cluster file
            validate_job = None

            if SUBCOMMAND.VALIDATE not in config.skip_stage:
                mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE))
                for task_ident in range(1, len(clustered_files) + 1):
                    mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{}'.format(pipeline.batch_id, task_ident)))
                args = validate_args(config, libconf)
                args['output'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
                args['inputs'] = os.path.join(cluster_output, '{}-${}.tab'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))

                script_name = os.path.join(base, SUBCOMMAND.VALIDATE, 'submit.sh')
                job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
                job_options['memory_limit'] = config.schedule.validation_memory
                if libconf.protocol == PROTOCOL.TRANS:
                    job_options['memory_limit'] = config.schedule.trans_validation_memory

                if isinstance(scheduler, LocalScheduler):
                    validate_job = LocalArrayJob(
                        stage=SUBCOMMAND.VALIDATE,
                        tasks=len(clustered_files),
                        output_dir=os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{{task_ident}}'.format(pipeline.batch_id)),
                        name='MV_{}_{}'.format(libconf.library, pipeline.batch_id),
                        args=args,
                        **job_options
                    )
                else:
                    validate_job = ArrayJob(
                        stage=SUBCOMMAND.VALIDATE,
                        tasks=len(clustered_files),
                        output_dir=os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{{task_ident}}'.format(pipeline.batch_id)),
                        script=script_name,
                        name='MV_{}_{}'.format(libconf.library, pipeline.batch_id),
                        **job_options
                    )
                    pipeline.write_submission_script(SUBCOMMAND.VALIDATE, validate_job, args)
                pipeline.validations.append(validate_job)

            # make an annotation job for each validation/cluster job/file
            mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE))
            for task_ident in range(1, len(clustered_files) + 1):
                mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-{}'.format(pipeline.batch_id, task_ident)))
            args = annotate_args(config, libconf)
            args['output'] = os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
            args['inputs'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT), VALIDATION_PASS)

            script_name = os.path.join(base, SUBCOMMAND.ANNOTATE, 'submit.sh')
            job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            job_options['memory_limit'] = config.schedule.annotation_memory
            if validate_job:
                job_options['dependencies'] = [validate_job]
            annotate_job = ArrayJob(
                stage=SUBCOMMAND.ANNOTATE,
                tasks=len(clustered_files),
                script=script_name,
                name='MA_{}_{}'.format(libconf.library, pipeline.batch_id),
                output_dir=os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-{{task_ident}}'.format(pipeline.batch_id)),
                **job_options
            )
            pipeline.write_submission_script(SUBCOMMAND.ANNOTATE, annotate_job, args)
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
        args['annotations'] = config.reference.annotations_filename
        mkdirp(args['output'])
        args['inputs'] = annotation_output_files

        script_name = os.path.join(config.output, SUBCOMMAND.PAIR, 'submit.sh')

        pipeline.pairing = Job(
            SUBCOMMAND.PAIR,
            script=script_name,
            output_dir=args['output'],
            name='MP_{}'.format(pipeline.batch_id),
            dependencies=pipeline.annotations,
            **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
        )
        pipeline.write_submission_script(SUBCOMMAND.PAIR, pipeline.pairing, args)
        # set up the summary job
        args = summary_args(config)
        args['output'] = os.path.join(config.output, SUBCOMMAND.SUMMARY)
        mkdirp(args['output'])
        args['inputs'] = os.path.join(config.output, SUBCOMMAND.PAIR, 'mavis_paired*.tab')
        script_name = os.path.join(args['output'], 'submit.sh')

        pipeline.summary = Job(
            stage=SUBCOMMAND.SUMMARY,
            name='MS_{}'.format(pipeline.batch_id),
            output_dir=args['output'],
            script=script_name,
            dependencies=[pipeline.pairing],
            **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
        )
        pipeline.write_submission_script(SUBCOMMAND.SUMMARY, pipeline.summary, args)
        return pipeline

    def _job_status(self, job, submit=False, log=DEVNULL):
        """
        report information regarding a particular job status
        """
        run_time = -1
        if not job.job_ident and submit:
            self.scheduler.submit(job)
        if job.job_ident:
            log('{} ({}) is {}'.format(job.name, job.job_ident, job.status))
        else:
            log('{} is {}'.format(job.name, job.status))
        if job.status == JOB_STATUS.COMPLETED:
            if isinstance(job, ArrayJob):
                for task in job.task_list:
                    if not os.path.exists(task.complete_stamp()):
                        with log.indent() as log:
                            log('complete stamp is expected but does not exist')
                            log(task.complete_stamp())
                    else:
                        run_time = max(run_time, parse_run_time(task.complete_stamp()))
            elif not os.path.exists(job.complete_stamp()):
                with log.indent() as log:
                    log('complete stamp is expected but does not exist')
                    log(job.complete_stamp())
            else:
                run_time = max(run_time, parse_run_time(job.complete_stamp()))
            if run_time >= 0:
                log('run time: {}'.format(run_time), indent_level=1)
        else:
            if isinstance(job, ArrayJob):
                tasks_by_status = {}
                for task in job.task_list:
                    tasks_by_status.setdefault(task.status, []).append(task)
                for status, tasks in tasks_by_status.items():
                    comments = set([t.status_comment for t in tasks if t.status_comment])
                    with log.indent() as log:
                        log('{} tasks are {}'.format(len(tasks), status))
                        for comment in comments:
                            log('comment:', comment, indent_level=1)
            elif job.status not in {JOB_STATUS.PENDING, JOB_STATUS.NOT_SUBMITTED, JOB_STATUS.SUBMITTED}:
                content = LogFile.parse(job.logfile())
                log('{}: {}'.format(content.status, content.message), indent_level=1)

        return run_time

    def check_status(self, submit=False, log=DEVNULL):
        """
        Check all jobs for completetion. Report any failures, etc.

        Args
            submit (bool): submit any pending jobs
        """
        # update the information for all jobs where possible
        total_run_time = 0
        jobs_not_complete = 0
        jobs_with_errors = 0

        for job in self.validations + self.annotations + [self.pairing, self.summary]:
            self.scheduler.update_info(job)
        log('validate', time_stamp=True)
        for job in self.validations:
            run_time = self._job_status(job, submit, log.indent())
            if job.status == JOB_STATUS.COMPLETED:
                if run_time >= 0:
                    total_run_time += run_time
            elif job.status in {JOB_STATUS.ERROR, JOB_STATUS.FAILED, JOB_STATUS.CANCELLED}:
                jobs_with_errors += 1
            else:
                jobs_not_complete += 1

        log('annotate', time_stamp=True)
        for job in self.annotations:
            self._job_status(job, submit, log.indent())
            if job.status == JOB_STATUS.COMPLETED:
                if run_time >= 0:
                    total_run_time += run_time
            elif job.status in {JOB_STATUS.ERROR, JOB_STATUS.FAILED, JOB_STATUS.CANCELLED}:
                jobs_with_errors += 1
            else:
                jobs_not_complete += 1

        log('pairing', time_stamp=True)
        run_time = self._job_status(self.pairing, submit, log.indent())
        if self.pairing.status == JOB_STATUS.COMPLETED:
            if run_time >= 0:
                total_run_time += run_time
        elif self.pairing.status in {JOB_STATUS.ERROR, JOB_STATUS.FAILED, JOB_STATUS.CANCELLED}:
            jobs_with_errors += 1
        else:
            jobs_not_complete += 1

        log('summary', time_stamp=True)
        run_time = self._job_status(self.summary, submit, log.indent())
        if self.summary.status == JOB_STATUS.COMPLETED:
            if run_time >= 0:
                total_run_time += run_time
        elif self.summary.status in {JOB_STATUS.ERROR, JOB_STATUS.FAILED, JOB_STATUS.CANCELLED}:
            jobs_with_errors += 1
        else:
            jobs_not_complete += 1

        if jobs_not_complete + jobs_with_errors == 0:
            log('parllel run time:', total_run_time)
            return 0
        elif not jobs_with_errors:
            return 1
        else:
            return 2

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
        cast = {'None': None, 'False': False, 'True': True}

        jobs = {}
        for sec in parser.sections():
            if sec != 'general':
                section = {}
                for attr, value in parser[sec].items():
                    if attr in ['dependencies', 'inputs', 'outputs'] and value:
                        section[attr] = re.split(r'[;\s]+', value)
                    elif value == 'None':
                        section[attr] = None
                    elif value in cast:
                        value = cast[value]
                    else:
                        section[attr] = value
                if 'tasks' in parser[sec]:
                    jobs[sec] = ArrayJob(**section)
                else:
                    jobs[sec] = Job(**section)

        for job in jobs.values():
            for i, prior_job_name in enumerate(job.dependencies):
                job.dependencies[i] = jobs[prior_job_name]

        pipeline = cls(
            output_dir=parser['general']['output_dir'],
            scheduler=SCHEDULERS_BY_NAME[parser['general']['scheduler']],
            batch_id=parser['general']['batch_id']
        )

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
            elif job.stage == SUBCOMMAND.CHECKER:
                if pipeline.checker:
                    raise ValueError('mavis pipeline expects a single checker job')
                pipeline.checker = job
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
            'output_dir': self.output_dir,
            'scheduler': self.scheduler.NAME
        }

        for job in [self.summary, self.pairing] + self.validations + self.annotations:
            parser[job.name] = {k: re.sub(r'\$', '$$', v) for k, v in job.flatten().items()}

        with open(filename, 'w') as configfile:
            parser.write(configfile)
