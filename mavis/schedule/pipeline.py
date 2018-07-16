from configparser import ConfigParser, ExtendedInterpolation
import os
import re
import shutil
import subprocess

from shortuuid import uuid

from ..cluster import constants as _CLUSTER
from ..constants import SUBCOMMAND, PROTOCOL, EXIT_ERROR, EXIT_OK, EXIT_INCOMPLETE
from ..tools import convert_tool_output
from ..util import mkdirp, output_tabbed_file, LOG, DEVNULL
from ..validate import constants as _VALIDATE
from ..annotate import constants as _ANNOTATE
from ..annotate import file_io as _file_io
from ..summary import constants as _SUMMARY
from .job import Job, ArrayJob, LogFile, TorqueArrayJob
from .scheduler import SlurmScheduler, TorqueScheduler, SgeScheduler, consecutive_ranges
from .local import LocalJob, LocalScheduler
from .constants import JOB_STATUS, STD_OPTIONS, OPTIONS, SCHEDULER

PROGNAME = shutil.which('mavis')
SHEBANG = '#!/bin/bash'
SCHEDULERS_BY_NAME = {sched.NAME: sched for sched in [SlurmScheduler, TorqueScheduler, LocalScheduler, SgeScheduler]}


def stringify_args_to_command(args):
    """
    takes a list of arguments and prepares them for writing to a bash script
    """
    command = []
    for argname, value in args.items():
        if isinstance(value, _file_io.ReferenceFile):
            value = value.name
        if isinstance(value, str):
            command.append('--{} "{}"'.format(argname, value))
        else:
            try:
                value = ' '.join([str(v) for v in value])
            except TypeError:
                pass
            command.append('--{} {}'.format(argname, value))
    return command


def parse_run_time(filename):
    """
    parses the run time listed at the end of a file following mavis conventions
    """
    with open(filename, 'r') as fh:
        content = fh.read().strip()
        for line in [l.strip() for l in content.split('\n')][::-1]:
            match = re.match(r'^\s*run time \(s\): (\d+)\s*$', line)  # older style complete stamp
            if match:
                return int(match.group(1))
            match = re.search(r'start:\s*(\d+)\s*end:\s*(\d+)', line)
            if match:
                return int(match.group(2)) - int(match.group(1))
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
    Pull arguments from the main config and library specific config to pass to validate

    Args:
        config (MavisConfig): the main program config
        libconf (LibraryConfig): library specific configuration
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
    ] + list(_VALIDATE.DEFAULTS.keys())

    # overwrite args in order of increasing specificity
    args = {}
    args.update(_VALIDATE.DEFAULTS.items())
    args.update({k: v.name for k, v in config.reference.items()})
    args.update(config.validate.items())
    args.update(libconf.items())
    args = {k: v for k, v in args.items() if k in allowed_args}
    return args


def annotate_args(config, libconf):
    """
    Pull arguments from the main config and library specific config to pass to annotate

    Args:
        config (MavisConfig): the main program config
        libconf (LibraryConfig): library specific configuration
    """
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
    ] + list(_ANNOTATE.DEFAULTS.keys())
    args = {}
    args.update(_ANNOTATE.DEFAULTS.items())
    args.update({k: v.name for k, v in config.reference.items()})
    args.update(config.cluster.items())
    args.update(config.illustrate.items())
    args.update(config.annotate.items())
    args.update(libconf.items())
    args = {k: v for k, v in args.items() if k in allowed_args}
    return args


def summary_args(config):
    """
    Pull arguments from the main config and library specific config to pass to summary

    Args:
        config (MavisConfig): the main program config
        libconf (LibraryConfig): library specific configuration
    """
    allowed_args = [
        'flanking_call_distance',
        'split_call_distance',
        'contig_call_distance',
        'spanning_call_distance',
        'dgv_annotation',
        'annotations'
    ] + list(_SUMMARY.DEFAULTS.keys())
    args = {}
    args.update({k: v.name for k, v in config.reference.items()})
    args.update(config.pairing.items())
    args.update(config.summary.items())
    args = {k: v for k, v in args.items() if k in allowed_args}
    return args


def cluster_args(config, libconf):
    """
    Pull arguments from the main config and library specific config to pass to cluster

    Args:
        config (MavisConfig): the main program config
        libconf (LibraryConfig): library specific configuration
    """
    allowed_args = [
        'masking',
        'annotations',
        'library',
        'protocol',
        'disease_status',
        'strand_specific'
    ] + list(_CLUSTER.DEFAULTS.keys())
    args = {}
    args.update(_CLUSTER.DEFAULTS.items())
    args.update({k: v.name for k, v in config.reference.items()})
    args.update(config.cluster.items())
    args.update(config.illustrate.items())
    args.update(config.annotate.items())
    args.update(libconf.items())
    args = {k: v for k, v in args.items() if k in allowed_args}
    return args


class Pipeline:
    ERROR_STATES = {JOB_STATUS.ERROR, JOB_STATUS.FAILED, JOB_STATUS.CANCELLED, JOB_STATUS.UNKNOWN, JOB_STATUS.NOT_SUBMITTED}

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
        Args:
            output_dir (str): path to main output directory for all mavis pipeline results
            scheduler (Scheduler): the class for interacting with a job scheduler
            validations (:class:`list` of :class:`Job`): list of validation jobs
            annotations (:class:`list` of :class:`Job`): list of annotation jobs
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

    def write_submission_script(self, subcommand, job, args, aligner_path=None):
        """
        Args:
            subcommand (SUBCOMMAND): the pipeline step this script will run
            job (Job): the job the script is for
            args (dict): arguments for the subcommand
        """
        LOG('writing:', job.script, time_stamp=True)
        with open(job.script, 'w') as fh:
            fh.write("""{shebang}
{aligner_path}
cd {cwd}
START_TIME=$(date +%s)\n\n""".format(
                shebang=SHEBANG,
                aligner_path='export PATH={}:$PATH'.format(os.path.dirname(aligner_path)) if aligner_path else '',
                cwd=os.getcwd()
            ))
            commands = [PROGNAME, subcommand] + stringify_args_to_command(args)
            fh.write(' \\\n\t'.join(commands) + '\n\n')
            fh.write("""
code=$?

if [ "$code" -ne "0" ]
then
    exit $code
fi

END_TIME=$(date +%s)

echo "start: $START_TIME end: $END_TIME" > {}/MAVIS-${}.COMPLETE

            """.format(
                args['output'],
                self.scheduler.ENV_JOB_IDENT if not isinstance(job, ArrayJob) else self.scheduler.ENV_ARRAY_IDENT))

    @classmethod
    def format_args(cls, subcommand, args):
        command = [subcommand]
        for arg, val in args.items():
            command.append('--{}'.format(arg))
            if isinstance(val, str):
                command.append(val)
            else:
                try:
                    command.extend(iter(val))
                except TypeError:
                    command.append(val)
        return [str(v) for v in command]

    @classmethod
    def build(cls, config):
        """
        Args:
            config (MavisConfig): the main configuration. Note this is the config after all reference inputs have been loaded
        Returns:
            Pipeline: the pipeline instance with job dependencies information etc.
        """
        from ..main import main as _main
        conversion_dir = mkdirp(os.path.join(config.output, 'converted_inputs'))
        config.output = os.path.abspath(config.output)
        if config.schedule.scheduler not in SCHEDULERS_BY_NAME:
            raise NotImplementedError('unsupported scheduler', config.schedule.scheduler, list(SCHEDULERS_BY_NAME.keys()))

        scheduler = SCHEDULERS_BY_NAME[config.schedule.scheduler](
            config.schedule.get('concurrency_limit', OPTIONS.concurrency_limit),
            remote_head_ssh=config.schedule.get('remote_head_ssh', OPTIONS.remote_head_ssh)
        )
        pipeline = Pipeline(output_dir=config.output, scheduler=scheduler)

        annotation_output_files = []
        for libconf in config.libraries.values():
            base = os.path.join(config.output, '{}_{}_{}'.format(libconf.library, libconf.disease_status, libconf.protocol))
            LOG('setting up the directory structure for', libconf.library, 'as', base)
            libconf.inputs = run_conversion(config, libconf, conversion_dir)

            # run the cluster stage
            cluster_output = mkdirp(os.path.join(base, SUBCOMMAND.CLUSTER))  # creates the clustering output dir
            args = cluster_args(config, libconf)
            args.update({'batch_id': pipeline.batch_id, 'output': cluster_output})
            args['split_only'] = SUBCOMMAND.CLUSTER in config.get('skip_stage', [])
            args['inputs'] = libconf.inputs
            LOG('clustering', '(split only)' if args['split_only'] else '', time_stamp=True)
            clustering_log = os.path.join(args['output'], 'MC_{}_{}.log'.format(libconf.library, pipeline.batch_id))
            LOG('writing:', clustering_log, time_stamp=True)
            args['log'] = clustering_log
            clustered_files = _main(cls.format_args(SUBCOMMAND.CLUSTER, args))

            # make a validation job for each cluster file
            validate_jobs = []

            if SUBCOMMAND.VALIDATE not in config.skip_stage:
                mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE))
                for task_ident in range(1, len(clustered_files) + 1):
                    mkdirp(os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{}'.format(pipeline.batch_id, task_ident)))
                args = validate_args(config, libconf)

                script_name = os.path.join(base, SUBCOMMAND.VALIDATE, 'submit.sh')
                job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
                job_options['memory_limit'] = config.schedule.validation_memory

                if libconf.protocol == PROTOCOL.TRANS:
                    job_options['memory_limit'] = config.schedule.trans_validation_memory

                if scheduler.NAME == SCHEDULER.LOCAL:
                    job_options['reference_genome'] = args['reference_genome']
                    if libconf.protocol == PROTOCOL.TRANS:
                        job_options['annotations'] = args['annotations']

                    for task_ident in range(1, len(clustered_files) + 1):
                        args['inputs'] = [os.path.join(cluster_output, '{}-{}.tab'.format(pipeline.batch_id, task_ident))]
                        args['output'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{}'.format(pipeline.batch_id, task_ident))
                        job_name = 'MV_{}_{}-{}'.format(libconf.library, pipeline.batch_id, task_ident)
                        args['log'] = os.path.join(args['output'], 'job-{name}-{job_ident}.log')
                        validate_job = LocalJob(
                            stage=SUBCOMMAND.VALIDATE,
                            output_dir=args['output'],
                            stdout=args['log'],
                            name=job_name,
                            args=cls.format_args(SUBCOMMAND.VALIDATE, args),
                            func=_main,
                            **job_options
                        )
                        pipeline.validations.append(validate_job)
                        validate_jobs.append(validate_job)
                else:
                    args['inputs'] = os.path.join(cluster_output, '{}-${}.tab'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
                    args['output'] = os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
                    aligner_path = shutil.which(args['aligner'].split(' ')[0])
                    job_class = ArrayJob if scheduler.NAME != SCHEDULER.TORQUE else TorqueArrayJob
                    validate_job = job_class(
                        stage=SUBCOMMAND.VALIDATE,
                        task_list=len(clustered_files),
                        output_dir=os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{{task_ident}}'.format(pipeline.batch_id)),
                        script=script_name,
                        name='MV_{}_{}'.format(libconf.library, pipeline.batch_id),
                        **job_options
                    )
                    pipeline.write_submission_script(SUBCOMMAND.VALIDATE, validate_job, args, aligner_path=aligner_path)
                    pipeline.validations.append(validate_job)
                    validate_jobs.append(validate_job)

            # make an annotation job for each validation/cluster job/file
            mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE))
            for task_ident in range(1, len(clustered_files) + 1):
                mkdirp(os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-{}'.format(pipeline.batch_id, task_ident)))
            args = annotate_args(config, libconf)

            script_name = os.path.join(base, SUBCOMMAND.ANNOTATE, 'submit.sh')
            job_options = {k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            job_options['memory_limit'] = config.schedule.annotation_memory

            if isinstance(scheduler, LocalScheduler):
                job_options['annotations'] = args['annotations']
                job_options['reference_genome'] = args['reference_genome']
                if args['template_metadata']:
                    job_options['template_metadata'] = args['template_metadata']
                for task_ident in range(1, len(clustered_files) + 1):
                    args['output'] = os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-{}'.format(pipeline.batch_id, task_ident))
                    # annotate 'clustered' files if the pipeline does not include the validation step
                    if SUBCOMMAND.VALIDATE not in config.skip_stage:
                        args['inputs'] = [os.path.join(base, SUBCOMMAND.VALIDATE, '{}-{}'.format(pipeline.batch_id, task_ident), _VALIDATE.PASS_FILENAME)]
                    else:
                        args['inputs'] = [os.path.join(cluster_output, '{}-{}.tab'.format(pipeline.batch_id, task_ident))]
                    job_name = 'MA_{}_{}-{}'.format(libconf.library, pipeline.batch_id, task_ident)
                    args['log'] = os.path.join(args['output'], 'job-{name}-{job_ident}.log')
                    annotate_job = LocalJob(
                        stage=SUBCOMMAND.ANNOTATE,
                        script=script_name,
                        name=job_name,
                        stdout=args['log'],
                        output_dir=args['output'],
                        args=cls.format_args(SUBCOMMAND.ANNOTATE, args),
                        func=_main,
                        **job_options
                    )
                    pipeline.annotations.append(annotate_job)
                    annotation_output_files.append(os.path.join(args['output'], _ANNOTATE.PASS_FILENAME))
                    if validate_jobs:
                        annotate_job.dependencies.append(validate_jobs[task_ident - 1])
            else:
                args['output'] = os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))
                # annotate 'clustered' files if the pipeline does not include the validation step
                if SUBCOMMAND.VALIDATE not in config.skip_stage:
                    args['inputs'] = [os.path.join(base, SUBCOMMAND.VALIDATE, '{}-${}'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT), _VALIDATE.PASS_FILENAME)]
                else:
                    args['inputs'] = [os.path.join(cluster_output, '{}-${}.tab'.format(pipeline.batch_id, scheduler.ENV_TASK_IDENT))]

                job_class = ArrayJob if scheduler.NAME != SCHEDULER.TORQUE else TorqueArrayJob
                annotate_job = job_class(
                    stage=SUBCOMMAND.ANNOTATE,
                    task_list=len(clustered_files),
                    script=script_name,
                    name='MA_{}_{}'.format(libconf.library, pipeline.batch_id),
                    output_dir=os.path.join(base, SUBCOMMAND.ANNOTATE, '{}-{{task_ident}}'.format(pipeline.batch_id)),
                    **job_options
                )
                pipeline.write_submission_script(SUBCOMMAND.ANNOTATE, annotate_job, args)
                pipeline.annotations.append(annotate_job)
                if validate_jobs:
                    annotate_job.dependencies.extend(validate_jobs)

                # add the expected output file names for input to pairing
                for taskid in range(1, len(clustered_files) + 1):
                    fname = os.path.join(args['output'], _ANNOTATE.PASS_FILENAME)
                    fname = re.sub(r'\${}'.format(scheduler.ENV_TASK_IDENT), str(taskid), fname)
                    annotation_output_files.append(fname)

        # set up the pairing job
        args = {}
        args.update(config.pairing.items())
        args['output'] = os.path.join(config.output, SUBCOMMAND.PAIR)
        args['annotations'] = config.reference.annotations
        mkdirp(args['output'])
        args['inputs'] = annotation_output_files
        job_name = 'MP_{}'.format(pipeline.batch_id)

        script_name = os.path.join(config.output, SUBCOMMAND.PAIR, 'submit.sh')

        if isinstance(scheduler, LocalScheduler):
            args['log'] = os.path.join(args['output'], 'job-{name}-{job_ident}.log')
            pipeline.pairing = LocalJob(
                stage=SUBCOMMAND.PAIR,
                script=script_name,
                output_dir=args['output'],
                stdout=args['log'],
                name=job_name,
                dependencies=pipeline.annotations,
                args=cls.format_args(SUBCOMMAND.PAIR, args),
                func=_main,
                **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            )
        else:
            pipeline.pairing = Job(
                SUBCOMMAND.PAIR,
                script=script_name,
                output_dir=args['output'],
                name=job_name,
                dependencies=pipeline.annotations,
                **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            )
            pipeline.write_submission_script(SUBCOMMAND.PAIR, pipeline.pairing, args)
        # set up the summary job
        args = summary_args(config)
        args['output'] = os.path.join(config.output, SUBCOMMAND.SUMMARY)
        mkdirp(args['output'])
        args['inputs'] = [os.path.join(config.output, SUBCOMMAND.PAIR, 'mavis_paired*.tab')]
        script_name = os.path.join(args['output'], 'submit.sh')
        job_name = 'MS_{}'.format(pipeline.batch_id)
        if isinstance(scheduler, LocalScheduler):
            args['log'] = os.path.join(args['output'], 'job-{name}-{job_ident}.log')
            pipeline.summary = LocalJob(
                stage=SUBCOMMAND.SUMMARY,
                name=job_name,
                output_dir=args['output'],
                stdout=args['log'],
                script=script_name,
                dependencies=[pipeline.pairing],
                args=cls.format_args(SUBCOMMAND.SUMMARY, args),
                func=_main,
                **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            )
        else:
            pipeline.summary = Job(
                stage=SUBCOMMAND.SUMMARY,
                name=job_name,
                output_dir=args['output'],
                script=script_name,
                dependencies=[pipeline.pairing],
                **{k: v for k, v in config.schedule.items() if k in STD_OPTIONS}
            )
            pipeline.write_submission_script(SUBCOMMAND.SUMMARY, pipeline.summary, args)
        return pipeline

    def _resubmit_job(self, job):
        """
        Given a failed job, cancel it and all of its dependencies and then resubmit them
        """
        # resubmit the job or all failed tasks for the job. Update any dependencies
        failed_tasks = set()
        try:
            for task in job.task_list:
                if task.status in self.ERROR_STATES:
                    failed_tasks.add(task.task_ident)
            if len(failed_tasks) == len(job.task_list):
                failed_tasks = []
        except AttributeError:  # non-array jobs
            pass
        # SGE cannot submit a task list that is non-consecutive so we will cancel the entire array
        if self.scheduler.NAME == SCHEDULER.SGE and len(consecutive_ranges(failed_tasks)) != 1:
            failed_tasks = set()

        if failed_tasks:  # resubmit failed tasks only and create a new job
            new_job = job.copy_with_tasks(failed_tasks)
            for task_ident in failed_tasks:
                # cancel and remove the failed task
                try:
                    self.scheduler.cancel(job, task_ident=task_ident)
                    job.remove_task(task_ident)
                except subprocess.CalledProcessError:  # ignore cancelling errors
                    pass
            self.scheduler.submit(new_job)
        else:
            # 'clean' the current job so that it is no longer 'submitted'
            self.scheduler.cancel(job)
            job.reset()
            try:
                for task in job.task_list:
                    task.status = JOB_STATUS.NOT_SUBMITTED
                    task.status_comment = ''
            except AttributeError:
                pass
            self.scheduler.submit(job)
            new_job = job

        if new_job.stage == SUBCOMMAND.VALIDATE:
            if new_job not in self.validations:
                self.validations.append(new_job)
            # cancel and resubmit annotate, pairing and summary jobs
            new_annotations = []
            for ajob in self.annotations:
                if ajob.dependencies == [job] and failed_tasks:  # only dependent on this job
                    try:
                        new_ajob = ajob.copy_with_tasks(failed_tasks)
                        new_annotations.append(new_ajob)
                        for task in failed_tasks:
                            self.scheduler.cancel(ajob, task_ident=task)
                            ajob.remove_task(task)
                        new_ajob.dependencies = [new_job]
                        self.scheduler.submit(new_ajob)
                    except AttributeError:
                        self.scheduler.cancel(ajob)
                        ajob.reset()
                        if new_job not in ajob.dependencies:
                            ajob.dependencies.append(new_job)
                        self.scheduler.submit(ajob)
                elif job in ajob.dependencies:
                    # dependent on multiple jobs
                    self.scheduler.cancel(ajob)
                    ajob.reset()
                    if new_job not in ajob.dependencies:
                        ajob.dependencies.append(new_job)
                    self.scheduler.submit(ajob)
                # ignore annotation jobs not related to the failed validation job
            self.annotations.extend(new_annotations)
        elif new_job.stage == SUBCOMMAND.ANNOTATE:
            if new_job not in self.annotations:
                self.annotations.append(new_job)

        if new_job.stage in {SUBCOMMAND.VALIDATE, SUBCOMMAND.ANNOTATE}:
            # cancel pairing
            self.scheduler.cancel(self.pairing)
            self.pairing.reset()
            self.pairing.dependencies = self.annotations[:]

        # all resubmissions result in cancelling summary
        self.scheduler.cancel(self.summary)
        self.summary.reset()

    def _job_status(self, job, submit=False, resubmit=False, log=DEVNULL):
        """
        report information regarding a particular job status
        """
        run_time = -1
        if not job.job_ident and (submit or resubmit):
            self.scheduler.submit(job)
        elif job.job_ident and resubmit and job.status in self.ERROR_STATES:
            self._resubmit_job(job)
        if job.job_ident:
            log('{} ({}) is {}'.format(job.name, job.job_ident, job.status))
        else:
            log('{} is {}'.format(job.name, job.status))
        if job.status == JOB_STATUS.COMPLETED:
            if isinstance(job, ArrayJob):
                for task in job.task_list:
                    if not os.path.exists(task.complete_stamp()):
                        log('complete stamp is expected but does not exist', indent_level=1)
                        log(task.complete_stamp(), indent_level=2)
                    else:
                        run_time = max(run_time, parse_run_time(task.complete_stamp()))
            elif not os.path.exists(job.complete_stamp()):
                with log.indent() as log:
                    log('complete stamp is expected but does not exist')
                    log(job.complete_stamp())
            else:
                run_time = max(run_time, parse_run_time(job.complete_stamp()))
            if run_time >= 0:
                if isinstance(job, ArrayJob):
                    log('{} {} COMPLETED'.format(job.tasks, 'task is' if job.tasks == 1 else 'tasks are'), indent_level=1)
                log('run time: {}'.format(run_time), indent_level=1)
        else:
            if isinstance(job, ArrayJob):
                tasks_by_status = {}
                for task in job.task_list:
                    tasks_by_status.setdefault(task.status, []).append(task)
                for status, tasks in tasks_by_status.items():
                    comments = set([t.status_comment for t in tasks if t.status_comment])
                    context = 'tasks are' if len(tasks) != 1 else 'task is'
                    LOG('{} {} {}'.format(len(tasks), context, status), indent_level=2)
                    for comment in comments:
                        LOG('comment:', comment, indent_level=3)
            elif job.status not in {JOB_STATUS.PENDING, JOB_STATUS.NOT_SUBMITTED, JOB_STATUS.SUBMITTED}:
                try:
                    content = LogFile.parse(job.logfile())
                    log('{}: {}'.format(content.status, content.message), indent_level=1)
                except FileNotFoundError:
                    log('missing log file:', job.logfile(), indent_level=1)

        return run_time

    def check_status(self, submit=False, resubmit=False, log=DEVNULL):
        """
        Check all jobs for completetion. Report any failures, etc.

        Args:
            submit (bool): submit any pending jobs
        """
        # update the information for all jobs where possible
        run_times = [[], [], [], []]
        jobs_not_complete = 0
        jobs_with_errors = 0

        for job in self.validations + self.annotations + [self.pairing, self.summary]:
            self.scheduler.update_info(job)
        log('validate', time_stamp=True)
        for job in self.validations:
            run_time = self._job_status(job, submit=submit, resubmit=resubmit, log=log.indent())
            if job.status == JOB_STATUS.COMPLETED:
                if run_time >= 0:
                    run_times[0].append(run_time)
        self.scheduler.wait()

        log('annotate', time_stamp=True)
        if not all([job.status == JOB_STATUS.COMPLETED for job in self.validations]) and self.scheduler.NAME == 'LOCAL' and (submit or resubmit):
            log('Stopping submission. Dependencies not complete', indent_level=1)
            submit = False
            resubmit = False

        for job in self.annotations:
            self._job_status(job, submit=submit, resubmit=resubmit, log=log.indent())
            if job.status == JOB_STATUS.COMPLETED:
                if run_time >= 0:
                    run_times[1].append(run_time)
        self.scheduler.wait()

        log('pairing', time_stamp=True)
        if not all([job.status == JOB_STATUS.COMPLETED for job in self.annotations]) and self.scheduler.NAME == 'LOCAL' and (submit or resubmit):
            log('Stopping submission. Dependencies not complete', indent_level=1)
            submit = False
            resubmit = False

        run_time = self._job_status(self.pairing, submit=submit, resubmit=resubmit, log=log.indent())
        if self.pairing.status == JOB_STATUS.COMPLETED:
            if run_time >= 0:
                run_times[2].append(run_time)
        self.scheduler.wait()

        log('summary', time_stamp=True)
        if self.pairing.status != JOB_STATUS.COMPLETED and self.scheduler.NAME == 'LOCAL' and (submit or resubmit):
            log('Stopping submission. Dependencies not complete', indent_level=1)
            submit = False
            resubmit = False

        run_time = self._job_status(self.summary, submit=submit, resubmit=resubmit, log=log.indent())
        if self.summary.status == JOB_STATUS.COMPLETED:
            if run_time >= 0:
                run_times[3].append(run_time)
        self.scheduler.wait()

        for job in self.validations + self.annotations + [self.pairing, self.summary]:
            if submit or resubmit and job.status != JOB_STATUS.COMPLETED:
                self.scheduler.update_info(job)
            if job.status in self.ERROR_STATES:
                jobs_with_errors += 1
            elif job.status != JOB_STATUS.COMPLETED:
                jobs_not_complete += 1

        if jobs_not_complete + jobs_with_errors == 0:
            if all([r for r in run_times]):
                log('parallel run time:', sum([max(r) for r in run_times]))
            return EXIT_OK
        elif not jobs_with_errors:
            return EXIT_INCOMPLETE
        else:
            return EXIT_ERROR

    @classmethod
    def read_build_file(cls, filepath):
        """
        read the configuration file which stored the build information concerning jobs and dependencies

        Args:
            filepath (str): path to the input config file
        """
        from ..main import main as _main

        if not os.path.exists(filepath):
            raise FileNotFoundError('File does not exist: {}'.format(filepath))
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser.read(filepath)
        cast = {'None': None, 'False': False, 'True': True}

        pipeline = cls(
            output_dir=parser['general']['output_dir'],
            scheduler=SCHEDULERS_BY_NAME[parser['general']['scheduler']](
                concurrency_limit=parser['general']['concurrency_limit'] if 'concurrency_limit' in parser['general'] else OPTIONS.concurrency_limit,
                remote_head_ssh=parser['general']['remote_head_ssh'] if 'remote_head_ssh' in parser['general'] else OPTIONS.remote_head_ssh
            ),
            batch_id=parser['general']['batch_id']
        )

        jobs = {}
        for sec in parser.sections():
            if sec != 'general':
                section = {}
                for attr, value in parser[sec].items():
                    if attr in ['dependencies', 'inputs', 'outputs', 'args', 'task_list'] and value:
                        section[attr] = [s.strip() for s in re.split(r'\n', value)]
                    elif value == 'None':
                        section[attr] = None
                    elif value in cast:
                        value = cast[value]
                    else:
                        section[attr] = value
                if pipeline.scheduler.NAME == SCHEDULER.LOCAL:
                    jobs[sec] = LocalJob(func=_main, **section)
                elif 'task_list' in section:
                    if pipeline.scheduler.NAME == SCHEDULER.TORQUE:
                        jobs[sec] = TorqueArrayJob(**section)
                    else:
                        jobs[sec] = ArrayJob(**section)
                else:
                    jobs[sec] = Job(**section)

        for job in jobs.values():
            for i, prior_job_name in enumerate(job.dependencies):
                job.dependencies[i] = jobs[prior_job_name]

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

        Args:
            filename (str): path to the output config file
        """
        parser = ConfigParser(interpolation=ExtendedInterpolation())
        parser['general'] = {
            'batch_id': self.batch_id,
            'output_dir': self.output_dir,
            'scheduler': self.scheduler.NAME,
            'remote_head_ssh': self.scheduler.remote_head_ssh,
            'concurrency_limit': str(self.scheduler.concurrency_limit)
        }

        for job in [self.summary, self.pairing] + self.validations + self.annotations:
            parser[job.display_name] = {k: re.sub(r'\$', '$$', v) for k, v in job.flatten().items()}

        with open(filename, 'w') as configfile:
            parser.write(configfile)
