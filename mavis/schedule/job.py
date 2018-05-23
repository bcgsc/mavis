from datetime import timedelta
import subprocess
import re

from ..constants import SUBCOMMAND, MavisNamespace
from .constants import JOB_STATUS, OPTIONS, STD_OPTIONS, SCHEDULER


class Scheduler:
    """
    Class responsible for methods interacting with the scheduler
    """
    ARRAY_DEPENDENCY = None
    JOB_DEPENDENCY = None
    ENV_TASK_IDENT = 'TASK_IDENT'

    @classmethod
    def submit(cls, job, task_ident=None, cascade=False):
        raise NotImplementedError('abstract method')

    @classmethod
    def status(cls, job):
        raise NotImplementedError('abstract method')

    @classmethod
    def cancel(cls, job):
        raise NotImplementedError('abstract method')

    @classmethod
    def acct(cls, job):
        raise NotImplementedError('abstract method')

    @classmethod
    def format_dependencies(cls, job, task_ident=None, cascade=False):
        """
        returns a string representing the dependency argument
        """
        if not cls.ARRAY_DEPENDENCY or not cls.JOB_DEPENDENCY:
            raise NotImplementedError('Abstract class. Implementing class must define class attributes ARRAY_DEPENDENCY and JOB_DEPENDENCY')

        if len(job.dependencies) == 1 and isinstance(job, ArrayJob) and isinstance(job.dependencies[0], ArrayJob) and task_ident is None:
            dependency = job.dependencies[0]
            if dependency.tasks != job.tasks:
                raise ValueError('An array job must be dependent only on another single array job with the same number of tasks', job, dependency)
            if not cascade and not dependency.job_ident:
                raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)
            elif not dependency.job_ident:
                cls.submit(dependency, task_ident, cascade=cascade)
            if task_ident is not None:
                return cls.ARRAY_DEPENDENCY.format('{}_{}'.format(dependency.job_ident, dependency.task_ident))
            return cls.ARRAY_DEPENDENCY.format(dependency.job_ident)
        for djob in job.dependencies:
            if not cascade and not djob.job_ident:
                raise ValueError('The dependencies must be submitted before the dependent job', job, djob)
            elif not djob.job_ident:
                cls.submit(djob, cascade=cascade)
        return cls.JOB_DEPENDENCY.format(','.join([d.job_ident for d in job.dependencies]))


class Job:

    def __init__(
        self,
        stage,
        job_ident=None,
        name=None,
        dependencies=None,
        stdout='job-{name}-{job_ident}.log',
        script=None,
        status=JOB_STATUS.UNKWOWN,
        **options
    ):
        """
        Args
            stage (str): the mavis pipleine stage this job belongs to
            job_ident (int): the job number/id according to the scheduler being used
            name (str): the job name according to the scheduler being used
            dependencies (list of Job): list of jobs which must complete for this job to run
            stdout (str): basename of the file to write std output to
            script (str): path to the script which contains the commands for the job
            options (**dict): override default options specified by OPTIONS
        """
        self.stage = SUBCOMMAND.enforce(stage)
        self.job_ident = job_ident
        self.name = name
        self.dependencies = dependencies if dependencies else []
        self.stdout = stdout
        self.script = script
        self.status = status

        self.status = JOB_STATUS.UNKWOWN

        # inputs to the function call should override the default values
        for option, value in [(o, OPTIONS[o]) for o in STD_OPTIONS]:
            setattr(self, option, options.get(option, value))

        # check that nothing weird was passed in the kwargs
        for option in options:
            if option not in STD_OPTIONS:
                raise AttributeError('unexpected attribute: {}'.format(option))

    def flatten(self):
        result = {}
        for attr, value in self.__dict__.items():
            if attr == 'dependencies':
                value = '\n'.join(sorted([str(j.name) for j in self.dependencies]))
            result[attr] = str(value)
        return result


class ArrayJob(Job):

    def __init__(self, stage, tasks, concurrency_limit=OPTIONS.concurrency_limit, **kwargs):
        Job.__init__(self, stage, **kwargs)
        self.concurrency_limit = concurrency_limit
        self.tasks = tasks


class JobProfile:
    def __init__(
        self,
        sumbitted_at=None,
        started_at=None,
        stopped_at=None,
        elapsed_time=None,
        max_memory_used=None
    ):
        self.sumbitted_at = sumbitted_at
        self.started_at = started_at
        self.stopped_at = stopped_at
        self.elapsed_time = elapsed_time
        self.max_memory_used = max_memory_used


class SlurmScheduler(Scheduler):
    NAME = SCHEDULER.SLURM
    ARRAY_DEPENDENCY = '--dependency=aftercorr:{}'
    JOB_DEPENDENCY = '--dependency=afterok:{}'
    ENV_TASK_IDENT = 'SLURM_ARRAY_TASK_ID'

    @staticmethod
    def format_stdout(job):
        if isinstance(job, ArrayJob):
            name = job.stdout.format(name='%x', job_ident='%A_%a')
        else:
            name = job.stdout.format(name='%x', job_ident='%j')
        return name

    @classmethod
    def submit(cls, job, task_ident=None, cascade=False):
        """
        runs a subprocess sbatch command

        Args
            job (Job): the job to be submitted
            task_ident (int): submit only a particular task of the current job
            cascade (bool): submit dependencies if they have not yet been submitted
        """
        command = ['sbatch']
        if job.job_ident:
            raise ValueError('Job has already been submitted and has the job number', job.job_ident)
        if job.queue:
            command.append('--partition={}'.format(job.queue))
        if job.memory_limit:
            command.extend(['--mem', str(job.memory_limit)])
        if job.time_limit:
            command.extend(['-t', str(timedelta(seconds=job.time_limit))])
        if job.import_env:
            command.append('--export=ALL')
        if job.dependencies:
            command.append(cls.format_dependencies(job, task_ident=task_ident, cascade=cascade))
        if job.name:
            command.extend(['-J', job.name])
        if job.stdout:
            command.extend(['-o', cls.format_stdout(job)])
        if job.mail_type:
            command.append('--mail-type={}'.format(job.mail_type))
        if job.mail_user:
            command.append('--mail-user={}'.format(job.mail_user))
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            concurrency_limit = '' if job.concurrency_limit is None else '%{}'.format(job.concurrency_limit)
            if task_ident is None:  # default to all
                command.append('--array=1-{}{}'.format(job.tasks, concurrency_limit))
            else:
                command.append('--array={}{}'.format(task_ident, concurrency_limit))


        command.append(job.script)
        content = subprocess.check_output(command).decode('utf8').strip()

        match = re.match(r'^submitted batch job (\d+)$', content, re.IGNORECASE)
        if not match:
            raise NotImplementedError('Error in retrieving the submitted job number. Did not match the expected pattern', content)
        job.job_ident = match.group(1)
        job.status = JOB_STATUS.SUBMITTED

    @classmethod
    def status(cls, job, task_ident=None):
        """
        runs a subprocess scontrol command to get job details
        """
        command = ['scontrol', 'show']
        if task_ident is not None:
            command.append('job={}_{}'.format(job.job_ident, task_ident))
        else:
            command.append('job={}'.format(job.job_ident))
        content = subprocess.check_output(command).decode('utf8').strip()
        # split by empty lines for parsing multiple job results
        status_hist = {}
        for job_content in re.split(r'\s*\n\n\s*', content):
            parsed = {}
            for line in re.split(r'\s+', job_content):
                attr, value = line.split('=', 1)
                parsed[attr] = value
            job_status = parsed['JobState']
            if job_status not in JOB_STATUS:
                raise NotImplementedError('Error parsing the unknown job status value', job_status)
            status_hist[job_status] = status_hist.get(job_status, 0) + 1
        return status_hist

    @classmethod
    def cancel(cls, job):
        command = ['scancel', job.job_ident]
        subprocess.check_output(command)
        job.job_ident = None
        job.status = JOB_STATUS.CANCELLED


def parse_fixed_width_table(string):
    string = string.strip()
    lines = [l for l in string.split('\n') if l.strip()]
    header = []
    for char in lines.pop(0):
        if not header or not re.match(r'\s', char):
            header.append(char)
        else:
            header[-1] = header[-1] + char
    column_sizes = [len(col) for col in header]
    header = [col.strip() for col in header]
    rows = []

    for line in lines:
        if re.match(r'^[\-]+$', line):
            continue  # ignore dashed separators
        row = {}
        pos = 0
        for col, size in zip(header, column_sizes):
            row[col] = line[pos:pos + size]
            pos += size
        rows.append(row)
    return rows


class SgeScheduler(Scheduler):
    NAME = SCHEDULER.SGE
    ARRAY_DEPENDENCY = '-hold_jid {}'
    JOB_DEPENDENCY = '-hold_jid_ad {}'
    ENV_TASK_IDENT = 'TASK_ID'
    ENV_JOB_IDENT = 'JOB_ID'
    ENV_JOB_NAME = 'JOB_NAME'

    @classmethod
    def format_stdout(cls, job):
        if isinstance(job, ArrayJob):
            name = job.stdout.format(
                name='\${}'.format(cls.ENV_JOB_NAME),
                job_ident='\${}-\${}'.format(
                    cls.ENV_JOB_IDENT,
                    cls.ENV_TASK_IDENT))
        else:
            name = job.stdout.format(
                name='\${}'.format(cls.ENV_JOB_NAME),
                job_ident='\${}'.format(cls.ENV_JOB_IDENT))
        return name

    @classmethod
    def submit(cls, job, task_ident=None, cascade=False):
        """
        runs a subprocess sbatch command

        Args
            job (Job): the job to be submitted
            task_ident (int): submit only a particular task of the current job
            cascade (bool): submit dependencies if they have not yet been submitted
        """
        command = ['qsub', '-j', 'y']  # always join output
        if job.job_ident:
            raise ValueError('Job has already been submitted and has the job number', job.job_ident)
        if job.queue:
            command.append('-q {}'.format(job.queue))
        if job.memory_limit:
            command.extend([
                '-l',
                'mem_free={0}M,mem_token={0}M,h_vmem={}M'.format(job.memory_limit)
            ])
        if job.time_limit:
            command.extend([
                '-l',
                'h_rt={}'.format(str(timedelta(seconds=job.time_limit)))])
        if job.import_env:
            command.append('-V')
        if job.dependencies:
            command.append(cls.format_dependencies(job, task_ident=task_ident, cascade=cascade))
        if job.name:
            command.extend(['-J', job.name])
        if job.stdout:
            command.extend(['-o', cls.format_stdout(job)])
        if job.mail_type:
            command.extend(['-m', job.mail_type])
        if job.mail_user:
            command.extend(['-M', job.mail_user])
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            if task_ident is None:  # default to all
                command.extend(['-t', '1-{}'.format(job.tasks)])
            else:
                command.append(['-t', str(task_ident)])

        command.append(job.script)
        content = subprocess.check_output(command).decode('utf8').strip()

        match = re.match(r'^submitted batch job (\d+)$', content, re.IGNORECASE)
        if not match:
            raise NotImplementedError('Error in retrieving the submitted job number. Did not match the expected pattern', content)
        job.job_ident = match.group(1)
        job.status = JOB_STATUS.SUBMITTED

    @classmethod
    def status(cls, job, task_ident=None):
        """
        runs a subprocess scontrol command to get job details
        """
        command = ['qstat']
        if job.queue:
            command.extend(['-q', job.queue])
        content = subprocess.check_output(command).decode('utf8').strip()

        for job in parse_fixed_width_table(content):
            job_status = row['state']
            if job_status not in JOB_STATUS:
                raise NotImplementedError('Error parsing the unknown job status value', job_status)
            status_hist[job_status] = status_hist.get(job_status, 0) + 1
        return status_hist

    @classmethod
    def cancel(cls, job):
        command = ['scancel', job.job_ident]
        subprocess.check_output(command)
        job.job_ident = None
        job.status = JOB_STATUS.CANCELLED

class PbsScheduler(Scheduler):
    NAME = SCHEDULER.PBS

class LocalScheduler(Scheduler):
    NAME = SCHEDULER.NONE


