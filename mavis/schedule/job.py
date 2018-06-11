from copy import copy as _copy
import os
import re
import time

from ..constants import SUBCOMMAND, MavisNamespace
from .constants import JOB_STATUS, OPTIONS, STD_OPTIONS


class LogFile:
    """
    stores information about the log status
    """
    STATUS = MavisNamespace('EMPTY', 'CRASH', 'INCOMPLETE', 'COMPLETE')
    """:class:`~mavis.constants.MavisNamespace`: The status of the job based on parsing of the logfile"""

    def __init__(self, filename, status, message=None):
        """
        Args:
            filename (str): path to the logfile
            status (LogFile.STATUS): the status of the logfile
            message (str): the message parsed from the logfile. Generally this is an error from the log
        """
        self.filename = filename
        self.status = self.STATUS.enforce(status)
        self.message = message.strip() if message is not None else None

    @classmethod
    def parse(cls, filename):
        """
        given a file parse to see if it looks like a complete log file (contains run time),
        was truncated, or reported an error
        """
        if not os.path.isfile(filename):
            raise FileNotFoundError('Log file does not exist', filename)
        log = None
        with open(filename, 'r') as fh:
            lines = [l.strip() for l in fh.readlines() if l.strip()]
            for line in lines[::-1]:
                line = line.strip().lower()
                if line and line[0] != '\x1b':  # ignore lines starting with terminal control characters
                    if re.search(r'(\b|^)((\S+)?error|fault|fatal|aborted|core dumped|killed|died|command not found)(\b|$)', line):
                        log = LogFile(filename, cls.STATUS.CRASH, line)
                    elif re.match(r'^\s*run time \(s\): (\d+)\s*$', line):
                        log = LogFile(filename, cls.STATUS.COMPLETE)
                    else:
                        log = LogFile(filename, cls.STATUS.INCOMPLETE, line)
                    return log
            return LogFile(filename, cls.STATUS.EMPTY)


class Job:

    def __init__(
        self,
        stage,
        output_dir,
        stdout=None,
        job_ident=None,
        name=None,
        dependencies=None,
        script=None,
        created_at=None,
        status=JOB_STATUS.NOT_SUBMITTED,
        status_comment='',
        **options
    ):
        """
        Args:
            stage (str): the mavis pipleine stage this job belongs to
            job_ident (int): the job number/id according to the scheduler being used
            output_dir (str): path to the output directory where logs/stamps for this job will be written
            name (str): the job name according to the scheduler being used
            dependencies (list of Job): list of jobs which must complete for this job to run
            stdout (str): basename of the file to write std output to
            script (str): path to the script which contains the commands for the job
            created_at (int): the time stamp for when the job was created (created != submitted)
            status (~mavis.schedule.constants.JOB_STATUS): The current (since last checked) status of the job
            status_comment (str): the comment which describes the status, generally this is used for reporting errors from the log file or failed dependencies (SLURM)
            options (**dict): override default options specified by OPTIONS
        """
        self.stage = SUBCOMMAND.enforce(stage)
        self.job_ident = job_ident
        self.name = name
        self.dependencies = dependencies if dependencies else []
        self.script = script
        self.status = JOB_STATUS.enforce(status)
        self.output_dir = output_dir
        self.stdout = os.path.join(output_dir, 'job-{name}-{job_ident}.log') if not stdout else stdout

        self.created_at = int(created_at if created_at else time.time())
        self.status = status
        self.status_comment = status_comment

        # inputs to the function call should override the default values
        for option, value in [(o, OPTIONS[o]) for o in STD_OPTIONS]:
            setattr(self, option, OPTIONS.type(option)(options.get(option, value)))

        # check that nothing weird was passed in the kwargs
        for option in options:
            if option not in STD_OPTIONS:
                raise AttributeError('unexpected attribute: {}'.format(option))

    @property
    def display_name(self):
        """
        Used for identifying this job in an ini config file
        """
        display_name = self.name if self.job_ident is None else '{}_{}'.format(self.name, self.job_ident)
        display_name = re.sub(r'[\[\]#;]', '_', display_name)
        return display_name

    def flatten(self):
        result = {}
        for attr, value in self.__dict__.items():
            if attr == 'dependencies':
                value = [j.display_name for j in value]
            try:
                if not isinstance(value, str):
                    value = '\n'.join([str(v) for v in value])
            except TypeError:
                pass
            result[attr] = str(value)
        return result

    def logfile(self):
        """
        returns the path to the logfile with job name and job id substituted into the stdout pattern
        """
        return self.stdout.format(name=self.name, job_ident=self.job_ident)

    def complete_stamp(self):
        """
        returns the path to the expected complete stamp
        """
        return os.path.join(self.output_dir, 'MAVIS-{job_ident}.COMPLETE').format(job_ident=self.job_ident, name=self.name)

    def reset(self):
        self.status = JOB_STATUS.NOT_SUBMITTED
        self.status_comment = ''
        self.job_ident = None


class ArrayJob(Job):
    """
    Class for dealing with array jobs. Jobs with many tasks
    """

    def __init__(self, stage, task_list, **kwargs):
        """
        Args:
            task_list (:class:`list` or :class:`int`): the ids of tasks in the job array
        """
        Job.__init__(self, stage, **kwargs)
        self.stdout = os.path.join(self.output_dir, 'job-{name}-{job_ident}-{task_ident}.log') if 'stdout' not in kwargs else kwargs['stdout']

        if isinstance(task_list, int):
            task_list = list(range(1, task_list + 1))
        self.task_list = [Task(self, n) for n in task_list]

    @property
    def tasks(self):
        return len(self.task_list)

    def get_task(self, task_ident):
        """
        returns a task by task id
        """
        task_ident = int(task_ident)
        for task in self.task_list:
            if task.task_ident == task_ident:
                return task
        raise KeyError('task id not found', task_ident, self.task_list)

    def has_task(self, task_ident):
        task_ident = int(task_ident)
        for task in self.task_list:
            if task.task_ident == task_ident:
                return True
        return False

    def remove_task(self, task_ident):
        self.task_list = [task for task in self.task_list if task.task_ident != task_ident]

    def logfile(self, task_ident):
        return self.stdout.format(name=self.name, job_ident=self.job_ident, task_ident=task_ident)

    def complete_stamp(self, task_ident):
        """
        returns the path to the expected complete stamp
        """
        return os.path.join(self.output_dir, 'MAVIS-{job_ident}.COMPLETE').format(job_ident=self.job_ident, name=self.name, task_ident=task_ident)

    def flatten(self):
        result = {k: v for k, v in Job.flatten(self).items() if k != 'task_list'}
        result['task_list'] = '\n'.join([str(t.task_ident) for t in self.task_list])
        return result

    def copy_with_tasks(self, task_list):
        copy = _copy(self)
        copy.task_list = [Task(self, n) for n in task_list]
        copy.dependencies = []
        copy.reset()
        return copy

    def reset(self):
        Job.reset(self)
        for task in self.task_list:
            task.reset()

    def __repr__(self):
        return '{}(job_ident={}, name={}, stage={}, status={})'.format(self.__class__.__name__, self.job_ident, self.name, self.stage, self.status)


class TorqueArrayJob(ArrayJob):

    def complete_stamp(self, task_ident):
        # example: MAVIS-136[1].torque01.bcgsc.ca.COMPLETE
        job_ident = re.sub(r'\[\]', '[{}]'.format(task_ident), self.job_ident)
        return os.path.join(self.output_dir, 'MAVIS-{job_ident}.COMPLETE').format(job_ident=job_ident, name=self.name, task_ident=task_ident)

    def logfile(self, task_ident):
        # example: job-MV_mock-A47933_batch-B9PE6YAtnHu4cHA2GrsEzX-1-136[1].torque01.bcgsc.ca-1.log-1
        name = '{}-{}'.format(self.name, task_ident)
        job_ident = re.sub(r'\[\]', '[{}]'.format(task_ident), self.job_ident)
        log = self.stdout.format(name=name, job_ident=job_ident, task_ident=task_ident)
        return '{}-{}'.format(log, task_ident)


class Task:

    def __init__(self, array_job, task_ident):
        self.array_job = array_job
        self.task_ident = int(task_ident)
        self.status = JOB_STATUS.NOT_SUBMITTED
        self.status_comment = ''

    def logfile(self):
        return self.array_job.logfile(self.task_ident)

    def complete_stamp(self):
        return self.array_job.complete_stamp(self.task_ident)

    def reset(self):
        self.status = JOB_STATUS.NOT_SUBMITTED
        self.status_comment = ''
