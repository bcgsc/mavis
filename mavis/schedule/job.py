from datetime import timedelta, datetime
import os
import subprocess
import re
import time

from ..constants import SUBCOMMAND, MavisNamespace
from .constants import JOB_STATUS, OPTIONS, STD_OPTIONS, SCHEDULER, cumulative_job_state, MAIL_TYPE


class LogFile:
    """
    stores information about the log status
    """
    STATUS = MavisNamespace('EMPTY', 'CRASH', 'INCOMPLETE', 'COMPLETE')

    def __init__(self, filename, status, message=None):
        self.filename = filename
        self.status = status
        self.message = message

    @classmethod
    def parse(cls, filename):
        """
        given a file parse to see if it looks like a complete log file (contains run time),
        was truncated, or reported an error
        """
        if not os.path.exists(filename):
            raise FileNotFoundError('Log file does not exist', filename)
        log = None
        with open(filename, 'r') as fh:
            lines = [l.strip() for l in fh.readlines() if l.strip()]
            if not lines:
                log = LogFile(filename, cls.STATUS.EMPTY)
            else:
                non_empty_line = lines[-1].lower()
                if re.search(r'(\b|^)((\S+)?error|fault|fatal|aborted|core dumped|killed|died|command not found)(\b|$)', non_empty_line):
                    log = LogFile(filename, cls.STATUS.CRASH, non_empty_line.strip())
                elif any([re.match(r'^\s*run time \(s\): (\d+)\s*$', line) for line in lines[-10:]]):
                    log = LogFile(filename, cls.STATUS.COMPLETE)
                else:
                    log = LogFile(filename, cls.STATUS.INCOMPLETE, lines[-1].strip())
        return log


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
        self.script = script
        self.status = status
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

    def flatten(self):
        result = {}
        for attr, value in self.__dict__.items():
            if attr == 'dependencies':
                value = [j.name for j in value]
            try:
                if not isinstance(value, str):
                    value = '\n'.join([str(v) for v in value])
            except TypeError:
                pass
            result[attr] = str(value)
        return result

    def logfile(self):
        return self.stdout.format(name=self.name, job_ident=self.job_ident)

    def complete_stamp(self):
        return os.path.join(self.output_dir, 'MAVIS.COMPLETE')


class ArrayJob(Job):

    def __init__(self, stage, tasks, concurrency_limit=OPTIONS.concurrency_limit, **kwargs):
        Job.__init__(self, stage, **kwargs)
        self.stdout = os.path.join(self.output_dir, 'job-{name}-{job_ident}-{task_ident}.log') if 'stdout' not in kwargs else kwargs['stdout']
        self.concurrency_limit = concurrency_limit
        self.tasks = int(tasks)
        self.task_list = [Task(self, n) for n in range(1, self.tasks + 1)]

    def logfile(self, task_ident):
        return self.stdout.format(name=self.name, job_ident=self.job_ident, task_ident=task_ident)

    def complete_stamp(self, task_ident):
        return Job.complete_stamp(self).format(task_ident=task_ident)

    def flatten(self):
        result = {k:v for k, v in Job.flatten(self).items() if k != 'task_list'}
        return result


class Task:

    def __init__(self, array_job, task_id):
        self.array_job = array_job
        self.task_id = task_id
        self.status = JOB_STATUS.NOT_SUBMITTED
        self.status_comment = ''

    def logfile(self):
        return self.array_job.logfile(self.task_id)

    def complete_stamp(self):
        return self.array_job.complete_stamp(self.task_id)
