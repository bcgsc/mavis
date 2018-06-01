from datetime import timedelta
import subprocess
import re
import logging
import socket

from ..util import LOG
from ..config import NullableType

from .job import ArrayJob
from .constants import SCHEDULER, JOB_STATUS, cumulative_job_state, MAIL_TYPE


class Scheduler:
    """
    Class responsible for methods interacting with the scheduler
    """
    ARRAY_DEPENDENCY = None
    """:class:`str`: string which takes format arguments and is used to add array dependency arguments for job submission"""
    JOB_DEPENDENCY = None
    """:class:`str`: string which takes format arguments and is used to add job dependency arguments for job submission"""
    ENV_TASK_IDENT = '{TASK_IDENT}'
    """:class:`str`: the expected pattern of environment variables which store the task id"""
    ENV_JOB_IDENT = '{JOB_IDENT}'
    """:class:`str`: the expected pattern of environment variables which store the job id"""
    DEPENDENCY_DELIM = ':'
    """:class:`str`: the delimiter to use between jobs of the same depenency type"""
    HEADER_PREFIX = '#'

    def __init__(self, concurrency_limit=None, remote_head_ssh='', remote_head_name=''):
        """
        Args:
            concurrency_limit (int): Size of t, mock.Mock(), the maximum allowed concurrent processes. Defaults to one less than the total number available
        """
        self.concurrency_limit = NullableType(int)(concurrency_limit)
        self.remote_head_ssh = remote_head_ssh
        self.remote_head_name = remote_head_name

    def command(self, command, shell=False):
        """
        Wrapper to deal with subprocess commands. If configured and not on the head node currently, will send the command through ssh

        Args:
            command (list or str): the command can be a list or a string and is passed to the subprocess to be run

        Returns:
            str: the content returns from stdout of the subprocess
        """
        if self.remote_head_ssh and self.remote_head_name != socket.gethostname():
            # ssh to remote head and run the command there
            if not isinstance(command, str):
                command = ' '.join(command)
            return subprocess.check_output(['ssh', str(self.remote_head_ssh), command]).decode('utf8').strip()
        return subprocess.check_output(command, shell=shell).decode('utf8').strip()

    def wait(self):
        pass

    def submit(self, job, resubmit=False):
        """
        submit a job to the scheduler
        """
        raise NotImplementedError('abstract method')

    def update_info(self, job):
        """
        update the information about the job from the scheduler
        """
        raise NotImplementedError('abstract method')

    def cancel(self, job, task_ident=None):
        raise NotImplementedError('abstract method')

    def format_dependencies(self, job):
        """
        returns a string representing the dependency argument
        """
        raise NotImplementedError('abstract method')


class SlurmScheduler(Scheduler):
    """
    Class for formatting commands to match a SLURM scheduler system
    SLURM docs can be found here https://slurm.schedmd.com
    """
    NAME = SCHEDULER.SLURM
    """:attr:`~mavis.schedule.constants.SCHEDULER`: the type of scheduler"""

    ENV_TASK_IDENT = 'SLURM_ARRAY_TASK_ID'
    ENV_JOB_IDENT = 'SLURM_JOB_ID'
    ENV_ARRAY_IDENT = 'SLURM_ARRAY_JOB_ID'

    def submit(self, job):
        """
        runs a subprocess sbatch command

        Args:
            job (Job): the job to be submitted
            task_ident (int): submit only a particular task of the current job
            resubmit (bool): indicates if this is a resubmission (will cancel and resubmit dependent jobs and submit individual tasks where appropriate)
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
            command.append(self.format_dependencies(job))
        if job.name:
            command.extend(['-J', job.name])
        if job.stdout:
            command.extend(['-o', job.stdout.format(
                name='%x',
                job_ident='%A' if isinstance(job, ArrayJob) else '%j',
                task_ident='%a'
            )])
        if job.mail_type and job.mail_user:
            command.append('--mail-type={}'.format(job.mail_type))
            command.append('--mail-user={}'.format(job.mail_user))
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            concurrency_limit = '' if self.concurrency_limit is None else '%{}'.format(self.concurrency_limit)
            command.append('--array={}{}'.format(','.join([str(task.task_ident) for task in job.task_list]), concurrency_limit))

        command.append(job.script)
        LOG('submitting', job.name)
        content = self.command(command)

        match = re.match(r'^submitted batch job (\d+)$', content, re.IGNORECASE)
        if not match:
            raise NotImplementedError('Error in retrieving the submitted job number. Did not match the expected pattern', content)
        job.job_ident = match.group(1)
        job.status = JOB_STATUS.SUBMITTED

        try:
            for task in job.task_list:
                task.status = job.status
                task.status_comment = job.status_comment
        except AttributeError:
            pass

    @classmethod
    def parse_sacct(cls, content):
        """
        parses content returned from the sacct command

        Args:
            content (str): the content returned from the sacct command
        """
        lines = content.strip().split('\n')
        header = lines[0].split('|')
        rows = []
        for line in lines[1:]:
            row = {col: val for col, val in zip(header, line.split('|'))}
            rows.append(row)
        # now combine the .batch split jobs
        results = {}
        for row in rows:
            jobid = re.sub(r'\.batch$', '', row['JobID'])
            if row['JobName'] != 'batch':
                results[jobid] = row
        for row in rows:
            jobid = re.sub(r'\.batch$', '', row['JobID'])
            if row['JobName'] == 'batch' and jobid in results:
                curr = results[jobid]
                for col, val in row.items():
                    if not curr[col]:
                        curr[col] = val
        rows = []
        for row in results.values():
            row['State'] = row['State'].split(' ')[0]
            if re.match(r'^\d+_\d+$', row['JobID']):
                job_ident, task_ident = row['JobID'].rsplit('_', 1)
                task_ident = int(task_ident)
            else:
                job_ident = row['JobID']
                task_ident = None
            rows.append({
                'job_ident': job_ident,
                'task_ident': task_ident,
                'name': row['JobName'],
                'status': row['State'],
                'status_comment': ''
            })

        return rows

    @classmethod
    def parse_scontrol_show(cls, content):
        """
        parse the content from the command: scontrol show job <JOBID>
        """
        rows = []
        for job_content in re.split(r'\n\s*\n', content):
            job_content = job_content.strip()
            if not job_content:  # ignore empty
                continue
            row = {}
            for pair in re.split(r'\s+', job_content):
                if '=' not in pair:
                    continue
                col, val = pair.split('=', 1)
                row[col] = val
            try:
                task_ident = int(row.get('ArrayTaskId', ''))
            except ValueError:
                task_ident = None
            rows.append({
                'job_ident': row['JobId'],
                'status': row['JobState'],
                'name': row['JobName'],
                'status_comment': row['Reason'] if row['Reason'].lower() != 'none' else '',
                'task_ident': task_ident
            })
        return rows

    def update_info(self, job):
        """
        Pull job information about status etc from the scheduler. Updates the input job
        """
        if not job.job_ident:
            return
        command = ['sacct', '-j', job.job_ident, '--long', '--parsable2']
        content = self.command(command)
        rows = self.parse_sacct(content)
        updated = False

        for row in rows:
            if row['job_ident'] == job.job_ident:
                if row['task_ident'] is not None:
                    if job.has_task(row['task_ident']):
                        task = job.get_task(row['task_ident'])
                        task.status = row['status']
                        task.status_comment = row['status_comment']
                else:
                    job.status = row['status']
                    job.status_comment = row['status_comment']
                    updated = True

        try:
            if not updated:
                job.status = cumulative_job_state([t.status for t in job.task_list])
        except AttributeError:
            pass

    def cancel(self, job, task_ident=None):
        """
        cancel a job
        """
        if not job.job_ident:
            return
        if task_ident is not None:
            self.command(['scancel', '{}_{}'.format(job.job_ident, task_ident)])
            job.get_task(task_ident).status = JOB_STATUS.CANCELLED
            LOG('cancelled task', job.name, job.job_ident, task_ident)
        else:
            self.command(['scancel', job.job_ident])
            job.status = JOB_STATUS.CANCELLED
            LOG('cancelled job', job.name, job.job_ident)

            try:
                for task in job.task_list:
                    task.status = JOB_STATUS.CANCELLED
            except AttributeError:
                pass

    def format_dependencies(self, job):
        """
        returns a string representing the dependency argument
        """
        try:
            if len(job.dependencies) == 1 and job.tasks == job.dependencies[0].tasks:
            # array job dependent on only another array job with the same number of tasks
                dependency = job.dependencies[0]
                if not dependency.job_ident:
                    raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)
                return '--dependency=aftercorr:{}'.format(dependency.job_ident)
        except AttributeError:
            pass

        dep_jobs = []
        for dependency in job.dependencies:
            if not dependency.job_ident:
                raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)
            try:
                for task in dependency.task_list:
                    dep_jobs.append('{}_{}'.format(dependency.job_ident, task.task_ident))
            except AttributeError:
                dep_jobs.append(str(dependency.job_ident))

        return '--dependency=afterok:{}'.format(':'.join(dep_jobs))


class SgeScheduler(Scheduler):
    """
    Class for managing interactions with the SGE scheduler
    """
    NAME = SCHEDULER.SGE
    """:attr:`~mavis.schedule.constants.SCHEDULER`: the type of scheduler"""
    ENV_TASK_IDENT = 'SGE_TASK_ID'
    ENV_JOB_IDENT = 'JOB_ID'
    """:class:`str`: expected pattern for environment variables which store the job id"""
    ENV_ARRAY_IDENT = ENV_JOB_IDENT
    ENV_JOB_NAME = 'JOB_NAME'
    """:class:`str`: expected pattern for environment variables which store the job name"""
    DEPENDENCY_DELIM = ','
    HEADER_PREFIX = '#$'

    STATE_MAPPING = {
        'q': JOB_STATUS.PENDING,
        'h': JOB_STATUS.PENDING,
        'R': JOB_STATUS.RUNNING,
        'r': JOB_STATUS.RUNNING,
        'd': JOB_STATUS.CANCELLED,
        's': JOB_STATUS.ERROR,
        'w': JOB_STATUS.PENDING,
        'E': JOB_STATUS.ERROR,
        'T': JOB_STATUS.ERROR,
        't': JOB_STATUS.RUNNING
    }
    """:class:`dict`: mapping from SGE job states to their MAVIS JOB_STATUS equivalent"""
    MAIL_TYPE_MAPPING = {
        MAIL_TYPE.BEGIN: 'b',
        MAIL_TYPE.NONE: 'n',
        MAIL_TYPE.FAIL: 'as',
        MAIL_TYPE.END: 'e',
        MAIL_TYPE.ALL: 'abes'
    }
    """:class:`dict`: mapping from MAVIS mail type options to SGE mail options"""

    @classmethod
    def parse_qacct(cls, content):
        """
        parses the information produced by qacct

        Args:
            content (str): the content returned from the qacct command

        Raises
            ValueError: when no job information is reported (this may happen due to a bad or too old job ID where information is no longer stored)
        """
        if re.match(r'^\s*Total System Usage.*', content):
            raise ValueError('Job information not found')
        rows = []
        for section in re.split(r'=+\n', content)[1:]:  # initial item will be empty
            row = {}
            for line in section.split('\n'):
                if re.match(r'^[\s=]*$', line):
                    continue
                col, val = re.split(r'\s+', line, 1)
                val = val.strip()
                if val == 'undefined':
                    val = None
                row[col] = val

            if row['exit_status'] == '0' and row['failed'] == '0':
                status = JOB_STATUS.COMPLETED
            elif '(Killed)' in row['exit_status']:
                status = JOB_STATUS.CANCELLED
            else:
                status = JOB_STATUS.FAILED
            if ':' in row['failed']:
                status_comment = row['failed'].split(':', 1)[1].strip()
            else:
                status_comment = ''
            rows.append({
                'name': row['jobname'],
                'job_ident': row['jobnumber'],
                'task_ident': row['taskid'],
                'status': status,
                'status_comment': status_comment
            })
        return rows

    @classmethod
    def parse_qstat(cls, content):
        """
        parses the qstat content into rows/dicts representing individual jobs

        Args:
            content (str): content returned from the qstat command
        """
        header = ['job-ID', 'prior', 'name', 'user', 'state', 'submit/start at', 'queue', 'slots', 'ja-task-ID']
        content = content.strip()
        if not content:
            return []
        lines = [l for l in content.split('\n') if l.strip()]
        column_sizes = []
        for col in header:
            match = re.search(col + r'\s*', lines[0])
            if not match:
                raise ValueError('Error in parsing the qstat content for the column from', col, lines[0])
            column_sizes.append(len(match.group(0)))
        rows = []

        for line in lines[1:]:
            if re.match(r'^[\-]+$', line):
                continue  # ignore dashed separators
            row = {}
            pos = 0
            for col, size in zip(header, column_sizes):
                row[col] = line[pos:pos + size].strip()
                pos += size
            task_ident = row['ja-task-ID']
            if not task_ident or set(task_ident) & set(',:-'):
                task_ident = None
            rows.append({
                'task_ident': task_ident,
                'job_ident': row['job-ID'],
                'name': row['name'],
                'status': cls.convert_state(row['state']),
                'status_comment': ''
            })
        return rows

    @classmethod
    def convert_state(cls, state):
        states = set()
        for char in state:
            states.add(cls.STATE_MAPPING[char])
        return cumulative_job_state(states)

    def submit(self, job):
        """
        runs a subprocess sbatch command

        Args:
            job (Job): the job to be submitted
        """
        command = ['qsub', '-j', 'y']  # always join output
        if job.job_ident:
            raise ValueError('Job has already been submitted and has the job number', job.job_ident)
        if job.queue:
            command.append('-q {}'.format(job.queue))
        if job.memory_limit:
            command.extend([
                '-l',
                'mem_free={0}M,mem_token={0}M,h_vmem={0}M'.format(job.memory_limit)
            ])
        if job.time_limit:
            command.extend([
                '-l',
                'h_rt={}'.format(str(timedelta(seconds=job.time_limit)))])
        if job.import_env:
            command.append('-V')
        if job.dependencies:
            command.append(self.format_dependencies(job))
        if job.name:
            command.extend(['-N', job.name])
        if job.mail_type and job.mail_user:
            command.extend(['-m', self.MAIL_TYPE_MAPPING[job.mail_type]])
            command.extend(['-M', job.mail_user])
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            command.extend(['-t', '1-{}'.format(job.tasks)])
        if job.stdout:
            command.extend(['-o', job.stdout.format(
                name='\\${}'.format(self.ENV_JOB_NAME),
                job_ident='\\${}'.format(self.ENV_JOB_IDENT),
                task_ident='\\$TASK_ID'
            )])

        command.append(job.script)
        command = ' '.join(command)
        LOG(command, level=logging.DEBUG)
        LOG('submitting', job.name)
        content = self.command(command, shell=True)

        # example: Your job-array 3760559.1-1:1 ("MV_mock-A36971_batch-E6aEZJnTQAau598tcsMjAE") has been submitted
        # example: Your job 3766949 ("MP_batch-TvkFvM52v3ncuNQZb2M9TD") has been submitted
        match = re.match(r'^Your job(-array)? (\d+)(\.\d+-\d+:1)? .* has been submitted$', content, re.IGNORECASE)
        if not match:
            raise NotImplementedError('Error in retrieving the submitted job number. Did not match the expected pattern', content)
        job.job_ident = match.group(2)
        job.status = JOB_STATUS.SUBMITTED

        try:
            for task in job.task_list:
                task.status = job.status
                task.status_comment = job.status_comment
        except AttributeError:
            pass

    def update_info(self, job):
        """
        runs a subprocess scontrol command to get job details and add them to the current job

        Raises
            ValueError: if the job information could not be retrieved
        """
        if not job.job_ident:
            return
        command = ['qstat']
        if job.queue:
            command.extend(['-q', job.queue])
        content = self.command(command)
        rows = [row for row in self.parse_qstat(content) if row['job_ident'] == job.job_ident]

        updated = False
        if not rows:
            # job no longer scheduled
            command = ['qacct', '-j', job.job_ident]
            content = self.command(command)
            rows = self.parse_qacct(content)
            # job is still on the scheduler
        for row in rows:
            if isinstance(job, ArrayJob) and row['task_ident']:
                task_ident = int(row['task_ident'])
                task = job.get_task(task_ident)
                task.status = row['status']
                task.status_comment = row['status_comment']
            else:
                job.status = row['status']
                job.status_comment = row['status_comment']
                updated = True

        try:
            if not updated:
                job.status = cumulative_job_state([task.status for task in job.task_list])
        except AttributeError:
            pass  # only applies to array jobs

    def cancel(self, job, task_ident=None):
        """
        cancel a job
        """
        if not job.job_ident:
            return
        if task_ident is not None:
            self.command(['qdel', job.job_ident, '-t', str(task_ident)])
            job.get_task(int(task_ident)).status = JOB_STATUS.CANCELLED
            LOG('cancelled task', job.name, job.job_ident, task_ident)
        else:
            self.command(['qdel', job.job_ident])
            job.status = JOB_STATUS.CANCELLED
            LOG('cancelled job', job.name, job.job_ident)

            try:
                for task in job.task_list:
                    task.status = JOB_STATUS.CANCELLED
            except AttributeError:
                pass

    def format_dependencies(self, job):
        """
        returns a string representing the dependency argument
        """
        # special case array dependency
        try:
            if len(job.dependencies) == 1 and job.tasks == job.dependencies[0].tasks:
                dependency = job.dependencies[0]
                if not dependency.job_ident:
                    raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)
                return '-hold_jid_ad {}'.format(dependency.job_ident)
        except AttributeError:
            pass
        for dependency in job.dependencies:
            if not dependency.job_ident:
                raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)

        return '-hold_jid {}'.format(','.join([d.job_ident for d in job.dependencies]))


class TorqueScheduler(SgeScheduler):
    """
    Class for managing interactions with the Torque scheduler
    """
    NAME = SCHEDULER.TORQUE
    """:attr:`~mavis.schedule.constants.SCHEDULER`: the type of scheduler"""
    DEPENDENCY_DELIM = ':'
    ARRAY_DEPENDENCY = '-W depend=afterokarray:{}'
    JOB_DEPENDENCY = '-W depend=afterok:{}'
    ENV_TASK_IDENT = 'PBS_ARRAYID'
    ENV_JOB_IDENT = 'PBS_JOBID'
    """:class:`str`: expected pattern for environment variables which store the job id"""
    ENV_ARRAY_IDENT = ENV_JOB_IDENT
    ENV_JOB_NAME = 'PBS_JOBNAME'
    """:class:`str`: expected pattern for environment variables which store the job name"""
    TAB_SIZE = 8
    MAIL_TYPE_MAPPING = {
        MAIL_TYPE.BEGIN: 'b',
        MAIL_TYPE.NONE: 'p',
        MAIL_TYPE.FAIL: 'fa',
        MAIL_TYPE.END: 'e',
        MAIL_TYPE.ALL: 'abef'
    }
    """:class:`dict`: mapping from MAVIS mail type options to Torque mail options"""
    STATE_MAPPING = {
        'C': JOB_STATUS.COMPLETED,
        'E': JOB_STATUS.RUNNING,
        'H': JOB_STATUS.PENDING,
        'Q': JOB_STATUS.PENDING,
        'T': JOB_STATUS.RUNNING,
        'W': JOB_STATUS.PENDING,
        'S': JOB_STATUS.ERROR,
        'R': JOB_STATUS.RUNNING
    }
    """:class:`dict`: mapping from Torque job states to their MAVIS JOB_STATUS equivalent"""

    def format_dependencies(self, job):
        """
        returns a string representing the dependency argument
        """
        arr_dependencies = []
        job_dependencies = []

        for dep in job.dependencies:
            if not dep.job_ident:
                raise ValueError('An array job must be dependent only on another single array job with the same number of tasks', job, dep)

            if isinstance(dep, ArrayJob):
                arr_dependencies.append(dep.job_ident)
            else:
                job_dependencies.append(dep.job_ident)

        result = []

        if arr_dependencies:
            result.append(self.ARRAY_DEPENDENCY.format(self.DEPENDENCY_DELIM.join(arr_dependencies)))
        if job_dependencies:
            result.append(self.JOB_DEPENDENCY.format(self.DEPENDENCY_DELIM.join(job_dependencies)))

        return ','.join(result)

    @classmethod
    def parse_qstat(cls, content):
        """
        parses the qstat content into rows/dicts representing individual jobs

        Args:
            content (str): content returned from the qstat command
        """
        content = re.sub(r'\t', ' ' * cls.TAB_SIZE, content)  # PBS  torque tab size is 8
        jobs = re.split(r'\s*\n\n\s*', content.strip())
        rows = []

        for job in jobs:
            if job.startswith('request_version') or not job:
                continue
            row = {}
            lines = job.split('\n')
            task_ident = None
            row['Job Id'] = lines[0].split(':', 1)[1].strip()
            match = re.match(r'^(\d+\[)(\d+)(\].*)$', row['Job Id'])
            if match:
                row['Job Id'] = match.group(1) + match.group(3)
                task_ident = match.group(2)
            tab_size = None
            columns = []
            values = []
            for line in lines[1:]:
                if not line.strip():
                    continue
                match = re.match(r'^(\s*)(\S.*)', line)
                curr_tab_size = len(match.group(1))
                if tab_size is None:
                    tab_size = curr_tab_size

                if curr_tab_size > tab_size or '=' not in line:
                    if not values:
                        raise NotImplementedError('Unexpected indentation prior to setting column', line)
                    values[-1] = values[-1] + line.strip()
                elif curr_tab_size == tab_size:
                    col, val = line.split('=', 1)
                    columns.append(col.strip())
                    values.append(val.strip())
                else:
                    raise NotImplementedError('Unexpected indentation', line)
            for col, val in zip(columns, values):
                row[col] = val
            status = cls.STATE_MAPPING[row['job_state']]
            if status == JOB_STATUS.COMPLETED:
                if 'exit_status' in row:
                    if row['exit_status'] != '0':
                        status = JOB_STATUS.FAILED
                else:
                    status = JOB_STATUS.CANCELLED
            rows.append({
                'job_ident': row['Job Id'],
                'name': row['Job_Name'],
                'status': status,
                'task_ident': task_ident,
                'status_comment': ''
            })
        return rows

    def submit(self, job, resubmit=False):
        """
        runs a subprocess sbatch command

        Args:
            job (Job): the job to be submitted
            resubmit (bool): if true the job will be submitted even if it already has a job_ident
        """
        command = ['qsub', '-j', 'oe']  # always join output as stdout
        if job.job_ident and not resubmit:
            raise ValueError('Job has already been submitted and has the job number', job.job_ident)
        if job.queue:
            command.append('-q {}'.format(job.queue))
        if job.memory_limit:
            command.extend([
                '-l',
                'mem={0}mb'.format(job.memory_limit)
            ])
        if job.time_limit:
            command.extend([
                '-l',
                'walltime={}'.format(str(timedelta(seconds=job.time_limit)))])
        if job.import_env:
            command.append('-V')
        if job.dependencies:
            command.append(self.format_dependencies(job))
        if job.name:
            command.extend(['-N', job.name])
        if job.stdout:
            command.extend(['-o', job.stdout.format(
                name='${}'.format(self.ENV_JOB_NAME),
                job_ident='${}'.format(self.ENV_JOB_IDENT),
                task_ident='${}'.format(self.ENV_TASK_IDENT)
            )])
        if job.mail_type and job.mail_user:
            command.extend(['-m', job.mail_type])
            command.extend(['-M', job.mail_user])
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            concurrency_limit = '' if self.concurrency_limit is None else '%{}'.format(self.concurrency_limit)
            command.extend(['-t', '{}{}'.format(','.join([t.task_ident for t in job.task_list]), concurrency_limit)])

        command.append(job.script)
        LOG('submitting', job.name)
        content = self.command(command)

        job.job_ident = content.strip()
        job.status = JOB_STATUS.SUBMITTED
        job.status_comment = ''

        # update task status
        try:
            for task in job.task_list:
                task.status = job.status
                task.status_comment = job.status_comment
        except AttributeError:
            pass

    def update_info(self, job):
        """
        runs a subprocess scontrol command to get job details and add them to the current job

        Raises
            ValueError: if the job information could not be retrieved
        """
        if job.job_ident is None:
            job.status = JOB_STATUS.NOT_SUBMITTED
            return
        command = ['qstat', '-f', job.job_ident]
        if isinstance(job, ArrayJob):
            command.append('-t')
        content = self.command(command)
        rows = self.parse_qstat(content)
        tasks_updated = False

        for row in rows:
            if row['job_ident'] != job.job_ident:
                continue
            if isinstance(job, ArrayJob) and row['task_ident']:
                task_ident = int(row['task_ident'])
                task = job.get_task(task_ident)
                task.status = row['status']
                task.status_comment = row['status_comment']
                tasks_updated = True
            else:
                job.status = row['status']
                job.status_comment = row['status_comment']

        if tasks_updated:
            job.status = cumulative_job_state([t.status for t in job.task_list])

    def cancel(self, job, task_ident=None):
        """
        cancel a job
        """
        if not job.job_ident:
            return
        if task_ident is not None:
            self.command(['qdel', job.job_ident, '-t', str(task_ident)])
            job.get_task(int(task_ident)).status = JOB_STATUS.CANCELLED
            LOG('cancelled task', job.name, job.job_ident, task_ident)
        else:
            self.command(['qdel', job.job_ident])
            job.status = JOB_STATUS.CANCELLED
            LOG('cancelled job', job.name, job.job_ident)

            try:
                for task in job.task_list:
                    task.status = JOB_STATUS.CANCELLED
            except AttributeError:
                pass
