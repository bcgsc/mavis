from datetime import timedelta, datetime
import subprocess
import re
import logging

from ..util import LOG

from .job import Job, ArrayJob, Task
from .constants import SCHEDULER, JOB_STATUS, cumulative_job_state, MAIL_TYPE


class Scheduler:
    """
    Class responsible for methods interacting with the scheduler
    """
    ARRAY_DEPENDENCY = None
    JOB_DEPENDENCY = None
    CORR_ARRAY_DEPENDENCY = None
    ENV_TASK_IDENT = 'TASK_IDENT'
    DEPENDENCY_DELIM = ':'
    HEADER_PREFIX = '#'

    @classmethod
    def submit(cls, job, task_ident=None, cascade=False):
        raise NotImplementedError('abstract method')

    @classmethod
    def update_info(cls, job):
        """
        update the information about the job from the scheduler
        """
        raise NotImplementedError('abstract method')

    @classmethod
    def cancel(cls, job):
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
        for dependency in job.dependencies:
            if not dependency.job_ident:
                if not cascade:
                    raise ValueError('The dependencies must be submitted before the dependent job', job, dependency)
                cls.submit(dependency, task_ident, cascade=cascade)
        return cls.JOB_DEPENDENCY.format(cls.DEPENDENCY_DELIM.join([d.job_ident for d in job.dependencies]))


class SlurmScheduler(Scheduler):
    """
    Class for formatting commands to match a SLURM scheduler system
    SLURM docs can be found here https://slurm.schedmd.com
    """
    NAME = SCHEDULER.SLURM
    ARRAY_DEPENDENCY = '--dependency=aftercorr:{}'
    JOB_DEPENDENCY = '--dependency=afterok:{}'
    ENV_TASK_IDENT = 'SLURM_ARRAY_TASK_ID'

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
    def parse_sacct(cls, content):
        """
        parses content returned from the sacct command
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
            if '_' in row['JobID']:
                job_ident, task_ident = row['JobID'].rsplit('_', 1)
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
            rows.append({
                'job_ident': row['JobId'],
                'status': row['JobState'],
                'name': row['JobName'],
                'status_comment': row['Reason'] if row['Reason'].lower() != 'none' else '',
                'task_ident': row.get('ArrayTaskId', None)
            })
        return rows

    @classmethod
    def get_job_ident(cls, name):
        command = ['sacct', '--name', name, '--parsable2']

    @classmethod
    def update_info(cls, job):
        if not job.job_ident:
            return
        command = ['sacct', '-j', job.job_ident, '--long', '--parsable2']
        #else:
        #    start_date = datetime.fromtimestamp(job.created_at).strftime('%Y-%m-%d')
        #    command = ['sacct', '--name', job.name, '--long', '--parsable2', '-S', start_date]
        content = subprocess.check_output(command).decode('utf8')
        rows = cls.parse_sacct(content)
        updated = False

        for row in rows:
            if row['job_ident'] == job.job_ident:
                if row['task_ident'] is not None:
                    job.task_list[int(row['task_ident']) - 1].status = row['status']
                    job.task_list[int(row['task_ident']) - 1].status_comment = row['status_comment']
                else:
                    job.status = row['status']
                    job.status_comment = row['status_comment']
                    updated = True
        try:
            if not updated:
                job.status = cumulative_job_state([t.status for t in job.task_list])
        except AttributeError:
            pass

    @classmethod
    def cancel(cls, job):
        command = ['scancel', job.job_ident]
        subprocess.check_output(command)
        job.job_ident = None
        job.status = JOB_STATUS.CANCELLED


class SgeScheduler(Scheduler):
    NAME = SCHEDULER.SGE
    ARRAY_DEPENDENCY = '-hold_jid_ad {}'
    JOB_DEPENDENCY = '-hold_jid {}'
    ENV_TASK_IDENT = 'SGE_TASK_ID'
    ENV_JOB_IDENT = 'JOB_ID'
    ENV_JOB_NAME = 'JOB_NAME'
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

    MAIL_TYPE_MAPPING = {
        MAIL_TYPE.BEGIN: 'b',
        MAIL_TYPE.NONE: 'n',
        MAIL_TYPE.FAIL: 'as',
        MAIL_TYPE.END: 'e',
        MAIL_TYPE.ALL: 'abes'
    }

    @classmethod
    def parse_qacct(cls, content):
        """
        parses the information produced by qacct

        Args
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
                col, val = re.split('\s+', line, 1)
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

        Args
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
                'mem_free={0}M,mem_token={0}M,h_vmem={0}M'.format(job.memory_limit)
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
            command.extend(['-N', job.name])
        if job.mail_type and job.mail_user:
            command.extend(['-m', cls.MAIL_TYPE_MAPPING[job.mail_type]])
            command.extend(['-M', job.mail_user])
        # options specific to job arrays
        if isinstance(job, ArrayJob):
            if task_ident is None:  # default to all
                command.extend(['-t', '1-{}'.format(job.tasks)])
            else:
                command.append(['-t', str(task_ident)])
        if job.stdout:
            command.extend(['-o', job.stdout.format(
                name='\${}'.format(cls.ENV_JOB_NAME),
                job_ident='\${}'.format(cls.ENV_JOB_IDENT),
                task_ident='\$TASK_ID'
            )])

        command.append(job.script)
        command = ' '.join(command)
        LOG(command, level=logging.DEBUG)
        content = subprocess.check_output(command, shell=True).decode('utf8').strip()

        # example: Your job-array 3760559.1-1:1 ("MV_mock-A36971_batch-E6aEZJnTQAau598tcsMjAE") has been submitted
        # example: Your job 3766949 ("MP_batch-TvkFvM52v3ncuNQZb2M9TD") has been submitted
        match = re.match(r'^Your job(-array)? (\d+)(\.\d+-\d+:1)? .* has been submitted$', content, re.IGNORECASE)
        if not match:
            raise NotImplementedError('Error in retrieving the submitted job number. Did not match the expected pattern', content)
        job.job_ident = match.group(2)
        job.status = JOB_STATUS.SUBMITTED

    @classmethod
    def update_info(cls, job):
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
        content = subprocess.check_output(command).decode('utf8').strip()
        rows = [row for row in cls.parse_qstat(content) if row['job_ident'] == job.job_ident]

        updated = False
        if not rows:
            # job no longer scheduled
            command = ['qacct', '-j', job.job_ident]
            content = subprocess.check_output(command).decode('utf8').strip()
            rows = cls.parse_qacct(content)
            # job is still on the scheduler
        for row in rows:
            if isinstance(job, ArrayJob) and row['task_ident']:
                task_ident = int(row['task_ident'])
                job.task_list[task_ident - 1].status = row['status']
                job.task_list[task_ident - 1].status_comment = row['status_comment']
            else:
                job.status = row['status']
                job.status_comment = row['status_comment']
                updated = True

        try:
            if not updated:
                job.status = cumulative_job_state([task.status for task in job.task_list])
        except AttributeError:
            pass  # only applies to array jobs


class TorqueScheduler(SgeScheduler):
    NAME = SCHEDULER.TORQUE
    DEPENDENCY_DELIM = ':'
    ARRAY_DEPENDENCY = '-W dependency=afterokarray:{}'
    JOB_DEPENDENCY = '-W dependency=afterok:{}'
    ENV_TASK_IDENT = 'PBS_ARRAYID'
    ENV_JOB_IDENT = 'PBS_JOBID'
    ENV_JOB_NAME = 'PBS_JOBNAME'
    MAIL_TYPE_MAPPING = {
        MAIL_TYPE.BEGIN: 'b',
        MAIL_TYPE.NONE: 'p',
        MAIL_TYPE.FAIL: 'fa',
        MAIL_TYPE.END: 'e',
        MAIL_TYPE.ALL: 'abef'
    }

    @classmethod
    def submit(cls, job, task_ident=None, cascade=False):
        """
        runs a subprocess sbatch command

        Args
            job (Job): the job to be submitted
            task_ident (int): submit only a particular task of the current job
            cascade (bool): submit dependencies if they have not yet been submitted
        """
        command = ['qsub', '-j', 'oe']  # always join output as stdout
        if job.job_ident:
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
            command.append(cls.format_dependencies(job, task_ident=task_ident, cascade=cascade))
        if job.name:
            command.extend(['-N', job.name])
        if job.stdout:
            command.extend(['-o', job.stdout.format(
                name='\${}'.format(cls.ENV_JOB_NAME),
                job_ident='\${}'.format(cls.ENV_JOB_IDENT),
                task_ident='\${}'.format(cls.ENV_TASK_IDENT)
            )])
        if job.mail_type and job.mail_user:
            command.extend(['-m', job.mail_type])
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
    def update_info(cls, job):
        """
        runs a subprocess scontrol command to get job details and add them to the current job

        Raises
            ValueError: if the job information could not be retrieved
        """
        command = ['qstat']
        if job.queue:
            command.extend(['-q', job.queue])
        content = subprocess.check_output(command).decode('utf8').strip()
        rows = [row for row in cls.parse_qstat(content) if row['job_ident'] == job.job_ident]

        updated = False
        if not rows:
            # job no longer scheduled
            command = ['qacct', '-j', job.job_ident]
            content = subprocess.check_output(command).decode('utf8').strip()
            rows = cls.parse_qacct(content)
            # job is still on the scheduler
        for row in rows:
            if isinstance(job, ArrayJob) and row['task_ident']:
                task_ident = int(row['task_ident'])
                job.task_list[task_ident - 1].status = row['status']
                job.task_list[task_ident - 1].status_comment = row['status_comment']
            else:
                job.status = row['status']
                job.status_comment = row['status_comment']
                updated = True

        try:
            if not updated:
                job.status = cumulative_job_state([task.status for task in job.task_list])
        except AttributeError:
            pass  # only applies to array jobs

    @classmethod
    def cancel(cls, job):
        command = ['scancel', job.job_ident]
        subprocess.check_output(command)
        job.job_ident = None
        job.status = JOB_STATUS.CANCELLED


class LocalScheduler(Scheduler):
    NAME = SCHEDULER.NONE


