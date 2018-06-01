import subprocess
import unittest
from unittest import mock

from mavis.schedule import scheduler as _scheduler
from mavis.schedule import constants as _constants
from mavis.schedule import job as _job
from mavis.constants import SUBCOMMAND

class TestParseQstat(unittest.TestCase):

    # TODO: single job running
    # TODO: batch job running
    # TODO: single job complete
    def test_single_job_complete(self):
        content = """
Job Id: 9.torque01.bcgsc.ca
    Job_Name = subtest.sh
    Job_Owner = creisle@torque01.bcgsc.ca
    resources_used.cput = 00:00:00
    resources_used.vmem = 346716kb
    resources_used.walltime = 00:01:00
    resources_used.mem = 3624kb
    resources_used.energy_used = 0
    job_state = C
    queue = batch
    server = torque01.bcgsc.ca
    Checkpoint = u
    ctime = Tue May 29 09:37:00 2018
    Error_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/workspa
        ce/creisle/temp/subtest.sh.e9
    exec_host = torque01.bcgsc.ca/0
    Hold_Types = n
    Join_Path = n
    Keep_Files = n
    Mail_Points = a
    mtime = Tue May 29 09:38:01 2018
    Output_Path = torque01.bcgsc.ca:/projects/trans_scratch/validations/worksp
        ace/creisle/temp/subtest.sh.o9
    Priority = 0
    qtime = Tue May 29 09:37:00 2018
    Rerunable = True
    Resource_List.walltime = 01:00:00
    Resource_List.nodes = 1
    Resource_List.nodect = 1
    session_id = 25438
    Variable_List = PBS_O_QUEUE=batch,PBS_O_HOME=/home/creisle,
        PBS_O_LOGNAME=creisle,
        PBS_O_PATH=/home/creisle/applications/node-v10.1.0-linux-x64/bin:/hom
        e/creisle/.npm-packages/bin:/home/creisle/bin:/home/creisle/applicatio
        ns/centos06/python-3.6.1/bin:/projects/tumour_char/analysis_scripts/bi
        n/pog:/gsc/software/linux-x86_64-centos6/git-2.12.0/bin/:/usr/local/bi
        n:/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/p
        rojects/trans_scratch/software/pipeline_commands/:/home/creisle/bin,
        PBS_O_MAIL=/var/spool/mail/creisle,PBS_O_SHELL=/bin/bash,
        PBS_O_LANG=en_US.UTF-8,
        PBS_O_WORKDIR=/projects/trans_scratch/validations/workspace/creisle/t
        emp,PBS_O_HOST=torque01.bcgsc.ca,PBS_O_SERVER=torque01.bcgsc.ca
    euser = creisle
    egroup = users
    queue_type = E
    comment = Job started on Tue May 29 at 09:37
    etime = Tue May 29 09:37:00 2018
    exit_status = 0
    submit_args = subtest.sh
    start_time = Tue May 29 09:37:01 2018
    start_count = 1
    fault_tolerant = False
    comp_time = Tue May 29 09:38:01 2018
    job_radix = 0
    total_runtime = 60.481239
    submit_host = torque01.bcgsc.ca
    init_work_dir = /projects/trans_scratch/validations/workspace/creisle/temp

    request_version = 1

        """
        rows = _scheduler.TorqueScheduler().parse_qstat(content)
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual(_constants.JOB_STATUS.COMPLETED, row['status'])
        self.assertEqual('9.torque01.bcgsc.ca', row['job_ident'])
        self.assertEqual('subtest.sh', row['name'])
        self.assertIs(None, row['task_ident'])
        self.assertEqual('', row['status_comment'])

    # TODO: batch job complete
    # TODO: single job error
    # TODO: batch job error
    # TODO: single job exiting
    # TODO: batch job exiting

class TestCancel(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_single_job(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_array_job(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        for task in job.task_list:
            self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_array_job_task(self, patcher):
        sched = _scheduler.TorqueScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job, task_ident='4')
        self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, job.status)
        for i, task in enumerate(job.task_list):
            if i == 3:
                self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
            else:
                self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, task.status)
        patcher.assert_called_with(['qdel', '1234', '-t', '4'])

    @mock.patch('mavis.schedule.scheduler.TorqueScheduler.command')
    def test_bad_command(self, patcher):
        patcher.side_effect = [subprocess.CalledProcessError(1, 'cmd')]
        sched = _scheduler.TorqueScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        with self.assertRaises(subprocess.CalledProcessError):
            sched.cancel(job)
        patcher.assert_called_with(['qdel', '1234'])