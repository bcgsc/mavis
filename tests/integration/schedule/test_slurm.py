import unittest
from unittest import mock
import configparser
import tempfile
import shutil
import os
import subprocess

from mavis.schedule import pipeline as _pipeline
from mavis.schedule import job as _job
from mavis.schedule import constants as _constants
from mavis.config import MavisConfig
from mavis.main import main

from ...util import get_data


class TestSubmit(unittest.TestCase):

    # TODO: test initial submission
    # TODO: test submit after failure
    # TODO: test reporting errors

    @mock.patch('subprocess.check_output')
    def test_single_job(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.Job(
            name='job1',
            stage='validate',
            script='submit.sh'
        )
        print(job)
        _job.SlurmScheduler.submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'job-%x-%j.log',
            '--mail-type=NONE',
            'submit.sh'
        ])

    @mock.patch('subprocess.check_output')
    def test_dependent_job(self, patch_check):
        patch_check.side_effect = ["Submitted batch job 1665695".encode('utf8')]
        job = _job.Job(
            name='job1',
            stage='validate',
            script='submit.sh',
            dependencies=[_job.Job(
                name='job2',
                stage='cluster',
                script='submit2.sh',
                job_ident='12345678'
            )]
        )
        print(job)
        _job.SlurmScheduler.submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000',
            '-t', '16:00:00',
            '--export=ALL',
            '--dependency=afterok:12345678',
            '-J', 'job1',
            '-o', 'job-%x-%j.log',
            '--mail-type=NONE',
            'submit.sh'
        ])

    @mock.patch('subprocess.check_output')
    def test_cascade(self, patch_check):
        patch_check.side_effect = [
            "Submitted batch job 12345678".encode('utf8'),
            "Submitted batch job 1665695".encode('utf8')
        ]
        job = _job.Job(
            name='job1',
            stage='validate',
            script='submit.sh',
            dependencies=[_job.Job(
                name='job2',
                stage='cluster',
                script='submit2.sh'
            )]
        )
        print(job)
        _job.SlurmScheduler.submit(job, cascade=True)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000',
            '-t', '16:00:00',
            '--export=ALL',
            '--dependency=afterok:12345678',
            '-J', 'job1',
            '-o', 'job-%x-%j.log',
            '--mail-type=NONE',
            'submit.sh'
        ])

    @mock.patch('subprocess.check_output')
    def test_no_cascade_error(self, patch_check):
        patch_check.side_effect = [
            "Submitted batch job 12345678".encode('utf8'),
            "Submitted batch job 1665695".encode('utf8')
        ]
        job = _job.Job(
            name='job1',
            stage='validate',
            script='submit.sh',
            dependencies=[_job.Job(
                name='job2',
                stage='cluster',
                script='submit2.sh'
            )]
        )
        print(job)
        with self.assertRaises(ValueError):
            _job.SlurmScheduler.submit(job, cascade=False)

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.ArrayJob(
            name='job1',
            stage='validate',
            script='submit.sh',
            tasks=10
        )
        print(job)
        _job.SlurmScheduler.submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'job-%x-%A_%a.log',
            '--mail-type=NONE',
            '--array=1-10',
            'submit.sh'
        ])

    @mock.patch('subprocess.check_output')
    def test_job_array_concurrency_limit(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.ArrayJob(
            name='job1',
            stage='validate',
            script='submit.sh',
            tasks=10,
            concurrency_limit=2
        )
        print(job)
        _job.SlurmScheduler.submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'job-%x-%A_%a.log',
            '--mail-type=NONE',
            '--array=1-10%2',
            'submit.sh'
        ])


class TestStatus(unittest.TestCase):
    # TODO: status of array job
    # TODO: status of single job
    # TODO: status of job waiting on dependency

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        content = """
JobId=1668211 ArrayJobId=1668211 ArrayTaskId=2 JobName=name
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=27 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:07 TimeLimit=3-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-22T19:27:50 EligibleTime=2018-05-22T19:27:51
   StartTime=2018-05-22T19:28:10 EndTime=2018-05-25T19:28:10 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=short AllocNode:Sid=n104:55950
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n305
   BatchHost=n305
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/job-%x-1668211_2.log
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/job-%x-1668211_2.log
   Power=

JobId=1668212 ArrayJobId=1668211 ArrayTaskId=1 JobName=name
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=27 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:07 TimeLimit=3-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-22T19:27:50 EligibleTime=2018-05-22T19:27:51
   StartTime=2018-05-22T19:28:10 EndTime=2018-05-25T19:28:10 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=short AllocNode:Sid=n104:55950
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n305
   BatchHost=n305
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/job-%x-1668211_1.log
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/job-%x-1668211_1.log
   Power=
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            job_ident='1668211',
            tasks=2,
            stage='validate'
        )
        hist = _job.SlurmScheduler.status(job)
        self.assertEqual(1, len(hist))
        self.assertEqual(_constants.JOB_STATUS.RUNNING, list(hist.keys())[0])
        self.assertEqual(2, list(hist.values())[0])

    @mock.patch('subprocess.check_output')
    def test_single_job(self, patch_check):
        content = """
JobId=1668282 JobName=subtest.sh
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=27 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:39 TimeLimit=90-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-22T19:47:04 EligibleTime=2018-05-22T19:47:04
   StartTime=2018-05-22T19:47:14 EndTime=2018-08-20T19:47:14 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=all AllocNode:Sid=n104:55950
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n114
   BatchHost=n114
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1668282.out
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1668282.out
   Power=

        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            job_ident='1668211',
            tasks=2,
            stage='validate'
        )
        hist = _job.SlurmScheduler.status(job)
        self.assertEqual(1, len(hist))
        self.assertEqual(_constants.JOB_STATUS.RUNNING, list(hist.keys())[0])
        self.assertEqual(1, list(hist.values())[0])


