import subprocess
import unittest
from unittest import mock

from mavis.schedule import job as _job
from mavis.schedule import constants as _constants
from mavis.schedule import scheduler as _scheduler
from mavis.constants import SUBCOMMAND


class TestSubmit(unittest.TestCase):

    # TODO: test initial submission
    # TODO: test submit after failure
    # TODO: test reporting errors

    @mock.patch('subprocess.check_output')
    def test_single_job(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.Job(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh'
        )
        print(job)
        _scheduler.SlurmScheduler().submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'temp/job-%x-%j.log',
            'submit.sh'
        ], shell=False)

    @mock.patch('subprocess.check_output')
    def test_partition(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.Job(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            queue='all'
        )
        print(job)
        _scheduler.SlurmScheduler().submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--partition=all',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'temp/job-%x-%j.log',
            'submit.sh'
        ], shell=False)

    @mock.patch('subprocess.check_output')
    def test_mail_options(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.Job(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL
        )
        print(job)
        _scheduler.SlurmScheduler().submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'temp/job-%x-%j.log',
            '--mail-type=ALL',
            '--mail-user=me@example.com',
            'submit.sh'
        ], shell=False)

    @mock.patch('subprocess.check_output')
    def test_dependent_job(self, patch_check):
        patch_check.side_effect = ["Submitted batch job 1665695".encode('utf8')]
        job = _job.Job(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            dependencies=[_job.Job(
                output_dir='temp',
                name='job2',
                stage='cluster',
                script='submit2.sh',
                job_ident='12345678'
            )]
        )
        print(job)
        _scheduler.SlurmScheduler().submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '--dependency=afterok:12345678',
            '-J', 'job1',
            '-o', 'temp/job-%x-%j.log',
            'submit.sh'
        ], shell=False)

    @mock.patch('subprocess.check_output')
    def test_dependency_error(self, patch_check):
        patch_check.side_effect = [
            "Submitted batch job 12345678".encode('utf8'),
            "Submitted batch job 1665695".encode('utf8')
        ]
        job = _job.Job(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            dependencies=[_job.Job(
                output_dir='temp',
                name='job2',
                stage='cluster',
                script='submit2.sh'
            )]
        )
        print(job)
        with self.assertRaises(ValueError):
            _scheduler.SlurmScheduler().submit(job)

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        patch_check.return_value = "Submitted batch job 1665695".encode('utf8')
        job = _job.ArrayJob(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            task_list=10
        )
        print(job)
        _scheduler.SlurmScheduler().submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        patch_check.assert_called_with([
            'sbatch',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'temp/job-%x-%A-%a.log',
            '--array=1-10',
            'submit.sh'
        ], shell=False)

    @mock.patch('subprocess.check_output')
    def test_job_array_concurrency_limit(self, patch_check):
        patch_check.side_effect = ["Submitted batch job 1665695".encode('utf8')]
        print(patch_check)
        job = _job.ArrayJob(
            output_dir='temp',
            name='job1',
            stage='validate',
            script='submit.sh',
            task_list=[1, 2, 3, 4, 5, 14, 16]
        )
        _scheduler.SlurmScheduler(concurrency_limit=2).submit(job)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        self.assertEqual('1665695', job.job_ident)
        exp = [
            'sbatch',
            '--mem', '16000M',
            '-t', '16:00:00',
            '--export=ALL',
            '-J', 'job1',
            '-o', 'temp/job-%x-%A-%a.log',
            '--array=1-5,14,16%2',
            'submit.sh'
        ]
        patch_check.assert_called_with(exp, shell=False)


class TestUpdate(unittest.TestCase):
    # TODO: status of array job
    # TODO: status of single job
    # TODO: status of job waiting on dependency

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        content = """
JobID|JobIDRaw|JobName|Partition|MaxVMSize|MaxVMSizeNode|MaxVMSizeTask|AveVMSize|MaxRSS|MaxRSSNode|MaxRSSTask|AveRSS|MaxPages|MaxPagesNode|MaxPagesTask|AvePages|MinCPU|MinCPUNode|MinCPUTask|AveCPU|NTasks|AllocCPUS|Elapsed|State|ExitCode|AveCPUFreq|ReqCPUFreqMin|ReqCPUFreqMax|ReqCPUFreqGov|ReqMem|ConsumedEnergy|MaxDiskRead|MaxDiskReadNode|MaxDiskReadTask|AveDiskRead|MaxDiskWrite|MaxDiskWriteNode|MaxDiskWriteTask|AveDiskWrite|AllocGRES|ReqGRES|ReqTRES|AllocTRES|
1671879_1|1671879|MV_mock-A36971_batch-tX8SW6tEiEfZ8ZLHDPDa83|short||||||||||||||||||1|00:00:00|FAILED|1:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671879_1.batch|1671879.batch|batch||||||||||||||||||1|1|00:00:00|FAILED|1:0||0|0|0|16000Mn|||||||||||||cpu=1,mem=16000M,node=1|
1671880_1|1671880|MV_mock-A47933_batch-tX8SW6tEiEfZ8ZLHDPDa83|short||||||||||||||||||1|00:00:00|FAILED|1:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671880_1.batch|1671880.batch|batch||||||||||||||||||1|1|00:00:00|FAILED|1:0||0|0|0|18000Mn|||||||||||||cpu=1,mem=18000M,node=1|
1671893_1|1671893|MV_mock-A36971_batch-tX8SW6tEiEfZ8ZLHDPDa83|short||||||||||||||||||1|00:00:01|FAILED|1:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671893_1.batch|1671893.batch|batch||||||||||||||||||1|1|00:00:01|FAILED|1:0||0|0|0|16000Mn|||||||||||||cpu=1,mem=16000M,node=1|
1671894_1|1671894|MV_mock-A47933_batch-tX8SW6tEiEfZ8ZLHDPDa83|short||||||||||||||||||1|00:00:00|FAILED|1:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671894_1.batch|1671894.batch|batch||||||||||||||||||1|1|00:00:00|FAILED|1:0||0|0|0|18000Mn|||||||||||||cpu=1,mem=18000M,node=1|
1671915_1|1671915|MV_mock-A36971_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:20|CANCELLED by 1365|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671915_1.batch|1671915.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:20|CANCELLED|0:15|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1671916_1|1671916|MV_mock-A47933_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:20|CANCELLED by 1365|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671916_1.batch|1671916.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:20|CANCELLED|0:15|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1671970_1|1671970|MV_mock-A36971_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:21|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671970_1.batch|1671970.batch|batch||125588K|n305|0|125588K|908K|n305|0|908K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:21|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1671971_1|1671971|MV_mock-A47933_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:20|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671971_1.batch|1671971.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:20|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1671974_1|1671974|MV_mock-A36971_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:11|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671974_1.batch|1671974.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:11|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1671975_1|1671975|MV_mock-A47933_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:10|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671975_1.batch|1671975.batch|batch||125588K|n305|0|125588K|908K|n305|0|908K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:10|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1671981_1|1671981|MV_mock-A36971_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:12|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671981_1.batch|1671981.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:12|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1671982_1|1671982|MV_mock-A47933_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:11|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1671982_1.batch|1671982.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:11|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1671983_1|1671983|MA_mock-A36971_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:05|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1671983_1.batch|1671983.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:05|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1671984_1|1671984|MA_mock-A47933_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:04|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1671984_1.batch|1671984.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:04|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1671985|1671985|MP_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:02|FAILED|2:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1671985.batch|1671985.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|FAILED|2:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1671986|1671986|MS_batch-ezPmnHmYjZjsj8gfCynbsX|short||||||||||||||||||1|00:00:00|CANCELLED by 1365|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1||
1672141|1672141|subtest.sh|all||||||||||||||||||1|00:01:01|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672141.batch|1672141.batch|batch||207004K|n106|0|207004K|1760K|n106|0|1760K|1K|n106|0|1K|00:00:00|n106|0|00:00:00|1|1|00:01:01|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672166|1672166|subtest.sh|all||||||||||||||||||1|00:01:03|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672166.batch|1672166.batch|batch||207004K|n106|0|207004K|1760K|n106|0|1760K|0|n106|0|0|00:00:00|n106|0|00:00:00|1|1|00:01:03|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_10|1672169|subtest.sh|all||||||||||||||||||1|00:01:02|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_10.batch|1672169.batch|batch||207004K|n130|0|207004K|1764K|n130|0|1764K|0|n130|0|0|00:00:00|n130|0|00:00:00|1|1|00:01:02|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_1|1672171|subtest.sh|all||||||||||||||||||1|00:01:00|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_1.batch|1672171.batch|batch||207004K|n106|0|207004K|1764K|n106|0|1764K|0|n106|0|0|00:00:00|n106|0|00:00:00|1|1|00:01:00|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_2|1672172|subtest.sh|all||||||||||||||||||1|00:01:00|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_2.batch|1672172.batch|batch||207004K|n106|0|207004K|1764K|n106|0|1764K|0|n106|0|0|00:00:00|n106|0|00:00:00|1|1|00:01:00|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_3|1672173|subtest.sh|all||||||||||||||||||1|00:01:01|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_3.batch|1672173.batch|batch||207004K|n106|0|207004K|1764K|n106|0|1764K|0|n106|0|0|00:00:00|n106|0|00:00:00|1|1|00:01:01|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_4|1672174|subtest.sh|all||||||||||||||||||1|00:01:01|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_4.batch|1672174.batch|batch||207004K|n106|0|207004K|1760K|n106|0|1760K|0|n106|0|0|00:00:00|n106|0|00:00:00|1|1|00:01:01|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n106|0|0.06M|0.00M|n106|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_5|1672175|subtest.sh|all||||||||||||||||||1|00:01:02|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_5.batch|1672175.batch|batch||207004K|n130|0|207004K|1764K|n130|0|1764K|0|n130|0|0|00:00:00|n130|0|00:00:00|1|1|00:01:02|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_6|1672176|subtest.sh|all||||||||||||||||||1|00:01:02|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_6.batch|1672176.batch|batch||207004K|n130|0|207004K|1764K|n130|0|1764K|1K|n130|0|1K|00:00:00|n130|0|00:00:00|1|1|00:01:02|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_7|1672177|subtest.sh|all||||||||||||||||||1|00:01:01|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_7.batch|1672177.batch|batch||207004K|n130|0|207004K|1756K|n130|0|1756K|1K|n130|0|1K|00:00:00|n130|0|00:00:00|1|1|00:01:01|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_8|1672178|subtest.sh|all||||||||||||||||||1|00:01:02|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_8.batch|1672178.batch|batch||207004K|n130|0|207004K|1764K|n130|0|1764K|0|n130|0|0|00:00:00|n130|0|00:00:00|1|1|00:01:02|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672169_9|1672179|subtest.sh|all||||||||||||||||||1|00:01:02|COMPLETED|0:0||Unknown|Unknown|Unknown|7000Mc||||||||||||cpu=1,mem=7000M,node=1|cpu=1,mem=7000M,node=1|
1672169_9.batch|1672179.batch|batch||207004K|n130|0|207004K|1760K|n130|0|1760K|1K|n130|0|1K|00:00:00|n130|0|00:00:00|1|1|00:01:02|COMPLETED|0:0|2.19M|0|0|0|7000Mc|0|0.06M|n130|0|0.06M|0.00M|n130|0|0.00M||||cpu=1,mem=7000M,node=1|
1672268_2|1672268|MV_mock-A36971_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:36|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672268_2.batch|1672268.batch|batch||2025280K|n305|0|2025280K|58776K|n305|0|58776K|0|n305|0|0|00:00:03|n305|0|00:00:03|1|1|00:00:36|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|13.01M|n305|0|13.01M|0.05M|n305|0|0.05M||||cpu=1,mem=16000M,node=1|
1672269_3|1672269|MV_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:33|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672269_3.batch|1672269.batch|batch||2019424K|n305|0|2019424K|55844K|n305|0|55844K|0|n305|0|0|00:00:03|n305|0|00:00:03|1|1|00:00:33|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|11.60M|n305|0|11.60M|0.05M|n305|0|0.05M||||cpu=1,mem=18000M,node=1|
1672270_2|1672270|MA_mock-A36971_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672270_2.batch|1672270.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672271_3|1672271|MA_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672271_3.batch|1672271.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672272|1672272|MP_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:02|FAILED|2:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672272.batch|1672272.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|FAILED|2:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1672273|1672273|MS_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:00|CANCELLED by 1365|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1||
1672270_1|1672274|MA_mock-A36971_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:04|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672270_1.batch|1672274.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:04|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672271_1|1672275|MA_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672271_1.batch|1672275.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672271_2|1672276|MA_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672271_2.batch|1672276.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672268_1|1672277|MV_mock-A36971_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:32|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672268_1.batch|1672277.batch|batch||2018912K|n305|0|2018912K|51952K|n305|0|51952K|0|n305|0|0|00:00:03|n305|0|00:00:03|1|1|00:00:32|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|10.79M|n305|0|10.79M|0.05M|n305|0|0.05M||||cpu=1,mem=16000M,node=1|
1672269_1|1672278|MV_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:32|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672269_1.batch|1672278.batch|batch||2018660K|n305|0|2018660K|57016K|n305|0|57016K|0|n305|0|0|00:00:03|n305|0|00:00:03|1|1|00:00:32|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|12.79M|n305|0|12.79M|0.05M|n305|0|0.05M||||cpu=1,mem=18000M,node=1|
1672269_2|1672279|M1673291_mock-A47933_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:33|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672269_2.batch|1671673291279.batch|batch||2019424K|n305|0|2019424K|54212K|n305|0|54212K|0|n305|0|0|00:00:03|n305|0|00:00:03|1|1|00:00:33|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|11.51M|n305|0|11.51M|0.05M|n305|0|0.05M||||cpu=1,mem=18000M,node=1|
1672454_2|1672454|M1673291_mock-A36971_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:29|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672454_2.batch|1671673291454.batch|batch||125588K|n305|0|125588K|908K|n305|0|908K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:29|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1672455_3|1672455|MV_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:26|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672455_3.batch|1672455.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:26|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1672456_2|1672456|MA_mock-A36971_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672456_2.batch|1672456.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672457_3|1672457|MA_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:03|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672457_3.batch|1672457.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:03|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672458|1672458|MP_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672458.batch|1672458.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1672456_1|1672459|MA_mock-A36971_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:04|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672456_1.batch|1672459.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:04|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672457_1|1672460|MA_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:03|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672457_1.batch|1672460.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:03|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672457_2|1672461|MA_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|12000Mn||||||||||||cpu=1,mem=12000M,node=1|cpu=1,mem=12000M,node=1|
1672457_2.batch|1672461.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|12000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=12000M,node=1|
1672454_1|1672462|MV_mock-A36971_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:25|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672454_1.batch|1672462.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:25|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
1672455_1|1672463|MV_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:25|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672455_1.batch|1672463.batch|batch||125588K|n305|0|125588K|904K|n305|0|904K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:25|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1672455_2|1672464|MV_mock-A47933_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:27|COMPLETED|0:0||Unknown|Unknown|Unknown|18000Mn||||||||||||cpu=1,mem=18000M,node=1|cpu=1,mem=18000M,node=1|
1672455_2.batch|1672464.batch|batch||125588K|n305|0|125588K|900K|n305|0|900K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:27|COMPLETED|0:0|2.19M|0|0|0|18000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=18000M,node=1|
1672465|1672465|MS_batch-uKEUyUuWbi2mgd75KjP4k5|short||||||||||||||||||1|00:00:02|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1|
1672465.batch|1672465.batch|batch||125588K|n305|0|125588K|896K|n305|0|896K|0|n305|0|0|00:00:00|n305|0|00:00:00|1|1|00:00:02|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|0|n305|65534|0|0|n305|65534|0||||cpu=1,mem=16000M,node=1|
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            output_dir='temp',
            job_ident='1672457',
            task_list=3,
            stage='validate'
        )
        _scheduler.SlurmScheduler().update_info(job)
        self.assertEqual(_constants.JOB_STATUS.COMPLETED, job.status)
        self.assertEqual(3, len(job.task_list))


class TestParseScontrolShow(unittest.TestCase):

    def test_pending_job(self):
        content = """
JobId=1673292 JobName=MP_batch-8PyNX8EN4cBdD9vQd9FrRG
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=31 Nice=0 Account=all QOS=normal
   JobState=PENDING Reason=DependencyNeverSatisfied Dependency=afterok:1673291_*
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=16:00:00 TimeMin=N/A
   SubmitTime=2018-05-24T11:32:44 EligibleTime=Unknown
   StartTime=Unknown EndTime=Unknown Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=short AllocNode:Sid=n104:47409
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=(null)
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=16000,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=16000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/output_slurm/pairing/submit.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/output_slurm/pairing/job-%x-1673292.log
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/output_slurm/pairing/job-%x-1673292.log
   Power=

        """
        rows = _scheduler.SlurmScheduler().parse_scontrol_show(content)
        self.assertEqual(1, len(rows))
        self.assertEqual({
            'job_ident': '1673292',
            'task_ident': None,
            'status': 'PENDING',
            'status_comment': 'DependencyNeverSatisfied',
            'name': 'MP_batch-8PyNX8EN4cBdD9vQd9FrRG'
        }, rows[0])

    def test_job_array(self):
        content = """
JobId=1673301 ArrayJobId=1673301 ArrayTaskId=3 JobName=subtest.sh
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=31 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:11 TimeLimit=90-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-24T11:38:28 EligibleTime=2018-05-24T11:38:28
   StartTime=2018-05-24T11:38:29 EndTime=2018-08-22T11:38:29 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=all AllocNode:Sid=n104:47409
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n245
   BatchHost=n245
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_3.out
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_3.out
   Power=

JobId=1673303 ArrayJobId=1673301 ArrayTaskId=2 JobName=subtest.sh
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=31 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:11 TimeLimit=90-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-24T11:38:28 EligibleTime=2018-05-24T11:38:28
   StartTime=2018-05-24T11:38:29 EndTime=2018-08-22T11:38:29 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=all AllocNode:Sid=n104:47409
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n235
   BatchHost=n235
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_2.out
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_2.out
   Power=

JobId=1673302 ArrayJobId=1673301 ArrayTaskId=1 JobName=subtest.sh
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=31 Nice=0 Account=all QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:11 TimeLimit=90-00:00:00 TimeMin=N/A
   SubmitTime=2018-05-24T11:38:28 EligibleTime=2018-05-24T11:38:28
   StartTime=2018-05-24T11:38:29 EndTime=2018-08-22T11:38:29 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=all AllocNode:Sid=n104:47409
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n137
   BatchHost=n137
   NumNodes=1 NumCPUs=1 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=7000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=7000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/subtest.sh
   WorkDir=/projects/trans_scratch/validations/workspace/creisle/temp
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_1.out
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/slurm-1673301_1.out
   Power=

        """
        rows = _scheduler.SlurmScheduler().parse_scontrol_show(content)
        self.assertEqual(3, len(rows))

    def test_cancelled_task(self):
        content = """

JobId=1697512 ArrayJobId=1697503 ArrayTaskId=1 JobName=MV_mock-A47933_batch-uwSwW68EW43XNdvq85NxJ7
   UserId=creisle(1365) GroupId=users(100) MCS_label=N/A
   Priority=42 Nice=0 Account=all QOS=normal
   JobState=CANCELLED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:15
   RunTime=00:00:02 TimeLimit=16:00:00 TimeMin=N/A
   SubmitTime=2018-05-31T20:01:46 EligibleTime=2018-05-31T20:01:49
   StartTime=2018-05-31T20:02:05 EndTime=2018-05-31T20:02:07 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=all AllocNode:Sid=n104:173998
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=n245
   BatchHost=n245
   NumNodes=1 NumCPUs=1 NumTasks=0 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=1,mem=18000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=18000M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/slurm/mock-A47933_diseased_transcriptome/validate/submit.sh
   WorkDir=/home/creisle
   StdErr=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/slurm/mock-A47933_diseased_transcriptome/validate/batch-uwSwW68EW43XNdvq85NxJ7-1/job-%x-1697503-1.log
   StdIn=/dev/null
   StdOut=/projects/trans_scratch/validations/workspace/creisle/temp/test_submission/slurm/mock-A47933_diseased_transcriptome/validate/batch-uwSwW68EW43XNdvq85NxJ7-1/job-%x-1697503-1.log
   Power=

        """
        rows = _scheduler.SlurmScheduler().parse_scontrol_show(content)
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, row['status'])


class TestParseSacctTable(unittest.TestCase):

    def test_basic_table(self):
        content = """
JobID|JobIDRaw|JobName|Partition|MaxVMSize|MaxVMSizeNode|MaxVMSizeTask|AveVMSize|MaxRSS|MaxRSSNode|MaxRSSTask|AveRSS|MaxPages|MaxPagesNode|MaxPagesTask|AvePages|MinCPU|MinCPUNode|MinCPUTask|AveCPU|NTasks|AllocCPUS|Elapsed|State|ExitCode|AveCPUFreq|ReqCPUFreqMin|ReqCPUFreqMax|ReqCPUFreqGov|ReqMem|ConsumedEnergy|MaxDiskRead|MaxDiskReadNode|MaxDiskReadTask|AveDiskRead|MaxDiskWrite|MaxDiskWriteNode|MaxDiskWriteTask|AveDiskWrite|AllocGRES|ReqGRES|ReqTRES|AllocTRES
1672273|1672273|MS_batch-iJUMYRdLFDsuu9eVzGmmKm|short||||||||||||||||||1|00:00:00|CANCELLED by 1365|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|
        """
        rows = _scheduler.SlurmScheduler().parse_sacct(content)
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, row['status'])
    # TODO: test empty header

    def test_cancelled_task(self):
        content = """
JobID|JobName|User|ReqMem|Elapsed|State|MaxRSS|AveRSS|Partition
1697503_3|MV_mock-A47933_batch-uwSwW68EW43XNdvq85NxJ7|creisle|18000Mn|00:00:10|COMPLETED|||all
1697503_3.batch|batch||18000Mn|00:00:10|COMPLETED|904K|904K|
1697503_1|MV_mock-A47933_batch-uwSwW68EW43XNdvq85NxJ7|creisle|18000Mn|00:00:02|CANCELLED by 1365|||all
1697503_1.batch|batch||18000Mn|00:00:02|CANCELLED|896K|896K|
1697503_2|MV_mock-A47933_batch-uwSwW68EW43XNdvq85NxJ7|creisle|18000Mn|00:00:10|COMPLETED|||all
1697503_2.batch|batch||18000Mn|00:00:10|COMPLETED|904K|904K|
        """
        rows = _scheduler.SlurmScheduler().parse_sacct(content)
        self.assertEqual(3, len(rows))
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, rows[1]['status'])
        self.assertEqual(_constants.JOB_STATUS.COMPLETED, rows[0]['status'])

    def test_pending_array(self):
        content = """
JobID|JobName|User|ReqMem|Elapsed|State|MaxRSS|AveRSS|Partition
1701003_[37-200]|MA_L1522785992-normal_batch-aUmErftiY7eEWvENfSeJwc|creisle|12000Mn|00:00:00|PENDING|||all
1701003_1|MA_L1522785992-normal_batch-aUmErftiY7eEWvENfSeJwc|creisle|12000Mn|00:05:00|RUNNING|||all
        """
        rows = _scheduler.SlurmScheduler().parse_sacct(content)
        self.assertEqual(2, len(rows))
        self.assertEqual(_constants.JOB_STATUS.PENDING, rows[0]['status'])
        self.assertEqual(_constants.JOB_STATUS.RUNNING, rows[1]['status'])
        self.assertIs(None, rows[0]['task_ident'])
        self.assertEqual(1, rows[1]['task_ident'])

    def test_resubmission_array(self):
        content = """
JobID|JobIDRaw|JobName|Partition|MaxVMSize|MaxVMSizeNode|MaxVMSizeTask|AveVMSize|MaxRSS|MaxRSSNode|MaxRSSTask|AveRSS|MaxPages|MaxPagesNode|MaxPagesTask|AvePages|MinCPU|MinCPUNode|MinCPUTask|AveCPU|NTasks|AllocCPUS|Elapsed|State|ExitCode|AveCPUFreq|ReqCPUFreqMin|ReqCPUFreqMax|ReqCPUFreqGov|ReqMem|ConsumedEnergy|MaxDiskRead|MaxDiskReadNode|MaxDiskReadTask|AveDiskRead|MaxDiskWrite|MaxDiskWriteNode|MaxDiskWriteTask|AveDiskWrite|AllocGRES|ReqGRES|ReqTRES|AllocTRES
1873472_162|1873671|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|10:18:26|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_162.batch|1873671.batch|batch||13703984K|n106|0|8708976K|11725204K|n106|0|6743424K|53K|n106|0|53K|10:06:31|n106|0|10:06:31|1|1|10:18:26|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|11767.09M|n106|0|11767.09M|29.74M|n106|0|29.74M||||cpu=1,mem=16000M,node=1
1873472_163|1873672|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|08:09:50|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_163.batch|1873672.batch|batch||13690948K|n106|0|8686468K|11712556K|n106|0|6721328K|45K|n106|0|45K|07:57:40|n106|0|07:57:40|1|1|08:09:50|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|13345.62M|n106|0|13345.62M|26.77M|n106|0|26.77M||||cpu=1,mem=16000M,node=1
1873472_164|1873673|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|12:26:33|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_164.batch|1873673.batch|batch||13730588K|n106|0|9577424K|11750552K|n106|0|6777084K|55K|n106|0|55K|12:13:52|n106|0|12:13:52|1|1|12:26:33|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|17065.30M|n106|0|17065.30M|34.39M|n106|0|34.39M||||cpu=1,mem=16000M,node=1
1873472_165|1873674|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|05:32:32|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_165.batch|1873674.batch|batch||13735224K|n106|0|9574752K|11756988K|n106|0|6773916K|52K|n106|0|52K|05:21:46|n106|0|05:21:46|1|1|05:32:32|COMPLETED|0:0|2.18M|0|0|0|16000Mn|0|15997.17M|n106|0|15997.17M|37.74M|n106|0|37.74M||||cpu=1,mem=16000M,node=1
1873472_166|1873675|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|07:30:37|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_166.batch|1873675.batch|batch||13722476K|n106|0|8669768K|11742400K|n106|0|6702776K|53K|n106|0|53K|07:18:31|n106|0|07:18:31|1|1|07:30:37|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|14716.82M|n106|0|14716.82M|21.39M|n106|0|21.39M||||cpu=1,mem=16000M,node=1
1873472_167|1873676|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|06:45:32|COMPLETED|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_167.batch|1873676.batch|batch||13686828K|n106|0|8596932K|11707132K|n106|0|6565180K|49K|n106|0|49K|06:35:26|n106|0|06:35:26|1|1|06:45:32|COMPLETED|0:0|2.19M|0|0|0|16000Mn|0|10274.82M|n106|0|10274.82M|39.37M|n106|0|39.37M||||cpu=1,mem=16000M,node=1
1873472_168|1873677|MV_P02300_batch-egprnnYFaJtPtnECYfGiKf|all||||||||||||||||||1|16:00:06|TIMEOUT|0:0||Unknown|Unknown|Unknown|16000Mn||||||||||||cpu=1,mem=16000M,node=1|cpu=1,mem=16000M,node=1
1873472_168.batch|1873677.batch|batch||13749848K|n106|0|8700272K|11771032K|n106|0|6734652K|46K|n106|0|46K|15:48:39|n106|0|15:48:39|1|1|16:00:07|CANCELLED|0:15|2.19M|0|0|0|16000Mn|0|10613.36M|n106|0|10613.36M|25.00M|n106|0|25.00M||||cpu=1,mem=16000M,node=1
        """
        rows = _scheduler.SlurmScheduler().parse_sacct(content)
        complete = [row['status'] for row in rows if row['status'] == _constants.JOB_STATUS.COMPLETED]
        fail = [row['status'] for row in rows if row['status'] == _constants.JOB_STATUS.CANCELLED]
        self.assertEqual(6, len(complete))
        self.assertEqual(1, len(fail))


class TestCancel(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.SlurmScheduler.command')
    def test_single_job(self, patcher):
        sched = _scheduler.SlurmScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        patcher.assert_called_with(['scancel', '1234'])

    @mock.patch('mavis.schedule.scheduler.SlurmScheduler.command')
    def test_array_job(self, patcher):
        sched = _scheduler.SlurmScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        for task in job.task_list:
            self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
        patcher.assert_called_with(['scancel', '1234'])

    @mock.patch('mavis.schedule.scheduler.SlurmScheduler.command')
    def test_array_job_task(self, patcher):
        sched = _scheduler.SlurmScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job, task_ident=4)
        self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, job.status)
        for i, task in enumerate(job.task_list):
            if i == 3:
                self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
            else:
                self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, task.status)
        patcher.assert_called_with(['scancel', '1234_4'])

    @mock.patch('mavis.schedule.scheduler.SlurmScheduler.command')
    def test_bad_command(self, patcher):
        patcher.side_effect = [subprocess.CalledProcessError(1, 'cmd')]
        sched = _scheduler.SlurmScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        with self.assertRaises(subprocess.CalledProcessError):
            sched.cancel(job)
        patcher.assert_called_with(['scancel', '1234'])
