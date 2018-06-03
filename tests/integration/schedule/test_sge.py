import subprocess
import unittest
from unittest import mock

from mavis.schedule import job as _job
from mavis.schedule import scheduler as _scheduler
from mavis.schedule import constants as _constants
from mavis.constants import SUBCOMMAND

QACCT_ARR3_OK = """
==============================================================
qname        merge.q
hostname     n601.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      arrtest
jobnumber    3757289
taskid       1
account      sge
priority     0
qsub_time    Thu May 24 10:54:05 2018
start_time   Thu May 24 10:54:12 2018
end_time     Thu May 24 10:55:12 2018
granted_pe   NONE
slots        1
failed       0
exit_status  0
ru_wallclock 60s
ru_utime     0.057s
ru_stime     0.087s
ru_maxrss    5.160KB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    20948
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   8
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     224
ru_nivcsw    59
cpu          0.144s
mem          0.000GBs
io           0.001GB
iow          0.000s
maxvmem      1.934MB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users
==============================================================
qname        merge.q
hostname     n602.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      arrtest
jobnumber    3757289
taskid       3
account      sge
priority     0
qsub_time    Thu May 24 10:54:05 2018
start_time   Thu May 24 10:54:12 2018
end_time     Thu May 24 10:55:12 2018
granted_pe   NONE
slots        1
failed       0
exit_status  0
ru_wallclock 60s
ru_utime     0.063s
ru_stime     0.079s
ru_maxrss    5.156KB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    20954
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   8
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     220
ru_nivcsw    65
cpu          0.142s
mem          0.000GBs
io           0.001GB
iow          0.000s
maxvmem      1.934MB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users
==============================================================
qname        merge.q
hostname     n604.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      arrtest
jobnumber    3757289
taskid       2
account      sge
priority     0
qsub_time    Thu May 24 10:54:05 2018
start_time   Thu May 24 10:54:17 2018
end_time     Thu May 24 10:55:17 2018
granted_pe   NONE
slots        1
failed       0
exit_status  0
ru_wallclock 60s
ru_utime     0.055s
ru_stime     0.086s
ru_maxrss    5.156KB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    20954
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   8
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     218
ru_nivcsw    66
cpu          0.141s
mem          0.000GBs
io           0.001GB
iow          0.000s
maxvmem      1.930MB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users
"""


class TestUpdate(unittest.TestCase):
    # TODO: status of array job
    # TODO: status of single job
    # TODO: status of job waiting on dependency

    @mock.patch('subprocess.check_output')
    def test_job_array_waiting(self, patch_check):
        content = """
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
3751935 0.00000 subtest.sh creisle      qw    05/23/2018 13:44:04                                    1 1-10:1
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            output_dir='temp',
            job_ident='3751935',
            task_list=10,
            stage='validate'
        )
        _scheduler.SgeScheduler().update_info(job)
        self.assertEqual(_constants.JOB_STATUS.PENDING, job.status)

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        content = """
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n601.numbers.bcgsc.ca      1 1
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n602.numbers.bcgsc.ca      1 2
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n604.numbers.bcgsc.ca      1 3
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n603.numbers.bcgsc.ca      1 4
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n601.numbers.bcgsc.ca      1 5
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n602.numbers.bcgsc.ca      1 6
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n604.numbers.bcgsc.ca      1 7
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n603.numbers.bcgsc.ca      1 8
3751935 0.50500 subtest.sh creisle      r     05/23/2018 13:44:12 merge.q@n601.numbers.bcgsc.ca      1 9
3751935 0.50500 subtest.sh creisle      qw    05/23/2018 13:44:12 merge.q@n602.numbers.bcgsc.ca      1 10
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            output_dir='temp',
            job_ident='3751935',
            task_list=10,
            stage='validate'
        )
        _scheduler.SgeScheduler().update_info(job)

        for task in job.task_list[:9]:
            self.assertEqual(_constants.JOB_STATUS.RUNNING, task.status)
        self.assertEqual(_constants.JOB_STATUS.PENDING, job.task_list[-1].status)
        self.assertEqual(_constants.JOB_STATUS.PENDING, job.status)

    @mock.patch('subprocess.check_output')
    def test_single_job(self, patch_check):
        content = """
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
 217940 1.50000 subtest.sh creisle      qw    05/22/2018 23:39:55                                    1
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.Job(
            output_dir='temp',
            job_ident='217940',
            stage='validate'
        )
        _scheduler.SgeScheduler().update_info(job)
        self.assertEqual(_constants.JOB_STATUS.PENDING, job.status)

    @mock.patch('subprocess.check_output')
    def test_completed_array(self, patch_check):
        patch_check.side_effect = [
            ''.encode('utf8'),
            QACCT_ARR3_OK.encode('utf8')
        ]
        job = _job.ArrayJob(
            output_dir='temp',
            job_ident='3757289',
            stage='validate',
            name='arrtest',
            task_list=3
        )
        _scheduler.SgeScheduler().update_info(job)
        self.assertEqual(_constants.JOB_STATUS.COMPLETED, job.status)
        for task in job.task_list:
            self.assertEqual(_constants.JOB_STATUS.COMPLETED, task.status)


class TestParseQacct(unittest.TestCase):

    def test_job_array(self):
        content = QACCT_ARR3_OK
        rows = _scheduler.SgeScheduler().parse_qacct(content)
        expected = {
            'job_ident': '3757289',
            'name': 'arrtest',
            'status': _constants.JOB_STATUS.COMPLETED,
            'status_comment': ''
        }
        for task_id, row in zip([1, 3, 2], rows):
            exp = {'task_ident': str(task_id)}
            exp.update(expected)
            self.assertEqual(exp, row)

    def test_passed(self):
        content = """
==============================================================
qname        transabyss.q
hostname     tac3n15.hpc.bcgsc.ca
group        users
owner        bioapps
project      NONE
department   defaultdepartment
jobname      A89009negative
jobnumber    3744253
taskid       40
account      sge
priority     0
qsub_time    Tue May 22 09:26:31 2018
start_time   Tue May 22 10:32:42 2018
end_time     Tue May 22 13:28:32 2018
granted_pe   openmpi
slots        8
failed       0
exit_status  0
ru_wallclock 10550s
ru_utime     42298.581s
ru_stime     34509.422s
ru_maxrss    2.608MB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    5382919
ru_majflt    978
ru_nswap     0
ru_inblock   14027520
ru_oublock   9259368
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     20635137
ru_nivcsw    14100587
cpu          76808.002s
mem          119.207KGBs
io           579.782GB
iow          0.000s
maxvmem      14.885GB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users -q transabyss.q -l h_vmem=3.85G,mem_free=3.85G,mem_token=3.85G -pe openmpi 8
        """
        rows = _scheduler.SgeScheduler().parse_qacct(content)
        self.assertEqual(1, len(rows))
        expected = {
            'job_ident': '3744253',
            'task_ident': '40',
            'name': 'A89009negative',
            'status': _constants.JOB_STATUS.COMPLETED,
            'status_comment': ''
        }
        self.assertEqual(expected, rows[0])

    def test_non_zero_exit(self):
        content = """
==============================================================
qname        merge.q
hostname     n603.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      error
jobnumber    3755560
taskid       undefined
account      sge
priority     0
qsub_time    Thu May 24 09:42:58 2018
start_time   Thu May 24 09:43:12 2018
end_time     Thu May 24 09:44:12 2018
granted_pe   NONE
slots        1
failed       0
exit_status  1
ru_wallclock 60s
ru_utime     0.054s
ru_stime     0.088s
ru_maxrss    5.148KB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    21134
ru_majflt    0
ru_nswap     0
ru_inblock   8
ru_oublock   16
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     228
ru_nivcsw    62
cpu          0.142s
mem          0.000GBs
io           0.001GB
iow          0.000s
maxvmem      1.926MB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users
        """

        rows = _scheduler.SgeScheduler().parse_qacct(content)
        self.assertEqual(1, len(rows))
        expected = {
            'job_ident': '3755560',
            'task_ident': None,
            'name': 'error',
            'status': _constants.JOB_STATUS.FAILED,
            'status_comment': ''
        }
        self.assertEqual(expected, rows[0])

    def test_failed(self):
        content = """
==============================================================
qname        merge.q
hostname     n603.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      MV_mock-A36971_batch-E6aEZJnTQAau598tcsMjAE
jobnumber    3760712
taskid       1
account      sge
priority     0
qsub_time    Thu May 24 13:35:02 2018
start_time   -/-
end_time     -/-
granted_pe   NONE
slots        1
failed       26  : opening input/output file
exit_status  0
ru_wallclock 0s
ru_utime     0.000s
ru_stime     0.000s
ru_maxrss    0.000B
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    0
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   0
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     0
ru_nivcsw    0
cpu          0.000s
mem          0.000GBs
io           0.000GB
iow          0.000s
maxvmem      0.000B
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users -l h_rt=57600,h_vmem=16000M,mem_free=16000M,mem_token=16000M
        """
        rows = _scheduler.SgeScheduler().parse_qacct(content)
        self.assertEqual(1, len(rows))
        expected = {
            'job_ident': '3760712',
            'task_ident': '1',
            'name': 'MV_mock-A36971_batch-E6aEZJnTQAau598tcsMjAE',
            'status': _constants.JOB_STATUS.FAILED,
            'status_comment': 'opening input/output file'
        }
        self.assertEqual(expected, rows[0])

    def test_cancelled(self):
        content = """
==============================================================
qname        merge.q
hostname     n603.numbers.bcgsc.ca
group        users
owner        creisle
project      NONE
department   defaultdepartment
jobname      arrtest
jobnumber    3757249
taskid       undefined
account      sge
priority     0
qsub_time    Thu May 24 10:50:27 2018
start_time   Thu May 24 10:50:45 2018
end_time     Thu May 24 10:51:09 2018
granted_pe   NONE
slots        1
failed       100 : assumedly after job
exit_status  137                  (Killed)
ru_wallclock 24s
ru_utime     0.052s
ru_stime     0.088s
ru_maxrss    5.160KB
ru_ixrss     0.000B
ru_ismrss    0.000B
ru_idrss     0.000B
ru_isrss     0.000B
ru_minflt    20737
ru_majflt    0
ru_nswap     0
ru_inblock   0
ru_oublock   8
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     215
ru_nivcsw    63
cpu          0.140s
mem          0.000GBs
io           0.001GB
iow          0.000s
maxvmem      1.934MB
arid         undefined
ar_sub_time  undefined
category     -U transabyss_users
        """
        rows = _scheduler.SgeScheduler().parse_qacct(content)
        self.assertEqual(1, len(rows))
        expected = {
            'job_ident': '3757249',
            'task_ident': None,
            'name': 'arrtest',
            'status': _constants.JOB_STATUS.CANCELLED,
            'status_comment': 'assumedly after job'
        }
        self.assertEqual(expected, rows[0])

    def test_job_not_found(self):
        content = """
Total System Usage
    WALLCLOCK         UTIME         STIME           CPU             MEMORY                 IO                IOW
================================================================================================================
   3786481073 6713770428.951 4374477378.582 11585461604.347   187237653407.317      156350319.140              0.000
        """
        with self.assertRaises(ValueError):
            _scheduler.SgeScheduler().parse_qacct(content)


class TestParseQstat(unittest.TestCase):
    def test_single_job(self):
        content = """
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
-----------------------------------------------------------------------------------------------------------------
 217940 1.50000 subtest.sh creisle      qw    05/22/2018 23:39:55                                    1
        """
        rows = _scheduler.SgeScheduler().parse_qstat(content)
        self.assertEqual(1, len(rows))
        expected = {
            'job_ident': '217940',
            'task_ident': None,
            'status': _constants.JOB_STATUS.PENDING,
            'name': 'subtest.sh',
            'status_comment': ''
        }
        self.assertEqual(expected, rows[0])

    def test_no_jobs_found(self):
        rows = _scheduler.SgeScheduler().parse_qstat("")
        self.assertEqual([], rows)


class TestCancel(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_single_job(self, patcher):
        sched = _scheduler.SgeScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job(self, patcher):
        sched = _scheduler.SgeScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.CANCELLED, job.status)
        for task in job.task_list:
            self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
        patcher.assert_called_with(['qdel', '1234'])

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job_task(self, patcher):
        sched = _scheduler.SgeScheduler()
        job = _job.ArrayJob(SUBCOMMAND.VALIDATE, 10, output_dir='', job_ident='1234')
        sched.cancel(job, task_ident=4)
        self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, job.status)
        for i, task in enumerate(job.task_list):
            if i == 3:
                self.assertEqual(_constants.JOB_STATUS.CANCELLED, task.status)
            else:
                self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, task.status)
        patcher.assert_called_with(['qdel', '1234', '-t', '4'])

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_bad_command(self, patcher):
        patcher.side_effect = [subprocess.CalledProcessError(1, 'command')]
        sched = _scheduler.SgeScheduler()
        job = _job.Job(SUBCOMMAND.VALIDATE, '', job_ident='1234')
        sched.cancel(job)
        self.assertEqual(_constants.JOB_STATUS.NOT_SUBMITTED, job.status)


class TestSubmit(unittest.TestCase):

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_job(self, patcher):
        patcher.side_effect = ['Your job 3891651 ("MV1") has been submitted']
        job = _job.Job(
            SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            script='script.sh',
            name='MV1',
            memory_limit=1
        )
        sched = _scheduler.SgeScheduler()
        sched.submit(job)
        self.assertEqual('3891651', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        patcher.assert_called_with(
            'qsub -j y -q all -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-N MV1 -o output_dir/job-\\$JOB_NAME-\\$JOB_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_job_with_array_dep(self, patcher):
        patcher.side_effect = ['Your job 3891651 ("MV1") has been submitted']
        job = _job.Job(
            SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            script='script.sh',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL
        )
        dep = _job.ArrayJob(job_ident='1234', task_list=10, output_dir='', stage=SUBCOMMAND.VALIDATE)
        job.dependencies.append(dep)
        sched = _scheduler.SgeScheduler()
        sched.submit(job)
        self.assertEqual('3891651', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        patcher.assert_called_with(
            'qsub -j y -q all -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-hold_jid 1234 -N MV1 -m abes -M me@example.com '
            '-o output_dir/job-\\$JOB_NAME-\\$JOB_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_job_with_job_dep(self, patcher):
        patcher.side_effect = ['Your job 3891651 ("MV1") has been submitted']
        job = _job.Job(
            SUBCOMMAND.VALIDATE,
            queue='all',
            output_dir='output_dir',
            script='script.sh',
            name='MV1',
            memory_limit=1,
            mail_user='me@example.com',
            mail_type=_constants.MAIL_TYPE.ALL
        )
        dep = _job.Job(job_ident='1234', output_dir='', stage=SUBCOMMAND.VALIDATE)
        job.dependencies.append(dep)
        sched = _scheduler.SgeScheduler()
        sched.submit(job)
        self.assertEqual('3891651', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)
        patcher.assert_called_with(
            'qsub -j y -q all -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-hold_jid 1234 -N MV1 -m abes -M me@example.com '
            '-o output_dir/job-\\$JOB_NAME-\\$JOB_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job(self, patcher):
        patcher.side_effect = ['Your job-array 3891657.2-4:1 ("MV1") has been submitted']
        job = _job.ArrayJob(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', script='script.sh', name='MV1', task_list=[2, 3, 4], memory_limit=1)
        sched = _scheduler.SgeScheduler(concurrency_limit=2)
        sched.submit(job)
        self.assertEqual('3891657', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)

        patcher.assert_called_with(
            'qsub -j y -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-N MV1 -t 2-4 -o output_dir/job-\\$JOB_NAME-\\$JOB_ID-\\$TASK_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job_with_job_dep(self, patcher):
        patcher.side_effect = ['Your job-array 3891657.2-4:1 ("MV1") has been submitted']
        job = _job.ArrayJob(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', script='script.sh', name='MV1', task_list=[2, 3, 4], memory_limit=1)
        sched = _scheduler.SgeScheduler(concurrency_limit=2)

        dep = _job.Job(job_ident='1234', output_dir='', stage=SUBCOMMAND.VALIDATE)
        job.dependencies.append(dep)

        sched.submit(job)
        self.assertEqual('3891657', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)

        patcher.assert_called_with(
            'qsub -j y -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-hold_jid 1234 '
            '-N MV1 -t 2-4 -o output_dir/job-\\$JOB_NAME-\\$JOB_ID-\\$TASK_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job_with_array_dep(self, patcher):
        patcher.side_effect = ['Your job-array 3891657.2-4:1 ("MV1") has been submitted']
        job = _job.ArrayJob(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', script='script.sh', name='MV1', task_list=[2, 3, 4], memory_limit=1)
        sched = _scheduler.SgeScheduler(concurrency_limit=2)

        dep = _job.ArrayJob(job_ident='1234', task_list=[2, 3, 4], output_dir='', stage=SUBCOMMAND.VALIDATE)
        job.dependencies.append(dep)

        sched.submit(job)
        self.assertEqual('3891657', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)

        patcher.assert_called_with(
            'qsub -j y -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-hold_jid_ad 1234 '
            '-N MV1 -t 2-4 -o output_dir/job-\\$JOB_NAME-\\$JOB_ID-\\$TASK_ID.log script.sh', shell=True)

    @mock.patch('mavis.schedule.scheduler.SgeScheduler.command')
    def test_array_job_with_diff_array(self, patcher):
        patcher.side_effect = ['Your job-array 3891657.2-4:1 ("MV1") has been submitted']
        job = _job.ArrayJob(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', script='script.sh', name='MV1', task_list=[2, 3, 4], memory_limit=1)
        sched = _scheduler.SgeScheduler(concurrency_limit=2)

        dep = _job.ArrayJob(job_ident='1234', task_list=[2, 3, 4, 5], output_dir='', stage=SUBCOMMAND.VALIDATE)
        job.dependencies.append(dep)

        sched.submit(job)
        self.assertEqual('3891657', job.job_ident)
        self.assertEqual(_constants.JOB_STATUS.SUBMITTED, job.status)

        patcher.assert_called_with(
            'qsub -j y -l mem_free=1M,mem_token=1M,h_vmem=1M -l h_rt=16:00:00 -V '
            '-hold_jid 1234 '
            '-N MV1 -t 2-4 -o output_dir/job-\\$JOB_NAME-\\$JOB_ID-\\$TASK_ID.log script.sh', shell=True)

    def test_array_job_non_consec_error(self):
        job = _job.ArrayJob(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', script='script.sh', name='MV1', task_list=[2, 3, 4, 7], memory_limit=1)
        sched = _scheduler.SgeScheduler(concurrency_limit=2)
        with self.assertRaises(ValueError):
            sched.submit(job)

    def test_already_submitted_error(self):
        job = _job.Job(stage=SUBCOMMAND.VALIDATE, output_dir='output_dir', job_ident='1')
        sched = _scheduler.SgeScheduler(concurrency_limit=2)
        with self.assertRaises(ValueError):
            sched.submit(job)
