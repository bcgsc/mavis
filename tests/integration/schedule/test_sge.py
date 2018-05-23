import unittest
from unittest import mock

class TestStatus(unittest.TestCase):
    # TODO: status of array job
    # TODO: status of single job
    # TODO: status of job waiting on dependency

    @mock.patch('subprocess.check_output')
    def test_job_array(self, patch_check):
        content = """
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
==============================================================
job_number:                 217940
exec_file:                  job_scripts/217940
submission_time:            Tue May 22 23:39:55 2018
owner:                      creisle
uid:                        1365
group:                      users
gid:                        100
sge_o_home:                 /home/creisle
sge_o_log_name:             creisle
sge_o_path:                 /home/creisle/bin
sge_o_shell:                /bin/bash
sge_o_workdir:              /projects/trans_scratch/validations/workspace/creisle/temp
sge_o_host:                 m0001
account:                    sge
cwd:                        /projects/trans_scratch/validations/workspace/creisle/temp
mail_list:                  creisle@m0001.hpc.bcgsc.ca
notify:                     FALSE
job_name:                   subtest.sh
jobshare:                   0
env_list:
script_file:                subtest.sh
scheduling info:            queue instance "clingen.q@qh32.hpc.bcgsc.ca" dropped because it is temporarily not available
                            queue instance "clingen.q@qh33.hpc.bcgsc.ca" dropped because it is temporarily not available
                            All queues dropped because of overload or full
        """.encode('utf8')
        patch_check.return_value = content
        job = _job.ArrayJob(
            job_ident='217940',
            tasks=2,
            stage='validate'
        )
        hist = _job.SlurmScheduler.status(job)
        self.assertEqual(1, len(hist))
        self.assertEqual(_constants.JOB_STATUS.RUNNING, list(hist.keys())[0])
        self.assertEqual(1, list(hist.values())[0])