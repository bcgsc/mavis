import os
import unittest

from mavis.submit import SCHEDULER_CONFIG, SubmissionScript, OPTIONS


class TestConstructor(unittest.TestCase):

    def test_memory_limit(self):
        script = SubmissionScript('', scheduler='SGE', memory_limit=16)
        self.assertEqual(16, script.memory_limit)

    def test_jobname(self):
        script = SubmissionScript('', scheduler='SGE', jobname='name')
        self.assertEqual('name', script.jobname)

    def test_queue(self):
        script = SubmissionScript('', scheduler='SGE', queue='queue')
        self.assertEqual('queue', script.queue)

    def test_unexpected_argument_error(self):
        with self.assertRaises(TypeError):
            SubmissionScript('', scheduler='SGE', blargh='thing')

    def test_stdout(self):
        script = SubmissionScript('', scheduler='SGE', stdout='stdout')
        self.assertEqual('stdout', script.stdout)

    def test_default_weak_override(self):
        script = SubmissionScript('', scheduler='SGE')
        self.assertEqual(OPTIONS.memory_limit, script.memory_limit)
        os.environ['MAVIS_MEMORY_LIMIT'] = '0'
        script = SubmissionScript('', scheduler='SGE')
        self.assertEqual(0, script.memory_limit)


class TestBuildHeader(unittest.TestCase):

    def test_no_options(self):
        script = SubmissionScript('', scheduler='SGE')
        header = script.build_header()
        exp = SCHEDULER_CONFIG.SGE.shebang
        self.assertTrue(exp in header)

    def test_memory_limit(self):
        script = SubmissionScript('', scheduler='SGE', memory_limit=6000)
        header = script.build_header()
        exp = '#$ -l mem_free={0}G,mem_token={0}G,h_vmem={0}G'.format(6)
        self.assertTrue(exp in header)

    def test_stdout(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing')
        header = script.build_header()
        exp = '#$ -o thing/sge-$JOB_NAME-$JOB_ID.log'
        self.assertTrue(exp in header)

    def test_import_env(self):
        script = SubmissionScript('', scheduler='SGE', import_env=True)
        header = script.build_header()
        exp = '#$ -V'
        self.assertTrue(exp in header)
