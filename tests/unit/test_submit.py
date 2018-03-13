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

    def test_sge_mail_type_no_user(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='ALL')
        header = script.build_header()
        exp = '#$ -m abes'
        self.assertTrue(exp not in header)

    def test_slurm_mail_type_no_user(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='ALL')
        header = script.build_header()
        exp = '#SBATCH --mail-type=ALL'
        self.assertTrue(exp not in header)

    def test_sge_mail_user_only(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_user='someone')
        header = script.build_header()
        exp = '#$ -M someone'
        self.assertTrue(exp in header)

    def test_slurm_mail_user_only(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_user='someone')
        header = script.build_header()
        exp = '#SBATCH --mail-user=someone'
        self.assertTrue(exp in header)

    def test_sge_mail(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='ALL', mail_user='someone')
        header = script.build_header()
        print(header)
        self.assertTrue('#$ -M someone' in header)
        self.assertTrue('#$ -m abes' in header)

    def test_slurm_mail(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='ALL', mail_user='someone')
        header = script.build_header()
        self.assertTrue('#SBATCH --mail-user=someone' in header)
        self.assertTrue('#SBATCH --mail-type=ALL' in header)

    def test_sge_mail_fail(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='FAIL', mail_user='someone')
        header = script.build_header()
        print(header)
        self.assertTrue('#$ -M someone' in header)
        self.assertTrue('#$ -m as' in header)

    def test_slurm_mail_fail(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='FAIL', mail_user='someone')
        header = script.build_header()
        self.assertTrue('#SBATCH --mail-user=someone' in header)
        self.assertTrue('#SBATCH --mail-type=FAIL' in header)

    def test_sge_mail_option_from_env(self):
        os.environ['MAVIS_MAIL_TYPE'] = 'ALL'
        os.environ['MAVIS_MAIL_USER'] = 'someone'
        script = SubmissionScript('', scheduler='SGE')
        header = script.build_header()
        self.assertTrue('#$ -M someone' in header)
        self.assertTrue('#$ -m abes' in header)

    def test_slurm_mail_option_from_env(self):
        os.environ['MAVIS_MAIL_TYPE'] = 'ALL'
        os.environ['MAVIS_MAIL_USER'] = 'someone'
        script = SubmissionScript('', scheduler='SLURM')
        header = script.build_header()
        self.assertTrue('#SBATCH --mail-user=someone' in header)
        self.assertTrue('#SBATCH --mail-type=ALL' in header)

    def test_slurm_multiple_mail_types(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='FAIL,ALL', mail_user='someone')
        header = script.build_header()
        self.assertTrue('#SBATCH --mail-user=someone' in header)
        self.assertTrue('#SBATCH --mail-type=ALL,FAIL' in header)

    def test_slurm_bad_mail_type(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='FAIL,BAD', mail_user='someone')
        with self.assertRaises(KeyError):
            header = script.build_header()

        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='BAD', mail_user='someone')
        with self.assertRaises(KeyError):
            header = script.build_header()

    def test_sge_bad_mail_type(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='BAD', mail_user='someone')
        with self.assertRaises(KeyError):
            header = script.build_header()

    def test_sge_multiple_mail_options(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='FAIL,ALL', mail_user='someone')
        header = script.build_header()
        self.assertTrue('#$ -M someone' in header)
        self.assertTrue('#$ -m abes' in header)

    def test_slurm_mail_type_none_mix(self):
        script = SubmissionScript('', scheduler='SLURM', stdout='thing', mail_type='ALL,NONE', mail_user='someone')
        with self.assertRaises(ValueError):
            header = script.build_header()

    def test_sge_mail_type_none_mix(self):
        script = SubmissionScript('', scheduler='SGE', stdout='thing', mail_type='ALL,NONE', mail_user='someone')
        with self.assertRaises(ValueError):
            header = script.build_header()

    def tearDown(self):
        try:
            del os.environ['MAVIS_MAIL_TYPE']
        except KeyError:
            pass
        try:
            del os.environ['MAVIS_MAIL_USER']
        except KeyError:
            pass
