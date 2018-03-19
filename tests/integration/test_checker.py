import os
import shutil
import tempfile
import unittest
from unittest.mock import mock_open, patch

from mavis import checker
from mavis.constants import COMPLETE_STAMP, SUBCOMMAND

DATA_PREFIX = os.path.join(os.path.dirname(__file__), 'data')
MOCK_GENOME = 'mock-A36971'
MOCK_TRANS = 'mock-A47933'


class TestParseLogFile(unittest.TestCase):

    def mock_log(self, content):
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen), patch('os.path.getmtime') as getmtime:
            getmtime.return_value = 1
            isfile.return_value = True
            return checker.LogDetails('log')

    def test_command_not_found_error(self):
        log = self.mock_log(
            "stty: standard input: Inappropriate ioctl for device\n"
            "/opt/slurm/spool/slurmd/job814329/slurm_script: line 9: mavis: command not found\n"
        )
        self.assertEqual(checker.LOGFILE_STATUS.CRASH, log.status)

    def test_python_index_error(self):
        content = """
Traceback (most recent call last):
    File "/home/creisle/git/mavis/venv/bin/mavis", line 11, in <module>
        load_entry_point('mavis===v0.1.0-220-g3f65e68', 'console_scripts', 'mavis')()
    File "/home/creisle/git/mavis/venv/lib/python3.6/site-packages/mavis-v0.1.0_220_g3f65e68-py3.6.egg/mavis/main.py", line 554, in main
        check_completion(args.output)
    File "/home/creisle/git/mavis/venv/lib/python3.6/site-packages/mavis-v0.1.0_220_g3f65e68-py3.6.egg/mavis/main.py", line 450, in check_completion
        cur_time = check_single_job(d)
    File "/home/creisle/git/mavis/venv/lib/python3.6/site-packages/mavis-v0.1.0_220_g3f65e68-py3.6.egg/mavis/main.py", line 429, in check_single_job
        check_log(max(log_files, key=os.path.getctime))
    File "/home/creisle/git/mavis/venv/lib/python3.6/site-packages/mavis-v0.1.0_220_g3f65e68-py3.6.egg/mavis/main.py", line 359, in check_log
        if 'error' in lines[-1].lower():
IndexError: list index out of range"""
        log = self.mock_log(content)
        self.assertEqual(checker.LOGFILE_STATUS.CRASH, log.status)

    def test_python_keyerror(self):
        content = "KeyError: ('cannot check membership column. column not found in header', 'protocol', {'break2_orientation', 'break1_chromosome', 'break1_orientation', 'tools', 'defuse_cluster_id', 'break1_position_end', 'event_type', 'defuse_split_read_count', 'break2_chromosome', 'break2_position_end', 'stranded', 'defuse_spanning_read_count', 'break2_strand', 'library', 'break1_position_start', 'defuse_probability', 'untemplated_seq', 'opposing_strands', 'break1_strand', 'break2_position_start'})"
        log = self.mock_log(content)
        self.assertEqual(checker.LOGFILE_STATUS.CRASH, log.status)

    def test_empty_log(self):
        log = self.mock_log("")
        self.assertEqual(checker.LOGFILE_STATUS.EMPTY, log.status)
        log = self.mock_log("\n\n")
        self.assertEqual(checker.LOGFILE_STATUS.EMPTY, log.status)

    def test_incomplete_log(self):
        log = self.mock_log("other\n")
        self.assertEqual(checker.LOGFILE_STATUS.INCOMPLETE, log.status)
        log = self.mock_log("thing")
        self.assertEqual(checker.LOGFILE_STATUS.INCOMPLETE, log.status)

    def test_parse_run_time(self):
        content = "[2018-03-06 15:25:46.153560] complete: MAVIS.COMPLETE\nrun time (hh/mm/ss): 0:06:41\nrun time (s): 401"
        log = self.mock_log(content)
        self.assertEqual(checker.LOGFILE_STATUS.COMPLETE, log.status)
        self.assertEqual(401, log.run_time)
        content = "[2018-03-06 15:25:46.153560] complete: MAVIS.COMPLETE\nrun time (hh/mm/ss): 0:06:41\nrun time (s): 1\n"
        log = self.mock_log(content)
        self.assertEqual(checker.LOGFILE_STATUS.COMPLETE, log.status)
        self.assertEqual(1, log.run_time)


class TestModule(unittest.TestCase):

    def test_parse_run_time_none(self):
        content = ""
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen), patch('os.path.getmtime') as getmtime:
            getmtime.return_value = 1
            isfile.return_value = True
            result = checker.parse_run_time('log')
        self.assertEqual(None, result)

    def test_parse_valid_run_time(self):
        content = "[2018-03-06 15:25:46.153560] complete: MAVIS.COMPLETE\nrun time (hh/mm/ss): 0:06:41\nrun time (s): 1\n"
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen), patch('os.path.getmtime') as getmtime:
            getmtime.return_value = 1
            isfile.return_value = True
            result = checker.parse_run_time('log')
        self.assertEqual(1, result)


class TestCompletion(unittest.TestCase):

    def mock_log(self, name):
        log = 'batch-mock-1.log'
        content = "[2018-03-06 15:25:46.153560] complete: MAVIS.COMPLETE\nrun time (hh/mm/ss): 0:06:41\nrun time (s): 1\n"
        stamp = COMPLETE_STAMP
        stamp_content = "run time (hh/mm/ss): 0:07:19\nrun time (s): 439\n"
        os.makedirs(name)
        with open(os.path.join(name, log), 'w') as f:
            f.write(content)
        with open(os.path.join(name, stamp), 'w') as f:
            f.write(stamp_content)

    def setUp(self):
        # create the temp output directory to store file outputs
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

        for lib in [MOCK_GENOME + '_diseased_genome', MOCK_TRANS + '_diseased_transcriptome']:
            os.makedirs(os.path.join(self.temp_output, lib))
            for subdir in [SUBCOMMAND.ANNOTATE, SUBCOMMAND.VALIDATE]:
                self.mock_log(os.path.join(self.temp_output, lib, subdir, 'batch-mock-1'))

            self.mock_log(os.path.join(self.temp_output, lib, SUBCOMMAND.CLUSTER))

        for subdir in [SUBCOMMAND.PAIR, SUBCOMMAND.SUMMARY]:
            self.mock_log(os.path.join(self.temp_output, subdir))

    def test_completion_valid_dir(self):
        result = checker.check_completion(self.temp_output)
        self.assertEqual(True, result)

    def test_completion_invalid_dir(self):
        result = checker.check_completion('')
        self.assertEqual(False, result)

    def test_completion_empty_dir(self):
        temp_output = tempfile.mkdtemp()
        result = checker.check_completion(temp_output)
        self.assertEqual(False, result)
        shutil.rmtree(temp_output)

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)
