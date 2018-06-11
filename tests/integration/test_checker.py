import errno
import os
import unittest
from unittest.mock import mock_open, patch

from mavis.schedule import job as _job
from mavis.schedule import pipeline as _pipeline

MOCK_GENOME = 'mock-A36971'
MOCK_TRANS = 'mock-A47933'
ERROR_MESSAGE = """Traceback (most recent call last):
  File "/home/dpaulino/gitrepo/mavis/venv/bin/mavis_run.py", line 6, in <module>
    exec(compile(open(__file__).read(), __file__, 'exec'))
  File "/home/dpaulino/gitrepo/mavis/bin/mavis_run.py", line 7, in <module>
    from mavis.annotate import load_reference_genes, load_reference_genome, load_masking_regions, load_templates
  File "/home/dpaulino/gitrepo/mavis/mavis/__init__.py", line 6, in <module>
    __version__ = get_version()
  File "/home/dpaulino/gitrepo/mavis/mavis/util.py", line 32, in get_version
    v = subprocess.check_output('cd {}; git describe'.format(os.path.dirname(__file__)), shell=True)
  File "/projects/tumour_char/analysis_scripts/python/centos06/python-3.6.0/lib/python3.6/subprocess.py", line 336, in check_output
    **kwargs).stdout
  File "/projects/tumour_char/analysis_scripts/python/centos06/python-3.6.0/lib/python3.6/subprocess.py", line 418, in run
    output=stdout, stderr=stderr)
subprocess.CalledProcessError: Command 'cd /home/dpaulino/gitrepo/mavis/mavis; git describe' returned non-zero exit status 127."""


def mkdirs(newdir, mode=0o777):
    """
    make directories and ignores if it already exists.
    """
    try:
        os.makedirs(newdir, mode)
    except OSError as err:
        # Reraise the error unless it's about an already existing directory
        if err.errno != errno.EEXIST or not os.path.isdir(newdir):
            raise err


class TestParseLogFile(unittest.TestCase):

    def mock_log(self, content):
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen):
            isfile.return_value = True
            return _job.LogFile.parse('log')

    def test_command_not_found_error(self):
        log = self.mock_log(
            "stty: standard input: Inappropriate ioctl for device\n"
            "/opt/slurm/spool/slurmd/job814329/slurm_script: line 9: mavis: command not found\n"
        )
        self.assertEqual(_job.LogFile.STATUS.CRASH, log.status)

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
        self.assertEqual(_job.LogFile.STATUS.CRASH, log.status)

    def test_python_keyerror(self):
        content = "KeyError: ('cannot check membership column. column not found in header', 'protocol', {'break2_orientation', 'break1_chromosome', 'break1_orientation', 'tools', 'defuse_cluster_id', 'break1_position_end', 'event_type', 'defuse_split_read_count', 'break2_chromosome', 'break2_position_end', 'stranded', 'defuse_spanning_read_count', 'break2_strand', 'library', 'break1_position_start', 'defuse_probability', 'untemplated_seq', 'opposing_strands', 'break1_strand', 'break2_position_start'})"
        log = self.mock_log(content)
        self.assertEqual(_job.LogFile.STATUS.CRASH, log.status)

    def test_empty_log(self):
        log = self.mock_log("")
        self.assertEqual(_job.LogFile.STATUS.EMPTY, log.status)
        log = self.mock_log("\n\n")
        self.assertEqual(_job.LogFile.STATUS.EMPTY, log.status)

    def test_incomplete_log(self):
        log = self.mock_log("other\n")
        self.assertEqual(_job.LogFile.STATUS.INCOMPLETE, log.status)
        log = self.mock_log("thing")
        self.assertEqual(_job.LogFile.STATUS.INCOMPLETE, log.status)


class TestModule(unittest.TestCase):

    def test_parse_run_time_none(self):
        content = ""
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen), patch('os.path.getmtime') as getmtime:
            getmtime.return_value = 1
            isfile.return_value = True
            result = _pipeline.parse_run_time('log')
        self.assertEqual(-1, result)

    def test_parse_valid_run_time(self):
        content = "[2018-03-06 15:25:46.153560] complete: MAVIS.COMPLETE\nrun time (hh/mm/ss): 0:06:41\nrun time (s): 1\n"
        mockopen = mock_open(read_data=content)
        with patch('builtins.open', mockopen), patch('os.path.isfile') as isfile, patch('__main__.open', mockopen), patch('os.path.getmtime') as getmtime:
            getmtime.return_value = 1
            isfile.return_value = True
            result = _pipeline.parse_run_time('log')
        self.assertEqual(1, result)
