import glob
import os
import shutil
import sys
import tempfile
import unittest


from mavis.constants import PIPELINE_STEP
from mavis.main import main
from mavis.tools import SUPPORTED_TOOL
from mavis.util import unique_exists
from mock import patch


DATA_PREFIX = os.path.join(os.path.dirname(__file__), 'data')
TEMP_OUTPUT = None


def setUpModule():
    global TEMP_OUTPUT
    # create the temp output directory to store file outputs
    TEMP_OUTPUT = tempfile.mkdtemp()
    print('output dir', TEMP_OUTPUT)


class TestConvert(unittest.TestCase):

    def run_main(self, inputfile, file_type, strand_specific=False):
        args = [
            'mavis', PIPELINE_STEP.CONVERT,
            '-o', TEMP_OUTPUT,
            '-n', inputfile,
            '--file_type', file_type,
            '--strand_specific', strand_specific
        ]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())
            filename = os.path.join(TEMP_OUTPUT, 'mavis_{}_converted.tab'.format(file_type))
            print('output', filename)
            self.assertTrue(unique_exists(filename))
            return filename

    def test_chimerascan(self):
        self.run_main(os.path.join(DATA_PREFIX, 'chimerascan_output.bedpe'), SUPPORTED_TOOL.CHIMERASCAN, False)

    def test_defuse(self):
        self.run_main(os.path.join(DATA_PREFIX, 'defuse_output.tsv'), SUPPORTED_TOOL.DEFUSE, False)

    def test_delly(self):
        self.run_main(os.path.join(DATA_PREFIX, 'delly_events.vcf'), SUPPORTED_TOOL.DELLY, False)

    def test_manta(self):
        self.run_main(os.path.join(DATA_PREFIX, 'manta_events.vcf'), SUPPORTED_TOOL.MANTA, False)

    def test_pindel(self):
        self.run_main(os.path.join(DATA_PREFIX, 'pindel_events.vcf'), SUPPORTED_TOOL.PINDEL, False)

    def test_transabyss(self):
        self.run_main(os.path.join(DATA_PREFIX, 'transabyss_indels_output.tab'), SUPPORTED_TOOL.TA, False)
        self.run_main(os.path.join(DATA_PREFIX, 'transabyss_events.tab'), SUPPORTED_TOOL.TA, False)


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(TEMP_OUTPUT)


if __name__ == '__main__':
    unittest.main()
