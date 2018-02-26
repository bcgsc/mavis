import glob
import os
import shutil
import sys
import tempfile
import unittest
import statistics
# sys.stderr = sys.stdout  # redirect so stderr is captured during testing

from mavis.constants import SUBCOMMAND
from mavis.main import main
from mavis.tools import SUPPORTED_TOOL
from mavis.util import unique_exists
from mock import patch


DATA_PREFIX = os.path.join(os.path.dirname(__file__), './../integration/data')


class TestConvert(unittest.TestCase):

    def setUp(self):
        if 'MAVIS_ANNOTATIONS' in os.environ:
            del os.environ['MAVIS_ANNOTATIONS']
        self.temp_output = tempfile.mkdtemp()
        # [--library <name> {genome,transcriptome} {diseased,normal} [strand_specific] [/path/to/bam/file]]
        self.genome = ['--library', 'mock_genome', 'genome', 'diseased']
        self.genome_bam = os.path.join(DATA_PREFIX, 'mock_reads_for_events.sorted.bam')
        self.trans = ['--library', 'mock_trans', 'transcriptome', 'diseased']
        self.trans_bam = os.path.join(DATA_PREFIX, 'mock_trans_reads_for_events.sorted.bam')
        self.annotations = os.path.join(DATA_PREFIX, 'mock_reference_annotations.json')
        self.args = [
            'mavis', SUBCOMMAND.CONFIG
        ]
        self.input = os.path.join(DATA_PREFIX, 'mock_sv_events.tsv')

    def run_main(self, exit_status=0):
        outputfile = os.path.join(self.temp_output, 'config.cfg')
        self.args.extend(['-w', outputfile])
        with patch.object(sys, 'argv', [str(a) for a in self.args]):
            print('sys.argv', sys.argv)
            try:
                return_code = main()
            except SystemExit as ex:
                return_code = ex.code
            self.assertEqual(exit_status, return_code)

    def test_no_libs_no_annotations(self):
        self.run_main()

    def test_skip_no_annotations(self):
        self.args.extend(self.trans + ['False', self.trans_bam, '--input', self.input, 'mock_trans', '--skip_stage', SUBCOMMAND.VALIDATE])
        self.run_main()

    def test_requires_annotations_trans(self):
        self.args.extend(self.trans + ['False', self.trans_bam, '--input', self.input, 'mock_trans'])
        self.run_main(2)

    def test_require_bam_noskip_error(self):
        self.args.extend(self.genome + ['--annotations', self.annotations, '--input', self.input, 'mock_genome'])
        self.run_main(2)

    def test_genome_only(self):
        # should be ok without the annotations file
        self.args.extend(self.genome + ['False', self.genome_bam, '--input', self.input, 'mock_genome'])
        self.run_main()

    def test_trans_with_annotations(self):
        self.args.extend(
            self.genome +
            [False, self.genome_bam] +
            self.trans +
            [True, self.trans_bam, '--input', self.input, 'mock_genome', 'mock_trans', '--annotations', self.annotations]
        )
        with self.assertRaises(statistics.StatisticsError):  # too few annotations to calc median
            self.run_main()

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)


if __name__ == '__main__':
    unittest.main()
