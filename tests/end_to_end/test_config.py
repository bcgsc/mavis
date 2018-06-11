import argparse
import glob
import itertools
import os
import shutil
import statistics
import sys
import tempfile
import unittest
from unittest import mock

from mavis.constants import SUBCOMMAND
from mavis.main import main
from mavis.tools import SUPPORTED_TOOL
from mavis.util import unique_exists

from ..util import get_data


ARGERROR_EXIT_CODE = 2


class TestConfig(unittest.TestCase):

    def setUp(self):
        if 'MAVIS_ANNOTATIONS' in os.environ:
            del os.environ['MAVIS_ANNOTATIONS']
        self.temp_output = tempfile.mkdtemp()
        # [--library <name> {genome,transcriptome} {diseased,normal} [strand_specific] [/path/to/bam/file]]
        self.genome = ['--library', 'mock_genome', 'genome', 'diseased']
        self.genome_bam = get_data('mock_reads_for_events.sorted.bam')
        self.trans = ['--library', 'mock_trans', 'transcriptome', 'diseased']
        self.trans_bam = get_data('mock_trans_reads_for_events.sorted.bam')
        self.annotations = get_data('mock_reference_annotations.json')
        self.args = [
            'mavis', SUBCOMMAND.CONFIG
        ]
        self.input = get_data('mock_sv_events.tsv')

    def run_main(self, exit_status=0):
        outputfile = os.path.join(self.temp_output, 'config.cfg')
        self.args.extend(['-w', outputfile])
        with mock.patch.object(sys, 'argv', [str(a) for a in self.args]):
            print('sys.argv', sys.argv)
            try:
                return_code = main()
            except SystemExit as ex:
                return_code = ex.code
            self.assertEqual(exit_status, return_code)

    def test_no_libs_no_annotations(self):
        self.run_main()

    def test_no_input_error(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_input_missing_library(self):
        self.args.extend(self.genome + ['False', self.genome_bam, '--input', self.input, 'mock_genome', 'bad_genome'])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_assign_missing_library(self):
        self.args.extend(self.genome + ['False', self.genome_bam, '--input', self.input, 'mock_genome', '--assign', 'bad_genome', self.input])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_skip_no_annotations(self):
        self.args.extend(self.trans + ['False', self.trans_bam, '--input', self.input, 'mock_trans', '--skip_stage', SUBCOMMAND.VALIDATE])
        self.run_main()

    def test_requires_annotations_trans(self):
        self.args.extend(self.trans + ['False', self.trans_bam, '--input', self.input, 'mock_trans'])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_require_bam_noskip_error(self):
        self.args.extend(self.genome + ['--annotations', self.annotations, '--input', self.input, 'mock_genome'])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_genome_only(self):
        # should be ok without the annotations file
        self.args.extend(self.genome + ['False', self.genome_bam, '--input', self.input, 'mock_genome'])
        self.run_main()

    def test_genome_include_defaults(self):
        # should be ok without the annotations file
        self.args.extend(self.genome + ['False', self.genome_bam, '--input', self.input, 'mock_genome', '--add_defaults'])
        self.run_main()

    def test_trans_with_annotations(self):
        self.args.extend(
            itertools.chain(
                self.genome,
                [False, self.genome_bam],
                self.trans,
                [True, self.trans_bam, '--input', self.input, 'mock_genome', 'mock_trans', '--annotations', self.annotations])
        )
        with self.assertRaises(statistics.StatisticsError):  # too few annotations to calc median
            self.run_main()

    def test_convert_multiple(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss_events.tab', 'transabyss_indels_output.tab', 'transabyss'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main()

    def test_convert_multiple_strand(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss_events.tab', 'transabyss_indels_output.tab', 'transabyss', 'False'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main()

    def test_convert_quoted(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss_{events,indels_output}.tab', 'transabyss'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main()

    def test_convert_quoted_strand(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss_{events,indels_output}.tab', 'transabyss', 'False'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main()

    def test_convert_argument_error(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss', 'False'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main(ARGERROR_EXIT_CODE)

    def test_convert_argument_error2(self):
        self.args.extend(self.genome + ['False', self.genome_bam])
        self.args.extend(['--convert', 'ta', 'transabyss'])
        self.args.extend(['--assign', 'mock_genome', 'ta'])
        self.run_main(ARGERROR_EXIT_CODE)

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)


if __name__ == '__main__':
    unittest.main()
