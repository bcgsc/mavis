import glob
import os
import shlex
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from mavis.constants import SUBCOMMAND
from mavis.main import main
from mavis.util import unique_exists

from . import glob_exists, data


CONFIG = data('pipeline_config.cfg')
BWA_CONFIG = data('bwa_pipeline_config.cfg')
CLEAN_CONFIG = data('clean_pipeline_config.cfg')
MOCK_GENOME = 'mock-A36971'
MOCK_TRANS = 'mock-A47933'


def convert_qsub_to_args(filename, sub=None):
    with open(filename, 'r') as fh:
        lines = [l.strip() for l in fh.readlines() if not l.startswith('#') and l]
        lines = ' '.join(lines)
        if sub:
            for original, replacement in sub:
                lines.replace(original, replacement)
        lex = shlex.split(lines)
        return [a.strip() for a in lex]


@unittest.skipIf(
    not int(os.environ.get('RUN_FULL', 1)),
    'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestPipeline(unittest.TestCase):

    def setUp(self):
        # create the temp output directory to store file outputs
        envs = [e for e in os.environ.keys() if e.startswith('MAVIS_')]
        for evar in envs:
            del os.environ[evar]
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

    def check_and_run_annotate(self, lib):
        # run annotation
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE))
        qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.ANNOTATE, '*-1/submit.sh'))
        args = convert_qsub_to_args(qsub)
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())
        # check the generated files
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', 'annotations.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', 'annotations.fusion-cdna.fa'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', 'drawings'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', 'drawings', '*svg', strict=False))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', 'drawings', '*json', strict=False))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE + '/*-1', '*.COMPLETE'))

    def check_and_run_validate(self, lib):
        # run validation
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
        qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.VALIDATE, '*-1/submit.sh'))
        args = convert_qsub_to_args(qsub)  # read the arguments from the file
        print(args)
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        for suffix in [
            'contigs.bam',
            'contigs.fa',
            'contigs.sorted.bam',
            'contigs.sorted.bam.bai',
            'evidence.bed',
            'igv.batch',
            'raw_evidence.bam',
            'raw_evidence.sorted.bam',
            'raw_evidence.sorted.bam.bai',
            'validation-failed.tab',
            'validation-passed.tab'
        ]:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.' + suffix))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.COMPLETE'))

    def check_aligner_output_files(self, lib, mem=False):
        if mem:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.contigs.bwa_mem.sam'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.contigs.bwa_mem.log'))
        else:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.contigs.blat_out.pslx'))

    def check_cluster(self, lib, skipped=False):
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
        if skipped:
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
        else:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

    def check_and_run_pairing(self):
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub)
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

    def check_and_run_summary(self, count=3):
        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub)
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=count))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

    def test_pipeline_with_bwa(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, BWA_CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            print(args)
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_and_run_validate(lib)
            self.check_aligner_output_files(lib, mem=True)
            self.check_and_run_annotate(lib)

        self.check_and_run_pairing()
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_error_on_bad_config(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, 'thing/that/doesnot/exist.cfg', '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            print(sys.argv)
            try:
                main()
                raise AssertionError('should have thrown error')
            except SystemExit as err:
                self.assertEqual(2, err.code)

    def test_error_on_bad_input_file(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, data('bad_input_file.cfg'), '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(FileNotFoundError):
                main()

    def test_missing_reference(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, data('missing_reference.cfg'), '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(OSError):
                main()

    def test_full_pipeline(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_and_run_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_and_run_annotate(lib)

        # now run the pairing
        self.check_and_run_pairing()

        # now run the summary
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_no_optional_files(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, data('no_opt_pipeline.cfg'), '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_and_run_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_and_run_annotate(lib)
        # now run the pairing
        self.check_and_run_pairing()

        # now run the summary
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_reference_from_env(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, data('reference_from_env.cfg'), '-o', self.temp_output]
        env = {k: v for k, v in os.environ.items()}
        env.update({
            'MAVIS_TEMPLATE_METADATA': 'tests/integration/data/cytoBand.txt',
            'MAVIS_ANNOTATIONS': 'tests/integration/data/mock_annotations.json',
            'MAVIS_MASKING': 'tests/integration/data/mock_masking.tab',
            'MAVIS_REFERENCE_GENOME': 'tests/integration/data/mock_reference_genome.fa',
            'MAVIS_ALIGNER_REFERENCE': 'tests/integration/data/mock_reference_genome.2bit',
            'MAVIS_DGV_ANNOTATION': 'tests/integration/data/mock_dgv_annotation.txt',
        })
        with patch.object(os, 'environ', env):
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            # check that the subdirectories were built
            for lib in [MOCK_GENOME + '_*']:
                self.check_cluster(lib)
                self.check_and_run_validate(lib)
                self.check_aligner_output_files(lib)
                self.check_and_run_annotate(lib)
            # now run the pairing
            self.check_and_run_pairing()

            # now run the summary
            self.check_and_run_summary(count=2)

            with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
                self.assertEqual(0, main())
            self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_clean_files(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CLEAN_CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)

            # run validation
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.VALIDATE, '*-1/submit.sh'))
            args = convert_qsub_to_args(qsub)  # read the arguments from the file
            print(args)
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            for suffix in [
                'evidence.bed',
                'validation-failed.tab',
                'validation-passed.tab'
            ]:
                self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.' + suffix))
            for suffix in [
                'contigs.bam',
                'contigs.blat_out.pslx',
                'contigs.fa',
                'contigs.sorted.bam',
                'contigs.sorted.bam.bai',
                'igv.batch',
                'raw_evidence.bam',
                'raw_evidence.sorted.bam',
                'raw_evidence.sorted.bam.bai',
            ]:
                self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.' + suffix))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', '*.COMPLETE'))

            self.check_and_run_annotate(lib)
        # now run the pairing
        self.check_and_run_pairing()

        # now run the summary
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_skip_clustering(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.CLUSTER]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib, skipped=True)
            self.check_and_run_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_and_run_annotate(lib)
        self.check_and_run_pairing()
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_skip_validation(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.VALIDATE]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            self.check_and_run_annotate(lib)
        self.check_and_run_pairing()
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_skip_cluster_and_validate(self):
        args = [
            'mavis', SUBCOMMAND.PIPELINE, CONFIG,
            '-o', self.temp_output,
            '--skip_stage', SUBCOMMAND.VALIDATE,
            '--skip_stage', SUBCOMMAND.CLUSTER
        ]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib, skipped=True)
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            self.check_and_run_annotate(lib)
        self.check_and_run_pairing()
        self.check_and_run_summary()

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)


if __name__ == '__main__':
    unittest.main()
