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

from . import glob_exists

DATA_PREFIX = os.path.join(os.path.dirname(__file__), 'data')
CONFIG = os.path.join(DATA_PREFIX, 'pipeline_config.cfg')
BWA_CONFIG = os.path.join(DATA_PREFIX, 'bwa_pipeline_config.cfg')
CLEAN_CONFIG = os.path.join(DATA_PREFIX, 'clean_pipeline_config.cfg')
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
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

    def test_pipeline_with_bwa(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, BWA_CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            print(args)
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

            # run validation
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.VALIDATE, '*-1/submit.sh'))
            args = convert_qsub_to_args(qsub)  # read the arguments from the file
            print(args)
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            for suffix in [
                'contigs.bam',
                'contigs.bwa_mem.sam',
                'contigs.bwa_mem.log',
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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_error_on_bad_config(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, 'thing/that/doesnot/exist.cfg', '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            with self.assertRaises(OSError):
                self.assertEqual(0, main())

    def test_full_pipeline(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

            # run validation
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.VALIDATE, '*-1/submit.sh'))
            args = convert_qsub_to_args(qsub)  # read the arguments from the file
            print(args)
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            for suffix in [
                'contigs.bam',
                'contigs.blat_out.pslx',
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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_clean_files(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CLEAN_CONFIG, '-o', self.temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_skip_clustering(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.CLUSTER]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

            # run validation
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            qsub = unique_exists(os.path.join(self.temp_output, lib, SUBCOMMAND.VALIDATE, '*-1/submit.sh'))
            args = convert_qsub_to_args(qsub)  # read the arguments from the file
            print(args)
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            for suffix in [
                'contigs.bam',
                'contigs.blat_out.pslx',
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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def test_skip_validation(self):
        args = ['mavis', SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.VALIDATE]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

            # validation
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))

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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

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
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, '*.COMPLETE'))

            # validation
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))

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
        # now run the pairing
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.PAIR, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        qsub = unique_exists(os.path.join(self.temp_output, SUBCOMMAND.SUMMARY, 'submit.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=3))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CHECKER, '-o', self.temp_output]):
            self.assertEqual(0, main())
        self.assertTrue(glob_exists(self.temp_output, 'submit_pipeline*.sh'))

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)


if __name__ == '__main__':
    unittest.main()
