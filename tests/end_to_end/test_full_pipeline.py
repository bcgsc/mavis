import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from mavis.constants import SUBCOMMAND, EXIT_OK, EXIT_ERROR
from mavis.main import main
from mavis.util import unique_exists

from . import glob_exists
from ..util import get_data


CONFIG = get_data('pipeline_config.cfg')
BWA_CONFIG = get_data('bwa_pipeline_config.cfg')
CLEAN_CONFIG = get_data('clean_pipeline_config.cfg')
MOCK_GENOME = 'mock-A36971'
MOCK_TRANS = 'mock-A47933'
ENV = {e: v for e, v in os.environ.items() if not e.startswith('MAVIS_')}
ENV.update({'MAVIS_SCHEDULER': 'LOCAL'})

@unittest.skipIf(
    not int(os.environ.get('RUN_FULL', 1)),
    'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestPipeline(unittest.TestCase):

    def setUp(self):
        # create the temp output directory to store file outputs
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

    def print_file_tree(self):
        for root, dirs, files in os.walk(self.temp_output):
            level = root.replace(self.temp_output, '').count(os.sep)
            indent = ' ' * 4 * (level)
            print('{}{}/'.format(indent, os.path.basename(root)))
            subindent = ' ' * 4 * (level + 1)
            for f in files:
                print('{}{}'.format(subindent, f))

    def check_annotate(self, lib):
        # run annotation
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.ANNOTATE))
        # check the generated files
        for filename in ['annotations.tab', 'annotations.fusion-cdna.fa', 'drawings', 'drawings/*svg', 'drawings/*json', 'MAVIS.COMPLETE']:
            filename = os.path.join(self.temp_output, lib, SUBCOMMAND.ANNOTATE, '*-1', filename)
            self.assertTrue(glob_exists(filename), msg=filename)

    def check_validate(self, lib):
        # run validation
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))

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
            'validation-passed.tab',
            'MAVIS.COMPLETE'
        ]:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', suffix), msg=suffix)

    def check_aligner_output_files(self, lib, mem=False):
        if mem:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', 'contigs.bwa_mem.sam'))
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', 'contigs.bwa_mem.log'))
        else:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', 'contigs.blat_out.pslx'))

    def check_cluster(self, lib, skipped=False):
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER))
        logfile = os.path.join(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'MC_{}*batch-*.log'.format(lib))
        self.assertTrue(glob_exists(logfile), msg=logfile)
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'batch-*-1.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'filtered_pairs.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'clusters.bed'))
        if skipped:
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
        else:
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'cluster_assignment.tab'))
        self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.CLUSTER, 'MAVIS.COMPLETE'))

    def check_pairing(self):
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.PAIR, 'MAVIS.COMPLETE'))

    def check_summary(self, count=3):
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'mavis_summary*.tab', n=count))
        self.assertTrue(glob_exists(self.temp_output, SUBCOMMAND.SUMMARY, 'MAVIS.COMPLETE'))

    @mock.patch('os.environ', ENV.copy())
    def test_pipeline_with_bwa(self):
        main([SUBCOMMAND.PIPELINE, BWA_CONFIG, '-o', self.temp_output])
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])
        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_validate(lib)
            self.check_aligner_output_files(lib, mem=True)
            self.check_annotate(lib)

        self.check_pairing()
        self.check_summary()

    @mock.patch('os.environ', ENV.copy())
    def test_error_on_bad_config(self):
        with self.assertRaises(SystemExit) as err:
            main([SUBCOMMAND.PIPELINE, 'thing/that/doesnot/exist.cfg', '-o', self.temp_output])
        self.assertEqual(2, err.exception.code)

    @mock.patch('os.environ', ENV.copy())
    def test_error_on_bad_input_file(self):
        with self.assertRaises(FileNotFoundError):
            main([SUBCOMMAND.PIPELINE, get_data('bad_input_file.cfg'), '-o', self.temp_output])

    @mock.patch('os.environ', ENV.copy())
    def test_missing_reference(self):
        with self.assertRaises(OSError):
            main([SUBCOMMAND.PIPELINE, get_data('missing_reference.cfg'), '-o', self.temp_output])

    @mock.patch('os.environ', ENV.copy())
    def test_full_pipeline(self):
        main([SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output])
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))

        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])
        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_annotate(lib)

        self.check_pairing()
        self.check_summary()

        retcode = main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output])
        self.assertEqual(EXIT_OK, retcode)

    @mock.patch('os.environ', ENV.copy())
    def test_no_optional_files(self):
        main([SUBCOMMAND.PIPELINE, get_data('no_opt_pipeline.cfg'), '-o', self.temp_output])
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))

        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])
        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.check_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary()

    @mock.patch('os.environ', ENV.copy())
    def test_reference_from_env(self):
        os.environ.update({
            'MAVIS_TEMPLATE_METADATA': get_data('cytoBand.txt'),
            'MAVIS_ANNOTATIONS': get_data('mock_annotations.json'),
            'MAVIS_MASKING': get_data('mock_masking.tab'),
            'MAVIS_REFERENCE_GENOME': get_data('mock_reference_genome.fa'),
            'MAVIS_ALIGNER_REFERENCE': get_data('mock_reference_genome.2bit'),
            'MAVIS_DGV_ANNOTATION': get_data('mock_dgv_annotation.txt'),
        })
        main([SUBCOMMAND.PIPELINE, get_data('reference_from_env.cfg'), '-o', self.temp_output])
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])
        # check that the subdirectories were built
        for lib in [MOCK_GENOME + '_*']:
            self.check_cluster(lib)
            self.check_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary(count=2)

    @mock.patch('os.environ', ENV.copy())
    def test_clean_files(self):
        main([SUBCOMMAND.PIPELINE, CLEAN_CONFIG, '-o', self.temp_output])

        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])

        # check that the subdirectories were built
        for lib in [MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))

            for suffix in [
                'evidence.bed',
                'validation-failed.tab',
                'validation-passed.tab'
            ]:
                self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', suffix))
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
                self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', suffix), msg=suffix)
            self.assertTrue(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE + '/*-1', 'MAVIS.COMPLETE'))

            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary()

    @mock.patch('os.environ', ENV.copy())
    def test_skip_clustering(self):
        main([SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.CLUSTER])
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib, skipped=True)
            self.check_validate(lib)
            self.check_aligner_output_files(lib)
            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary()

    @mock.patch('os.environ', ENV.copy())
    def test_skip_validation(self):
        main([SUBCOMMAND.PIPELINE, CONFIG, '-o', self.temp_output, '--skip_stage', SUBCOMMAND.VALIDATE])

        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib)
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary()

    @mock.patch('os.environ', ENV.copy())
    def test_skip_cluster_and_validate(self):
        args = [
            SUBCOMMAND.PIPELINE, CONFIG,
            '-o', self.temp_output,
            '--skip_stage', SUBCOMMAND.VALIDATE,
            '--skip_stage', SUBCOMMAND.CLUSTER
        ]
        main(args)
        self.assertTrue(glob_exists(self.temp_output, 'build.cfg'))
        main([SUBCOMMAND.SCHEDULE, '-o', self.temp_output, '--submit'])

        # check that the subdirectories were built
        for lib in[MOCK_GENOME + '_*', MOCK_TRANS + '_*']:
            self.check_cluster(lib, skipped=True)
            self.assertFalse(glob_exists(self.temp_output, lib, SUBCOMMAND.VALIDATE))
            self.check_annotate(lib)
        self.check_pairing()
        self.check_summary()

    def tearDown(self):
        # remove the temp directory and outputs
        self.print_file_tree()
        shutil.rmtree(self.temp_output)


if __name__ == '__main__':
    unittest.main()
