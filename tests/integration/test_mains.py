import glob
import os
import re
import shutil
from tempfile import mkdtemp
import unittest
from unittest import mock

from mavis.annotate.file_io import load_reference_genes, load_reference_genome, load_templates, ReferenceFile, load_annotations
from mavis.annotate.main import main as annotate_main
from mavis.cluster.main import main as cluster_main
from mavis.constants import DISEASE_STATUS, PROTOCOL
from mavis.validate.main import main as validate_main
import pysam

from . import RUN_FULL
from ..util import get_data

annotations = None
reference_genome = None
template_metadata = None
trans_bam_fh = None
genome_bam_fh = None
masking = mock.Mock(content={})  # do not mask


def setUpModule():
    global annotations, reference_genome, template_metadata, genome_bam_fh, trans_bam_fh, masking
    print('setup start')
    annotations = ReferenceFile('annotations', get_data('mock_annotations.json'))
    reference_genome = ReferenceFile('reference_genome', get_data('mock_reference_genome.fa'), eager_load=True)
    template_metadata = ReferenceFile('template_metadata', get_data('cytoBand.txt'), eager_load=True)
    genome_bam_fh = pysam.AlignmentFile(get_data('mock_reads_for_events.sorted.bam'))
    trans_bam_fh = pysam.AlignmentFile(get_data('mock_trans_reads_for_events.sorted.bam'))
    print('setup loading is complete')


def tearDownModule():
    trans_bam_fh.close()
    genome_bam_fh.close()


@unittest.skipIf(not RUN_FULL, 'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.output = mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.output)

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_mains(self):
        # test the clustering
        cluster_files = cluster_main(
            [get_data('mock_sv_events.tsv')], self.output, False, 'mock-A36971', PROTOCOL.GENOME, DISEASE_STATUS.DISEASED,
            limit_to_chr=[None], log_args=True,
            masking=masking, cluster_clique_size=15, cluster_radius=20,
            uninformative_filter=True, max_proximity=5000,
            annotations=annotations, min_clusters_per_file=5, max_files=1
        )
        self.assertGreaterEqual(100, len(cluster_files))
        self.assertLessEqual(1, len(cluster_files))
        # next test the validate runs without errors
        validate_main(
            [cluster_files[0]], self.output, genome_bam_fh, False, 'mock-A36971', PROTOCOL.GENOME,
            median_fragment_size=427, stdev_fragment_size=106, read_length=150,
            reference_genome=reference_genome, annotations=annotations, masking=masking,
            aligner_reference=ReferenceFile('aligner_reference', get_data('mock_reference_genome.2bit'))
        )
        for suffix in [
            'validation-passed.tab',
            'validation-failed.tab',
            'raw_evidence.bam',
            'raw_evidence.sorted.bam',
            'raw_evidence.sorted.bam.bai',
            'contigs.sorted.bam',
            'contigs.sorted.bam.bai',
            'contigs.bam',
            'igv.batch'
        ]:
            self.assertTrue(os.path.exists(os.path.join(self.output, suffix)))

        # test the annotation
        annotate_main(
            [os.path.join(self.output, 'validation-passed.tab')], self.output, 'mock-A36971', PROTOCOL.GENOME,
            reference_genome, annotations, template_metadata,
            min_domain_mapping_match=0.95, min_orf_size=300, max_orf_cap=3,
        )
        self.assertTrue(os.path.exists(os.path.join(self.output, 'annotations.tab')))
        self.assertTrue(os.path.exists(os.path.join(self.output, 'annotations.fusion-cdna.fa')))
        drawings_dir = os.path.join(self.output, 'drawings')
        self.assertTrue(os.path.exists(drawings_dir))
        self.assertLessEqual(1, len(glob.glob(os.path.join(drawings_dir, '*.svg'))))
        self.assertLessEqual(1, len(glob.glob(os.path.join(drawings_dir, '*.legend.json'))))
