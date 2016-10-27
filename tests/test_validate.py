from structural_variant.breakpoint import Breakpoint
from structural_variant.validate import *
from structural_variant.annotate import load_reference_genome, Gene, Transcript
from structural_variant.constants import *
import pysam
import unittest

HUMAN_REFERENCE_GENOME = None


class MockRead(pysam.AlignedSegment):
    """
    utility class to clean up code for generating a mock read
    for testing validate functions
    """

    def __init__(self, **kwargs):
        pysam.AlignedSegment.__init__(self)
        self.reference_start = kwargs.pop('reference_start')
        self.reference_id = kwargs.pop('reference_id')
        self.query_name = kwargs.pop(
            'query_name', 'read-{0}:{1}'.format(self.reference_id, self.reference_start))
        self.query_sequence = kwargs.pop('query_sequence')
        self.cigar = kwargs.pop('cigar')


def setUpModule():
    global HUMAN_REFERENCE_GENOME
    #HUMAN_REFERENCE_GENOME = load_reference_genome('/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa')


class TestEvidence(unittest.TestCase):

    def test_generate_window_orient_ns(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        w = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(440, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_window_orient_left(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        w = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(440, w[0])
        self.assertEqual(1110, w[1])
        self.assertEqual(671, len(w))

    def test_generate_window_orient_right(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        w = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(890, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_transcriptome_window_before_start(self):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1400, 1500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        w2 = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_after_end(self):
        b = Breakpoint(chr='1', start=1600, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1400, 1500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        w2 = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_exonic_long_exon(self):
        b = Breakpoint(chr='1', start=1200, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 2000), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        w2 = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_exon(self):
        b = Breakpoint(chr='1', start=900, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 2000), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        w2 = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_intron(self):
        b = Breakpoint(chr='1', start=1200, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(b, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        w2 = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_short_exon(self):
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (1400, 1500)])
        ann = {'1': [gene]}
        w = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(Interval(1040, 1808), w)

    def test_generate_transcriptome_window_multiple_transcripts(self):
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = Gene(name='KRAS', strand='+', chr='1')
        Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (1400, 1500)])
        Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (2100, 2200)])
        ann = {'1': [gene]}
        w = Evidence.generate_transcriptome_window(b, ann, read_length=100, median_insert_size=250, call_error=10, mm_isize_error=50)
        self.assertEqual(Interval(1040, 2508), w)

if __name__ == "__main__":
    unittest.main()
