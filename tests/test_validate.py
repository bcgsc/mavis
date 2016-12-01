from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.validate import *
from structural_variant.annotate import load_reference_genome, Gene, Transcript
from structural_variant.constants import ORIENT, READ_PAIR_TYPE
from tests import MockRead
import unittest
from tests import REFERENCE_GENOME_FILE

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestEvidence(unittest.TestCase):

    def test__call_by_flanking_reads(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 100, orient=ORIENT.LEFT),
                Breakpoint('fake', 500, orient=ORIENT.RIGHT),
                opposing_strands=False
            ),
            None, None
        )
        ev.flanking_reads[0].add(MockRead(reference_start=20, reference_end=60, next_reference_start=600))
        ev.flanking_reads[0].add(MockRead(reference_start=40, reference_end=80, next_reference_start=650))
        bpp = Evidence._call_by_flanking_reads(ev, SVTYPE.DEL)


    def test_generate_window_orient_ns(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(440, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_window_orient_left(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(440, w[0])
        self.assertEqual(1110, w[1])
        self.assertEqual(671, len(w))

    def test_generate_window_orient_right(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(890, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_transcriptome_window_before_start(self):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1400, 1500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_after_end(self):
        b = Breakpoint(chr='1', start=1600, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1400, 1500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_exonic_long_exon(self):
        b = Breakpoint(chr='1', start=1200, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 2000), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_exon(self):
        b = Breakpoint(chr='1', start=900, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 2000), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_intron(self):
        b = Breakpoint(chr='1', start=1200, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (2400, 2500)])
        ann = {'1': [gene]}
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_short_exon(self):
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        t1 = Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (1400, 1500)])
        ann = {'1': [gene]}
        w = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1040, 1808), w)

    def test_generate_transcriptome_window_multiple_transcripts(self):
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = Gene('1', 1, 9999, name='KRAS', strand='+')
        Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (1400, 1500)])
        Transcript(gene=gene, cds_start=1, cds_end=100, exons=[(1000, 1100), (1200, 1300), (2100, 2200)])
        ann = {'1': [gene]}
        w = Evidence.generate_transcriptome_window(
            b, ann, read_length=100, median_insert_size=250, call_error=10, stdev_isize=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1040, 2508), w)

    def test_read_pair_type_LR(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.LR, Evidence.read_pair_type(r))

    def test_read_pair_type_LR_reverse(self):
        r = MockRead(
            reference_id=1,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.LR, Evidence.read_pair_type(r))

    def test_read_pair_type_LL(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.LL, Evidence.read_pair_type(r))

    def test_read_pair_type_RR(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.RR, Evidence.read_pair_type(r))

    def test_read_pair_type_RL(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.RL, Evidence.read_pair_type(r))

    def test_read_pair_type_RL_reverse(self):
        r = MockRead(
            reference_id=1,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.RL, Evidence.read_pair_type(r))

class MockEvidence:
    def __init__(self, ref=None):
        self.HUMAN_REFERENCE_GENOME = ref


if __name__ == "__main__":
    unittest.main()
