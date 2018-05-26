import unittest

from mavis.annotate.file_io import load_reference_genome
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CIGAR, ORIENT, reverse_complement, STRAND
from mavis.interval import Interval
from mavis.validate.evidence import TranscriptomeEvidence
from mavis.validate.constants import DEFAULTS
from functools import partial

from . import MockRead, MockObject, get_example_genes
from ..util import get_data

REFERENCE_GENOME = None
REF_CHR = 'fake'


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME[REF_CHR].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestNetSizeTransEGFR(unittest.TestCase):

    def setUp(self):
        self.evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts=set(get_example_genes()['EGFR'].transcripts)
        )
        setattr(self.evidence, '_select_transcripts', lambda *pos: self.evidence.overlapping_transcripts)
        setattr(self.evidence, 'distance', partial(TranscriptomeEvidence.distance, self.evidence))

    def egfr_distance(self, pos1, pos2):
        return TranscriptomeEvidence.distance(self.evidence, pos1, pos2)

    def test_deletion_in_exon(self):
        bpp = BreakpointPair(
            Breakpoint('7', 55238890, orient=ORIENT.LEFT), Breakpoint('7', 55238899, orient=ORIENT.RIGHT),
            untemplated_seq=''
        )
        self.assertEqual(Interval(-8), bpp.net_size(self.egfr_distance))

        bpp = BreakpointPair(
            Breakpoint('7', 55238890, orient=ORIENT.LEFT), Breakpoint('7', 55238899, orient=ORIENT.RIGHT),
            untemplated_seq='GTAC'
        )
        self.assertEqual(Interval(-4), bpp.net_size(self.egfr_distance))

    def test_deletion_across_intron(self):
        # 55240539_55240621  55323947_55324313
        bpp = BreakpointPair(
            Breakpoint('7', 55240610, orient=ORIENT.LEFT), Breakpoint('7', 55323950, orient=ORIENT.RIGHT),
            untemplated_seq='GTAC'
        )
        self.assertEqual(Interval(-10), bpp.net_size(self.egfr_distance))
        # 55210998_55211181 55218987_55219055
        bpp = BreakpointPair(
            Breakpoint('7', 55211180, orient=ORIENT.LEFT), Breakpoint('7', 55218990, orient=ORIENT.RIGHT),
            untemplated_seq=''
        )
        self.assertEqual(Interval(-4 + -135, -4), bpp.net_size(self.egfr_distance))

    def test_insertion_at_exon_start_mixed(self):
        # EXON 15: 55232973-55233130
        # EXON 16: 55238868-55238906
        # EXON 17: 55240676-55240817
        bpp = BreakpointPair(
            Breakpoint('7', 55238867, orient=ORIENT.LEFT), Breakpoint('7', 55238868, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(6), bpp.net_size(self.egfr_distance))

    def test_insertion_at_exon_start(self):
        # 55238868_55238906
        bpp = BreakpointPair(
            Breakpoint('7', 55233130, orient=ORIENT.LEFT), Breakpoint('7', 55238868, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(6), bpp.net_size(self.egfr_distance))

    def test_insertion_at_exon_end_mixed(self):
        # 55238868_55238906
        bpp = BreakpointPair(
            Breakpoint('7', 55238905, orient=ORIENT.LEFT), Breakpoint('7', 55238906, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(6), bpp.net_size(self.egfr_distance))

    def test_insertion_at_exon_end(self):
        # 55238868_55238906
        bpp = BreakpointPair(
            Breakpoint('7', 55238906, orient=ORIENT.LEFT), Breakpoint('7', 55240676, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(6), bpp.net_size(self.egfr_distance))

    def test_insertion_in_intron(self):
        # 55238868_55238906
        bpp = BreakpointPair(
            Breakpoint('7', 5523750, orient=ORIENT.LEFT), Breakpoint('7', 5523751, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(6), bpp.net_size(self.egfr_distance))

    def test_indel_in_intron(self):
        # 55238868_55238906
        bpp = BreakpointPair(
            Breakpoint('7', 5523700, orient=ORIENT.LEFT), Breakpoint('7', 5523751, orient=ORIENT.RIGHT),
            untemplated_seq='TTATCG'
        )
        self.assertEqual(Interval(-44), bpp.net_size(self.egfr_distance))


class TestLt(unittest.TestCase):
    def test_break1(self):
        bpp1 = BreakpointPair(Breakpoint('1', 1, 10, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.LEFT), untemplated_seq='')
        bpp2 = BreakpointPair(Breakpoint('1', 1, 9, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.LEFT), untemplated_seq='')
        self.assertTrue(bpp2 < bpp1)

    def test_useq(self):
        bpp1 = BreakpointPair(Breakpoint('1', 1, 10, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.LEFT), untemplated_seq='')
        bpp2 = BreakpointPair(Breakpoint('1', 1, 10, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.LEFT), untemplated_seq=None)
        self.assertTrue(bpp2 > bpp1)

    def test_break2(self):
        bpp1 = BreakpointPair(Breakpoint('1', 1, 10, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.RIGHT), untemplated_seq='')
        bpp2 = BreakpointPair(Breakpoint('1', 1, 10, orient=ORIENT.LEFT), Breakpoint('2', 1, orient=ORIENT.LEFT), untemplated_seq=None)
        self.assertTrue(bpp2 < bpp1)


class TestBreakpointSequenceHomology(unittest.TestCase):

    def test_left_pos_right_pos(self):
        b1 = Breakpoint(REF_CHR, 157, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1788, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CAATGC', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

        b1 = Breakpoint(REF_CHR, 589, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 704, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('TTAA', 'ATAGC'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_left_pos_left_neg(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.NEG, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_left_neg_left_pos(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.NEG, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_right_pos_right_neg(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_right_neg_right_pos(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_close_del(self):
        # ....TT|TT....
        b1 = Breakpoint(REF_CHR, 1001, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1002, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_close_dup(self):
        # ....GATACATTTCTTCTTGAAAA...
        # -------------<=============
        # ===============>-----------
        # -------------CT-CT--------- first break homology
        # ------------T--T----------- second break homology
        b1 = Breakpoint(REF_CHR, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CT', 'TT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_non_specific_error(self):
        b1 = Breakpoint(REF_CHR, 740, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        with self.assertRaises(AttributeError):
            bpp.breakpoint_sequence_homology(REFERENCE_GENOME)
