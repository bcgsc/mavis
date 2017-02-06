import unittest
from structural_variant.constants import *
from structural_variant.breakpoint import *
from structural_variant.error import *
from structural_variant.annotate import load_reference_genome
from tests import MockRead
from tests import REFERENCE_GENOME_FILE

REFERENCE_GENOME = None
REF_CHR = None


def setUpModule():
    global REFERENCE_GENOME, REF_CHR
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    REF_CHR = list(REFERENCE_GENOME.keys())[0]
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME[REF_CHR].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestBreakpoint(unittest.TestCase):

    def test___eq__(self):
        self.assertNotEqual(Breakpoint('1', 1), None)
        self.assertEqual(Breakpoint('1', 1), Breakpoint('1', 1))

    def test___hash__(self):
        b = Breakpoint('1', 1, 2)
        c = Breakpoint('1', 1, 2)
        d = Breakpoint('1', 1, 1)

        temp = set()
        temp.add(b)
        temp.add(c)
        temp.add(d)
        self.assertEqual(2, len(temp))

        temp = dict()
        temp[b] = None
        temp[c] = None
        temp[d] = None
        self.assertEqual(2, len(temp.keys()))

    def test___len__(self):
        with self.assertRaises(AttributeError):
            Breakpoint('11', 87042760, 87041922, orient=ORIENT.LEFT, strand=STRAND.NS)

    def test_inherited_interval_methods(self):
        b = Breakpoint('1', 1, 10)
        self.assertEqual(1, b[0])
        self.assertEqual(10, b[1])
        self.assertEqual(10, len(b))

    def test_breakpoint_constructor(self):
        b = Breakpoint('1', 10, 50)
        self.assertEqual(10, b[0])
        self.assertEqual(50, b[1])
        self.assertTrue(Interval.overlaps((1, 10), b))
        self.assertTrue(Interval.overlaps((50, 55), b))
        self.assertFalse(Interval.overlaps((1, 9), b))


class TestBreakpointPair(unittest.TestCase):

    def test___eq__(self):
        b = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        c = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        self.assertFalse(b is c)
        self.assertEqual(b, c)
        d = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_sequence='')
        self.assertNotEqual(b, d)
        self.assertNotEqual(b, None)

    def test___hash__(self):
        b = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        c = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        d = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_sequence='')
        self.assertFalse(b is c)
        temp = dict()
        temp[b] = None
        temp[d] = None
        temp[c] = None
        self.assertEqual(2, len(temp.keys()))

        temp = set()
        temp.add(b)
        temp.add(c)
        temp.add(d)
        self.assertEqual(2, len(temp))

    def test___init__swap_break_order(self):
        b1 = Breakpoint('1', 1)
        b2 = Breakpoint('1', 50)
        bpp = BreakpointPair(b1, b2, opposing_strands=True)
        self.assertEqual(bpp.break1, b1)
        self.assertEqual(bpp.break2, b2)
        bpp = BreakpointPair(b2, b1, opposing_strands=True)
        self.assertEqual(bpp.break1, b1)
        self.assertEqual(bpp.break2, b2)

    def test___init__opstrand_conflict(self):
        with self.assertRaises(AttributeError):
            BreakpointPair(
                Breakpoint('1', 1, strand=STRAND.POS),
                Breakpoint('1', 2, strand=STRAND.POS),
                opposing_strands=True
            )

    def test___init__opstrand_indv_not_specified(self):
        bpp = BreakpointPair(Breakpoint('test', 1), Breakpoint('test', 10), opposing_strands=True)
        self.assertTrue(bpp.opposing_strands)
        bpp = BreakpointPair(Breakpoint('test', 1), Breakpoint('test', 10), opposing_strands=False)
        self.assertFalse(bpp.opposing_strands)

    def test___init__opstrand_not_specified(self):
        with self.assertRaises(NotSpecifiedError):
            BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 2))

    def test___init__stranded(self):
        with self.assertRaises(NotSpecifiedError):
            BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 2), stranded=True, opposing_strands=True)

    def test___get_item__(self):
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertEqual(bpp[0], bp1)
        self.assertEqual(bpp[1], bp2)
        with self.assertRaises(IndexError):
            bpp['?']
        with self.assertRaises(IndexError):
            bpp[2]

    def test_interchromosomal(self):
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertTrue(bpp.interchromosomal)
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(1, 7, 8, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertFalse(bpp.interchromosomal)

    def test___init__invalid_intra_RPRP(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_RNRN(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_RPLN(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
            )

    def test___init__invalid_intra_LPRN(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_RNLP(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
            )

    def test___init__invalid_intra_LNRP(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_inter_RL_opp(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.RIGHT),
                Breakpoint(2, 1, 2, ORIENT.LEFT),
                opposing_strands=True
            )

    def test___init__invalid_inter_LR_opp(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.LEFT),
                Breakpoint(2, 1, 2, ORIENT.RIGHT),
                opposing_strands=True
            )

    def test_classify_inverted_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.LEFT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=True
        )
        BreakpointPair.classify(b)

    def test_classify_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.RIGHT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=False
        )
        BreakpointPair.classify(b)

    def test_classify_inversion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_duplication(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
        )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
        )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
        )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_deletion_or_insertion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]),
                         sorted(BreakpointPair.classify(b)))

    def test_call_breakpoint_pair_single_one_event(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(18, bpp.break2.start)
        self.assertEqual(18, bpp.break2.end)
        self.assertEqual('GGG', bpp.untemplated_sequence)

    def test_call_breakpoint_pair_single_multi_events(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
                           'CTC'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(20, bpp.break1.start)
        self.assertEqual(20, bpp.break1.end)
        self.assertEqual(28, bpp.break2.start)
        self.assertEqual(28, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_sequence)
        self.assertEqual('ACTGAATCGTGGGTAGCTGCTAG', bpp.break1.seq)
        self.assertEqual('GGGCCTGCTC', bpp.break2.seq)

        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.D, 3), (CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(23, bpp.break1.start)
        self.assertEqual(23, bpp.break1.end)
        self.assertEqual(31, bpp.break2.start)
        self.assertEqual(31, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_sequence)

        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10), (CIGAR.D, 3), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(18, bpp.break2.start)
        self.assertEqual(18, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_sequence)

    def test_call_breakpoint_pair_two_indel(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT 1-30     1-?
        # r1  AAATTTCCCgggaattccggatcgatcgat 1-9      1-9
        # r2  aaatttcccgggaattccggaTCGATCGAT 22-30    100-108
        # i   ---------GGGAATTCCGGA--------- 10-21    n/a
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 9), (CIGAR.S, 21)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 21), (CIGAR.M, 9)],
            query_sequence=seq,
            is_reverse=False
        )
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('GGGAATTCCGGA', bpp.untemplated_sequence)
        self.assertEqual(9, bpp.break1.start)
        self.assertEqual(100, bpp.break2.start)
        self.assertEqual('AAATTTCCC', bpp.break1.seq)
        self.assertEqual('TCGATCGAT', bpp.break2.seq)

    def test_call_breakpoint_pair_two_del(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccggaTCGATCGAT
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 21), (CIGAR.M, 9)],
            query_sequence=seq,
            is_reverse=False
        )
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_sequence)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(100, bpp.break2.start)

    def test_call_breakpoint_pair_two_del_overlap(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccGGATCGATCGAT

        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 18), (CIGAR.M, 12)],
            query_sequence=seq,
            is_reverse=False
        )
        self.assertEqual(21, r1.reference_end)
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_sequence)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(103, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATTCCGGA', bpp.break1.seq)
        self.assertEqual('TCGATCGAT', bpp.break2.seq)

    def test_call_breakpoint_pair_two_inversion_overlap(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ------------------GGA---------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.M, 12), (CIGAR.S, 18)],
            query_sequence=reverse_complement(seq),
            is_reverse=True
        )
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_sequence)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(108, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATTCCGGA', bpp.break1.seq)
        self.assertEqual(reverse_complement('TCGATCGAT'), bpp.break2.seq)

    def test_call_breakpoint_pair_two_inversion_gap(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTccggatcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ----------------CC------------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 16), (CIGAR.S, 14)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.M, 12), (CIGAR.S, 18)],
            query_sequence=reverse_complement(seq),
            is_reverse=True
        )
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('CC', bpp.untemplated_sequence)
        self.assertEqual(16, bpp.break1.start)
        self.assertEqual(111, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATT', bpp.break1.seq)
        self.assertEqual(reverse_complement('GGATCGATCGAT'), bpp.break2.seq)

    def test_breakpoint_sequence_homology_LPRP(self):
        b1 = Breakpoint(REF_CHR, 157, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1788, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CAATGC', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

        b1 = Breakpoint(REF_CHR, 589, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 704, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('TTAA', 'ATAGC'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_LPLN(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.NEG, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_LNLP(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.NEG, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_RPRN(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_RNRP(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_close_del(self):
        # ....TT|TT....
        b1 = Breakpoint(REF_CHR, 1001, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1002, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_close_dup(self):
        # ....GATACATTTCTTCTTGAAAA...
        # -------------<=============
        # ===============>-----------
        # -------------CT-CT--------- first break homology
        # ------------T--T----------- second break homology
        b1 = Breakpoint(REF_CHR, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CT', 'TT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_breakpoint_sequence_homology_non_specific_error(self):
        b1 = Breakpoint(REF_CHR, 740, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        with self.assertRaises(AttributeError):
            bpp.breakpoint_sequence_homology(REFERENCE_GENOME)
