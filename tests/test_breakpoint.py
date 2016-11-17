import unittest
from structural_variant.constants import *
from structural_variant.breakpoint import *
from structural_variant.error import *
from tests import MockRead


class TestBreakpoint(unittest.TestCase):

    def test___eq__(self):
        self.assertNotEqual(Breakpoint('1', 1), None)
        self.assertEqual(Breakpoint('1', 1), Breakpoint('1', 1))


class TestBreakpointPair(unittest.TestCase):

    def test___get_item__(self):
        bp1 = Breakpoint(1, 1, 2, STRAND.NS, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, STRAND.NS, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertEqual(bpp[0], bp1)
        self.assertEqual(bpp[1], bp2)
        with self.assertRaises(IndexError):
            bpp['?']
        with self.assertRaises(IndexError):
            bpp[2]

    def test_interchromosomal(self):
        bp1 = Breakpoint(1, 1, 2, STRAND.NS, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, STRAND.NS, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertTrue(bpp.interchromosomal)
        bp1 = Breakpoint(1, 1, 2, STRAND.NS, ORIENT.LEFT)
        bp2 = Breakpoint(1, 7, 8, STRAND.NS, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        self.assertFalse(bpp.interchromosomal)

    def test_classify_invalid_rearrangement_error_intrachromosomal(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )

    def test_classify_invalid_rearrangement_error_interchromosomal(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, STRAND.NS, ORIENT.RIGHT),
                Breakpoint(2, 1, 2, STRAND.NS, ORIENT.LEFT),
                opposing_strands=True
            )

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, STRAND.NS, ORIENT.LEFT),
                Breakpoint(2, 1, 2, STRAND.NS, ORIENT.RIGHT),
                opposing_strands=True
            )

    def test_classify_inverted_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, STRAND.NS, ORIENT.LEFT),
            Breakpoint(2, 1, 2, STRAND.NS, ORIENT.LEFT),
            opposing_strands=True
        )
        BreakpointPair.classify(b)

    def test_classify_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, STRAND.NS, ORIENT.RIGHT),
            Breakpoint(2, 1, 2, STRAND.NS, ORIENT.LEFT),
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
            reference_start=1,
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
            reference_start=1,
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
            reference_start=1,
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
            reference_start=1,
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
            reference_start=1,
            cigar=[(CIGAR.M, 9), (CIGAR.S, 21)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=100,
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
            reference_start=1,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=100,
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
            reference_start=1,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=100,
            cigar=[(CIGAR.S, 18), (CIGAR.M, 12)],
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
            reference_start=1,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=100,
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
            reference_start=1,
            cigar=[(CIGAR.M, 16), (CIGAR.S, 14)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=100,
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
