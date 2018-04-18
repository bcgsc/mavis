import unittest
from unittest.mock import Mock

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import COLUMNS, ORIENT, STRAND, SVTYPE
from mavis.error import InvalidRearrangement, NotSpecifiedError
from mavis.interval import Interval
from mavis.util import read_bpp_from_input_file


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
        d = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_seq='')
        self.assertNotEqual(b, d)
        self.assertNotEqual(b, None)

    def test___hash__(self):
        b = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        c = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        d = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_seq='')
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
        with self.assertRaises(AssertionError):
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

    def test___init__invalid_intra_rprp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
                opposing_strands=False
            )

    def test___init__invalid_intra_rnrn(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                opposing_strands=False
            )

    def test___init__invalid_intra_rpln(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT),
                opposing_strands=True
            )

    def test___init__invalid_intra_lprn(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                opposing_strands=True
            )

    def test___init__invalid_intra_rnlp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT),
                opposing_strands=True
            )

    def test___init__invalid_intra_lnrp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
                opposing_strands=True
            )

    def test___init__invalid_inter_rl_opp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.RIGHT),
                Breakpoint(2, 1, 2, ORIENT.LEFT),
                opposing_strands=True
            )

    def test___init__invalid_inter_lr_opp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.LEFT),
                Breakpoint(2, 1, 2, ORIENT.RIGHT),
                opposing_strands=True
            )

    def test_accessing_data_attributes(self):
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        bpp.data['a'] = 1
        print(bpp.data)
        self.assertEqual(1, bpp.a)
        with self.assertRaises(AttributeError):
            bpp.random_attr

        with self.assertRaises(AttributeError):
            bpp.call_method

        bpp.data[COLUMNS.call_method] = 1
        print(bpp.data)
        self.assertEqual(1, bpp.call_method)

        COLUMNS.call_method = 'bbreak2_call_method'
        bpp.data[COLUMNS.call_method] = 2
        self.assertEqual(2, bpp.call_method)


class TestClassifyBreakpointPair(unittest.TestCase):

    def test_inverted_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.LEFT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=True
        )
        BreakpointPair.classify(b)

    def test_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.RIGHT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=False
        )
        BreakpointPair.classify(b)

    def test_inversion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual({SVTYPE.INV}, BreakpointPair.classify(b))

    def test_duplication(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual({SVTYPE.DUP}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
        )
        self.assertEqual({SVTYPE.DUP}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
        )
        self.assertEqual({SVTYPE.DUP}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
        )
        self.assertEqual({SVTYPE.DUP}, BreakpointPair.classify(b))

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
        )
        self.assertEqual({SVTYPE.DUP}, BreakpointPair.classify(b))

    def test_deletion_or_insertion(self):
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

    def test_insertion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 2, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        self.assertEqual(sorted([SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_no_type(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 2, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False, untemplated_seq=''
        )
        self.assertEqual(set(), BreakpointPair.classify(b))

    def test_deletion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 3, 3, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False, untemplated_seq=''
        )
        self.assertEqual(sorted([SVTYPE.DEL]), sorted(BreakpointPair.classify(b)))

    def test_deletion_with_useq(self):
        bpp = BreakpointPair(Breakpoint('1', 6964, orient='L'), Breakpoint('1', 7040, orient='R'), opposing=False, untemplated_seq='CCCT')
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(bpp)))

        def distance(x, y):
            return Interval(abs(x - y))
        net_size = BreakpointPair.net_size(bpp, distance)
        self.assertEqual(Interval(-71), net_size)
        self.assertEqual(sorted([SVTYPE.DEL]), sorted(BreakpointPair.classify(bpp, distance)))

    def test_deletion_no_distance_error(self):
        bpp = BreakpointPair(Breakpoint('1', 7039, orient='L'), Breakpoint('1', 7040, orient='R'), opposing=False)
        self.assertEqual(sorted([SVTYPE.INS]), sorted(BreakpointPair.classify(bpp)))


class TestNetSize(unittest.TestCase):

    def test_indel(self):
        bpp = BreakpointPair(
            Breakpoint('1', 13, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='TTT'
        )
        self.assertEqual(Interval(1), bpp.net_size())

    def test_large_indel(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 101, orient=ORIENT.RIGHT),
            untemplated_seq='TTT'
        )
        self.assertEqual(Interval(-87), bpp.net_size())

    def test_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 11, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='T'
        )
        self.assertEqual(Interval(1), bpp.net_size())

        bpp = BreakpointPair(
            Breakpoint('1', 11, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='TT'
        )
        self.assertEqual(Interval(2), bpp.net_size())

    def test_duplication_with_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.RIGHT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq='TTT'
        )
        self.assertEqual(Interval(9), bpp.net_size())

    def test_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.RIGHT),
            untemplated_seq=''
        )
        self.assertEqual(Interval(-4), bpp.net_size())

    def test_inversion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq=''
        )
        self.assertEqual(Interval(0), bpp.net_size())

    def test_inversion_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq='TT'
        )
        self.assertEqual(Interval(2), bpp.net_size())


class TestUntemplatedShift(unittest.TestCase):

    def test_indel(self):
        ref = {'1': Mock(seq='AGAAAAAAAAAACAGAGTCTATTAAGGCATCTTCTATGGTCAGATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAATGCTTTCTTTTATTAGTTGCCCTGTTTCAGATTCAGCTTTGTATCTATATCACCTGTTAATATGTGTGGACTCACAGAAATGATCATTGAGGGAATGCACCCTGTTTGGGTG')}
        bpp = BreakpointPair(
            Breakpoint('1', 40237990 - 40237846, orient=ORIENT.LEFT),
            Breakpoint('1', 40237991 - 40237846, orient=ORIENT.RIGHT),
            untemplated_seq='GTATATATATATATATAT'
        )
        result = bpp.untemplated_shift(ref)
        print(result)
        self.assertEqual((0, 1), result)
