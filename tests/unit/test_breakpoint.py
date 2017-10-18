import unittest

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import COLUMNS, ORIENT, STRAND, SVTYPE
from mavis.error import InvalidRearrangement, NotSpecifiedError
from mavis.interval import Interval
from mavis.util import read_bpp_from_input_file

from .mock import Mock


class MockFileHandle(Mock):
    def __init__(self, lines):
        Mock.__init__(self, lines=lines)

    def readlines(self):
        return self.lines


class TestReadBreakpointPairsFromFile(unittest.TestCase):

    def build_filehandle(self, row):
        header = [c for c in row]
        line = [row[c] for c in header]
        lines = ['\t'.join(header), '\t'.join([str(v) for v in line])]
        return MockFileHandle(lines)

    def test_break1_strand_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        with self.assertRaises(AssertionError) as err:
            print(err)
            bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=False)

        bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

    def test_break2_strand_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })

        with self.assertRaises(AssertionError) as err:
            print(err)
            bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=False)

        bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

    def test_both_break_strand_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        with self.assertRaises(AssertionError) as err:
            print(err)
            bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=False)

        bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=True)
        self.assertEqual(2, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

    def test_break1_orient_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.NS,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)

    def test_break2_orient_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.NS,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)
        raise unittest.SkipTest('TODO')

    def test_both_break_orient_ns(self):
        raise unittest.SkipTest('TODO')

    def test_base_case(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.RIGHT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NEG,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: True
        })
        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.RIGHT, bpps[0].break1.orient)
        self.assertEqual(True, bpps[0].opposing_strands)

    def test_unstranded_with_strand_calls(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.RIGHT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NEG,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: True
        })
        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=False, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, explicit_strand=True, expand_ns=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[0].break2.strand)


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
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_rnrn(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_rpln(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
            )

    def test___init__invalid_intra_lprn(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )

    def test___init__invalid_intra_rnlp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
            )

    def test___init__invalid_intra_lnrp(self):
        with self.assertRaises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
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

    def test_duplication(self):
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
