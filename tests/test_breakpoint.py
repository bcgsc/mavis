import unittest
from structural_variant.constants import *
from structural_variant.breakpoint import *
from structural_variant.error import *


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
