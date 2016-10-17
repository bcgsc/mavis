import unittest
from structural_variant.constants import *
from structural_variant.breakpoint import *
from structural_variant.error import *


class TestBreakpointPair(unittest.TestCase):

    def test_classify_invalid_rearrangement_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
            )
            BreakpointPair.classify(b)

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )
            BreakpointPair.classify(b)

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
            )
            BreakpointPair.classify(b)

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
            )
            BreakpointPair.classify(b)

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
            )
            BreakpointPair.classify(b)

        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
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
