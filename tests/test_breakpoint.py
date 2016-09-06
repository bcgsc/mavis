import unittest
from structural_variant.constants import *
from structural_variant.breakpoint import *
from structural_variant.error import *

class TestBreakpoint(unittest.TestCase):
    pass

class TestBreakpointPair(unittest.TestCase):
    def test_classify_intra_rr_pp_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_rr_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_rr_pu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))
    

    def test_classify_intra_rr_np(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_rr_nn_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_rr_nu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))
    

    def test_classify_intra_rr_up_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))
    

    def test_classify_intra_rr_un_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_rr_uu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_rl_pp_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_rl_pn_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
            BreakpointPair.classify(b)
    

    def test_classify_intra_rl_pu_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_rl_np_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
            BreakpointPair.classify(b)
    

    def test_classify_intra_rl_nn_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_rl_nu_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))
    

    def test_classify_intra_rl_up_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))
    

    def test_classify_intra_rl_un_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))
    

    def test_classify_intra_rl_uu_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))
    

    def test_classify_intra_ru_pp_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))
    

    def test_classify_intra_ru_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ru_pu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ru_np_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ru_nn_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_ru_nu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ru_up_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ru_un_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ru_uu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_pp_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_pn_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_lr_pu_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_np_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_lr_nn_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_nu_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_up_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_un_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lr_uu_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ll_pp_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_ll_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_pu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_np_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_nn_error(self):
        with self.assertRaises(InvalidRearrangement):
            b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
            BreakpointPair.classify(b)

    def test_classify_intra_ll_nu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_up_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_un_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ll_uu_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_lu_pp_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_lu_pu_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_np(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_lu_nn_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_nu_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_up_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_un_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_lu_uu(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.NS)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_pp_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ur_pu_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_np_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ur_nn_del_ins(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_nu_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_up_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_un_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ur_uu_del_ins_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT)
                    )
        self.assertEqual(sorted([SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ul_pp_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_ul_pn_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ul_pu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ul_np_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.INV], BreakpointPair.classify(b))

    def test_classify_intra_ul_nn_dup(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual([SVTYPE.DUP], BreakpointPair.classify(b))

    def test_classify_intra_ul_nu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ul_up_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ul_un_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_ul_uu_dup_inv(self):
        b = BreakpointPair(
                    Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.NS),
                    Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.LEFT)
                    )
        self.assertEqual(sorted([SVTYPE.DUP, SVTYPE.INV]), sorted(BreakpointPair.classify(b)))

    def test_classify_intra_uu_pp_del_dup_ins(self):
        pass

    def test_classify_intra_uu_pn_inv(self):
        pass

    def test_classify_intra_uu_pu_del_dup_ins_inv(self):
        pass

    def test_classify_intra_uu_np_inv(self):
        pass

    def test_classify_intra_uu_nn_del_dup_ins(self):
        pass

    def test_classify_intra_uu_nu_del_dup_ins_inv(self):
        pass

    def test_classify_intra_uu_up_del_dup_ins_inv(self):
        pass

    def test_classify_intra_uu_un_del_dup_ins_inv(self):
        pass

    def test_classify_intra_uu_uu_del_dup_ins_inv(self):
        pass

    def test_classify_inter_rr_pp_error(self):
        pass

    def test_classify_inter_rr_pn(self):
        pass

    def test_classify_inter_rr_pu(self):
        pass

    def test_classify_inter_rr_np(self):
        pass

    def test_classify_inter_rr_nn_error(self):
        pass

    def test_classify_inter_rr_nu(self):
        pass

    def test_classify_inter_rr_up(self):
        pass

    def test_classify_inter_rr_un(self):
        pass

    def test_classify_inter_rr_uu(self):
        pass

    def test_classify_inter_rl_pp(self):
        pass

    def test_classify_inter_rl_pn(self):
        pass

    def test_classify_inter_rl_pu(self):
        pass

    def test_classify_inter_rl_np(self):
        pass

    def test_classify_inter_rl_nn(self):
        pass

    def test_classify_inter_rl_nu(self):
        pass

    def test_classify_inter_rl_up(self):
        pass

    def test_classify_inter_rl_un(self):
        pass

    def test_classify_inter_rl_uu(self):
        pass

    def test_classify_inter_ru_pp(self):
        pass

    def test_classify_inter_ru_pn(self):
        pass

    def test_classify_inter_ru_pu(self):
        pass

    def test_classify_inter_ru_np(self):
        pass

    def test_classify_inter_ru_nn(self):
        pass

    def test_classify_inter_ru_nu(self):
        pass

    def test_classify_inter_ru_up(self):
        pass

    def test_classify_inter_ru_un(self):
        pass

    def test_classify_inter_ru_uu(self):
        pass

    def test_classify_inter_lr_pp(self):
        pass

    def test_classify_inter_lr_pn(self):
        pass

    def test_classify_inter_lr_pu(self):
        pass

    def test_classify_inter_lr_np(self):
        pass

    def test_classify_inter_lr_nn(self):
        pass

    def test_classify_inter_lr_nu(self):
        pass

    def test_classify_inter_lr_up(self):
        pass

    def test_classify_inter_lr_un(self):
        pass

    def test_classify_inter_lr_uu(self):
        pass

    def test_classify_inter_ll_pp_error(self):
        pass

    def test_classify_inter_ll_pn(self):
        pass

    def test_classify_inter_ll_pu(self):
        pass

    def test_classify_inter_ll_np(self):
        pass

    def test_classify_inter_ll_nn_error(self):
        pass

    def test_classify_inter_ll_nu(self):
        pass

    def test_classify_inter_ll_up(self):
        pass

    def test_classify_inter_ll_un(self):
        pass

    def test_classify_inter_ll_uu(self):
        pass

    def test_classify_inter_lu_pp(self):
        pass

    def test_classify_inter_lu_pn(self):
        pass

    def test_classify_inter_lu_pu(self):
        pass

    def test_classify_inter_lu_np(self):
        pass

    def test_classify_inter_lu_nn(self):
        pass

    def test_classify_inter_lu_nu(self):
        pass

    def test_classify_inter_lu_up(self):
        pass

    def test_classify_inter_lu_un(self):
        pass

    def test_classify_inter_lu_uu(self):
        pass

    def test_classify_inter_ur_pp(self):
        pass

    def test_classify_inter_ur_pn(self):
        pass

    def test_classify_inter_ur_pu(self):
        pass

    def test_classify_inter_ur_np(self):
        pass

    def test_classify_inter_ur_nn(self):
        pass

    def test_classify_inter_ur_nu(self):
        pass

    def test_classify_inter_ur_up(self):
        pass

    def test_classify_inter_ur_un(self):
        pass

    def test_classify_inter_ur_uu(self):
        pass

    def test_classify_inter_ul_pp(self):
        pass

    def test_classify_inter_ul_pn(self):
        pass

    def test_classify_inter_ul_pu(self):
        pass

    def test_classify_inter_ul_np(self):
        pass

    def test_classify_inter_ul_nn(self):
        pass

    def test_classify_inter_ul_nu(self):
        pass

    def test_classify_inter_ul_up(self):
        pass

    def test_classify_inter_ul_un(self):
        pass

    def test_classify_inter_ul_uu(self):
        pass

    def test_classify_inter_uu_pp(self):
        pass

    def test_classify_inter_uu_pn(self):
        pass

    def test_classify_inter_uu_pu(self):
        pass

    def test_classify_inter_uu_np(self):
        pass

    def test_classify_inter_uu_nn(self):
        pass

    def test_classify_inter_uu_nu(self):
        pass

    def test_classify_inter_uu_up(self):
        pass

    def test_classify_inter_uu_un(self):
        pass

    def test_classify_inter_uu_uu(self):
        pass
