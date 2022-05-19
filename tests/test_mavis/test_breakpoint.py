from unittest.mock import Mock

import pytest
from mavis.breakpoint import Breakpoint, BreakpointPair, classify_breakpoint_pair
from mavis.constants import ORIENT, STRAND, SVTYPE
from mavis.error import InvalidRearrangement, NotSpecifiedError
from mavis.interval import Interval


class TestBreakpoint:
    def test___eq__(self):
        assert Breakpoint('1', 1) != None  # noqa: E711
        assert Breakpoint('1', 1) == Breakpoint('1', 1)

    def test___hash__(self):
        b = Breakpoint('1', 1, 2)
        c = Breakpoint('1', 1, 2)
        d = Breakpoint('1', 1, 1)

        temp = set()
        temp.add(b)
        temp.add(c)
        temp.add(d)
        assert len(temp) == 2

        temp = dict()
        temp[b] = None
        temp[c] = None
        temp[d] = None
        assert len(temp.keys()) == 2

    def test___len__(self):
        with pytest.raises(AttributeError):
            Breakpoint('11', 87042760, 87041922, orient=ORIENT.LEFT, strand=STRAND.NS)

    def test_inherited_interval_methods(self):
        b = Breakpoint('1', 1, 10)
        assert b[0] == 1
        assert b[1] == 10
        assert len(b) == 10

    def test_breakpoint_constructor(self):
        b = Breakpoint('1', 10, 50)
        assert b[0] == 10
        assert b[1] == 50
        assert Interval.overlaps((1, 10), b)
        assert Interval.overlaps((50, 55), b)
        assert not Interval.overlaps((1, 9), b)


class TestBreakpointPair:
    def test___eq__(self):
        b = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        c = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        assert b is not c
        assert c == b
        d = BreakpointPair(
            Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_seq=''
        )
        assert d != b
        assert None != b  # noqa: E711

    def test___hash__(self):
        b = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        c = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True)
        d = BreakpointPair(
            Breakpoint('1', 1), Breakpoint('1', 3), opposing_strands=True, untemplated_seq=''
        )
        assert b is not c
        temp = dict()
        temp[b] = None
        temp[d] = None
        temp[c] = None
        assert len(temp.keys()) == 2

        temp = set()
        temp.add(b)
        temp.add(c)
        temp.add(d)
        assert len(temp) == 2

    def test___init__swap_break_order(self):
        b1 = Breakpoint('1', 1)
        b2 = Breakpoint('1', 50)
        bpp = BreakpointPair(b1, b2, opposing_strands=True)
        assert b1 == bpp.break1
        assert b2 == bpp.break2
        bpp = BreakpointPair(b2, b1, opposing_strands=True)
        assert b1 == bpp.break1
        assert b2 == bpp.break2

    def test___init__opstrand_conflict(self):
        with pytest.raises(AssertionError):
            BreakpointPair(
                Breakpoint('1', 1, strand=STRAND.POS),
                Breakpoint('1', 2, strand=STRAND.POS),
                opposing_strands=True,
            )

    def test___init__opstrand_indv_not_specified(self):
        bpp = BreakpointPair(Breakpoint('test', 1), Breakpoint('test', 10), opposing_strands=True)
        assert bpp.opposing_strands
        bpp = BreakpointPair(Breakpoint('test', 1), Breakpoint('test', 10), opposing_strands=False)
        assert not bpp.opposing_strands

    def test___init__opstrand_not_specified(self):
        with pytest.raises(NotSpecifiedError):
            BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 2))

    def test___init__stranded(self):
        with pytest.raises(NotSpecifiedError):
            BreakpointPair(
                Breakpoint('1', 1), Breakpoint('1', 2), stranded=True, opposing_strands=True
            )

    def test___get_item__(self):
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        assert bp1 == bpp[0]
        assert bp2 == bpp[1]
        with pytest.raises(IndexError):
            bpp['?']
        with pytest.raises(IndexError):
            bpp[2]

    def test_interchromosomal(self):
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(2, 1, 2, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        assert bpp.interchromosomal
        bp1 = Breakpoint(1, 1, 2, ORIENT.LEFT)
        bp2 = Breakpoint(1, 7, 8, ORIENT.LEFT)
        bpp = BreakpointPair(bp1, bp2, opposing_strands=True)
        assert not bpp.interchromosomal

    def test___init__invalid_intra_rprp(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
                opposing_strands=False,
            )

    def test___init__invalid_intra_rnrn(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                opposing_strands=False,
            )

    def test___init__invalid_intra_rpln(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT),
                opposing_strands=True,
            )

    def test___init__invalid_intra_lprn(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                opposing_strands=True,
            )

    def test___init__invalid_intra_rnlp(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT),
                opposing_strands=True,
            )

    def test___init__invalid_intra_lnrp(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
                Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
                opposing_strands=True,
            )

    def test___init__invalid_inter_rl_opp(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.RIGHT),
                Breakpoint(2, 1, 2, ORIENT.LEFT),
                opposing_strands=True,
            )

    def test___init__invalid_inter_lr_opp(self):
        with pytest.raises(InvalidRearrangement):
            BreakpointPair(
                Breakpoint(1, 1, 2, ORIENT.LEFT),
                Breakpoint(2, 1, 2, ORIENT.RIGHT),
                opposing_strands=True,
            )


class TestClassifyBreakpointPair:
    def test_inverted_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.LEFT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=True,
        )
        classify_breakpoint_pair(b)

    def test_translocation(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, ORIENT.RIGHT),
            Breakpoint(2, 1, 2, ORIENT.LEFT),
            opposing_strands=False,
        )
        classify_breakpoint_pair(b)

    def test_inversion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.INV}

    def test_duplication(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.DUP}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.LEFT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.DUP}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.LEFT),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.DUP}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.NS),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.DUP}

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.NS),
        )
        assert classify_breakpoint_pair(b) == {SVTYPE.DUP}

    def test_deletion_or_insertion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.POS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NEG, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.POS, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NEG, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        b = BreakpointPair(
            Breakpoint(1, 1, 2, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 10, 11, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL, SVTYPE.INS])

    def test_insertion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 2, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.INS])

    def test_no_type(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 2, 2, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
            untemplated_seq='',
        )
        assert classify_breakpoint_pair(b) == set()

    def test_deletion(self):
        b = BreakpointPair(
            Breakpoint(1, 1, 1, strand=STRAND.NS, orient=ORIENT.LEFT),
            Breakpoint(1, 3, 3, strand=STRAND.NS, orient=ORIENT.RIGHT),
            opposing_strands=False,
            untemplated_seq='',
        )
        assert sorted(classify_breakpoint_pair(b)) == sorted([SVTYPE.DEL])

    def test_deletion_with_useq(self):
        bpp = BreakpointPair(
            Breakpoint('1', 6964, orient='L'),
            Breakpoint('1', 7040, orient='R'),
            opposing=False,
            untemplated_seq='CCCT',
        )
        assert sorted(classify_breakpoint_pair(bpp)) == sorted([SVTYPE.DEL, SVTYPE.INS])

        def distance(x, y):
            return Interval(abs(x - y))

        net_size = BreakpointPair.net_size(bpp, distance)
        assert net_size == Interval(-71)
        assert sorted(classify_breakpoint_pair(bpp, distance)) == sorted([SVTYPE.DEL])

    def test_deletion_no_distance_error(self):
        bpp = BreakpointPair(
            Breakpoint('1', 7039, orient='L'), Breakpoint('1', 7040, orient='R'), opposing=False
        )
        assert sorted(classify_breakpoint_pair(bpp)) == sorted([SVTYPE.INS])


class TestNetSize:
    def test_indel(self):
        bpp = BreakpointPair(
            Breakpoint('1', 13, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='TTT',
        )
        assert bpp.net_size() == Interval(1)

    def test_large_indel(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 101, orient=ORIENT.RIGHT),
            untemplated_seq='TTT',
        )
        assert bpp.net_size() == Interval(-87)

    def test_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 11, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='T',
        )
        assert bpp.net_size() == Interval(1)

        bpp = BreakpointPair(
            Breakpoint('1', 11, orient=ORIENT.RIGHT),
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            untemplated_seq='TT',
        )
        assert bpp.net_size() == Interval(2)

    def test_duplication_with_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.RIGHT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq='TTT',
        )
        assert bpp.net_size() == Interval(9)

    def test_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.RIGHT),
            untemplated_seq='',
        )
        assert bpp.net_size() == Interval(-4)

    def test_inversion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq='',
        )
        assert bpp.net_size() == Interval(0)

    def test_inversion_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 10, orient=ORIENT.LEFT),
            Breakpoint('1', 15, orient=ORIENT.LEFT),
            untemplated_seq='TT',
        )
        assert bpp.net_size() == Interval(2)


class TestUntemplatedShift:
    def test_indel(self):
        ref = {
            '1': Mock(
                seq='AGAAAAAAAAAACAGAGTCTATTAAGGCATCTTCTATGGTCAGATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAATGCTTTCTTTTATTAGTTGCCCTGTTTCAGATTCAGCTTTGTATCTATATCACCTGTTAATATGTGTGGACTCACAGAAATGATCATTGAGGGAATGCACCCTGTTTGGGTG'
            )
        }
        bpp = BreakpointPair(
            Breakpoint('1', 40237990 - 40237846, orient=ORIENT.LEFT),
            Breakpoint('1', 40237991 - 40237846, orient=ORIENT.RIGHT),
            untemplated_seq='GTATATATATATATATAT',
        )
        result = bpp.untemplated_shift(ref)
        print(result)
        assert result == (0, 1)
