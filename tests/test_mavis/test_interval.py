import pytest

from mavis.interval import Interval, IntervalMapping


class TestInterval:
    def test___init__error(self):
        with pytest.raises(AttributeError):
            Interval(4, 3)
        with pytest.raises(AttributeError):
            Interval(3, 4, 0)

    def test___contains__(self):
        assert Interval(1, 2) in Interval(1, 7)
        assert not Interval(1, 7) in Interval(1, 2)
        assert Interval(1.0, 2) in Interval(1.0, 7)
        assert not Interval(1, 7) in Interval(1, 2)
        assert 1 in Interval(1, 7)
        assert 0 not in Interval(1, 7)

    def test_eq(self):
        assert Interval(1, 2) == Interval(1, 2)
        assert Interval(1, 2) == Interval(1, 2)

    def test_ne(self):
        assert Interval(1, 2) != Interval(1, 3)
        assert Interval(1, 2) != Interval(1, 3)

    def test___get_item__(self):
        temp = Interval(1, 2, 3)
        assert temp[0] == 1
        assert temp[1] == 2
        with pytest.raises(IndexError):
            temp[3]
        with pytest.raises(IndexError):
            temp[-1]
        with pytest.raises(IndexError):
            temp['1b']

    def test___gt__(self):
        assert Interval(10) > Interval(1)
        assert not Interval(1) > Interval(10)
        assert Interval(10) > Interval(1)
        assert not Interval(1) > Interval(1.01)

    def test_overlaps(self):
        left = Interval(-4, 1)
        middle = Interval(0, 10)
        right = Interval(5, 12)
        assert not Interval.overlaps(left, right)
        assert not Interval.overlaps(right, left)
        assert Interval.overlaps(left, middle)
        assert Interval.overlaps(right, middle)
        assert Interval.overlaps(middle, left)
        assert Interval.overlaps(middle, right)
        assert Interval.overlaps((1, 2), (2, 5))
        left = Interval(1148432, 1149343)
        right = Interval(1149493, 1150024)
        assert not Interval.overlaps(left, right)

        left = Interval(-4, 0.1)
        middle = Interval(0, 10)
        right = Interval(0.11, 12)
        assert not Interval.overlaps(left, right)
        assert not Interval.overlaps(right, left)
        assert Interval.overlaps(left, middle)
        assert Interval.overlaps(right, middle)
        assert Interval.overlaps(middle, left)
        assert Interval.overlaps(middle, right)

    def test___len__(self):
        assert len(Interval(1, 5)) == 5
        with pytest.raises(TypeError):
            len(Interval(1, 5.0))
        assert Interval(1, 5.0).length() == 4.0

    def test___lt__(self):
        assert Interval(1) < Interval(10)
        assert not Interval(10) < Interval(1)

    def test___and__(self):
        assert Interval(1, 1) & Interval(2) is None

    def test___sub__(self):
        # x in y
        assert Interval(0, 10) - Interval(5, 6) == [Interval(0, 4), Interval(7, 10)]
        # x overlaps the start of y
        assert Interval(0, 10) - Interval(-1, 6) == [Interval(7, 10)]
        # x overlaps the end of y
        assert Interval(0, 10) - Interval(5, 11) == [Interval(0, 4)]
        # x overlaps all of y
        assert Interval(0, 10) - Interval(-1, 11) == []
        # x does not overlap y
        assert Interval(0, 10) - Interval(11, 15) == [Interval(0, 10)]

    def test___xor__(self):
        # x in y
        assert Interval(0, 10) ^ Interval(0, 10) == []
        # x overlaps the start of y
        assert Interval(0, 10) ^ Interval(-1, 6) == [Interval(7, 10), Interval(-1, -1)]
        # x overlaps the end of y
        assert Interval(0, 10) ^ Interval(5, 11) == [Interval(0, 4), Interval(11, 11)]
        # x overlaps all of y
        assert Interval(0, 10) ^ Interval(-1, 11) == [Interval(-1, -1), Interval(11, 11)]
        # x does not overlap y
        assert Interval(0, 10) ^ Interval(11, 15) == [Interval(0, 10), Interval(11, 15)]

    def test_center(self):
        assert Interval(1, 5).center == 3
        assert Interval(2, 5).center == 3.5

    def test_position_in_range(self):
        pos = (12, 12)
        assert Interval.position_in_range([(1, 2), (3, 6), (7, 15)], pos) == (2, False)
        assert Interval.position_in_range([(1, 2), (3, 6), (7, 10), (14, 16)], pos) == (3, True)
        assert Interval.position_in_range([(1, 2), (3, 6), (7, 10)], pos) == (3, False)
        assert Interval.position_in_range([(15, 16), (17, 19)], pos) == (0, True)

        with pytest.raises(AttributeError):
            Interval.position_in_range([], 1)

    def test_convert_pos(self):
        mapping = {(1, 10): (101, 110), (21, 30): (201, 210), (41, 50): (301, 310)}

        assert Interval.convert_pos(mapping, 5) == 105
        assert Interval.convert_pos(mapping, 1) == 101
        assert Interval.convert_pos(mapping, 50) == 310

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 15)

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 0)

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 80)

    def test_convert_pos_forward_to_reverse(self):
        mapping = {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)}

        assert Interval.convert_pos(mapping, 5) == 306
        assert Interval.convert_pos(mapping, 41) == 110
        assert Interval.convert_pos(mapping, 21) == 210
        assert Interval.convert_pos(mapping, 1) == 310
        assert Interval.convert_pos(mapping, 2) == 309

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 15)

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 51)

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 0)

        with pytest.raises(IndexError):
            Interval.convert_pos(mapping, 31)

    def test_convert_pos_input_errors(self):
        # test input errors
        with pytest.raises(AttributeError):  # unequal length
            Interval.convert_pos({(1, 10): (4, 5)}, 3)

        with pytest.raises(AttributeError):  # overlapping ranges
            Interval.convert_pos({(1, 10): (11, 20), (5, 14): (21, 30)}, 6)

        with pytest.raises(AttributeError):  # range not increasing or decreasing
            mapping = {(1, 2): (1, 2), (3, 4): (4, 5), (5, 6): (3, 3)}
            Interval.convert_pos(mapping, 10)

        with pytest.raises(AttributeError):  # range not increasing or decreasing
            mapping = {(1, 2): (4, 5), (3, 4): (1, 2), (5, 6): (3, 3)}
            Interval.convert_pos(mapping, 10)

    def test_convert_pos_one_to_one(self):
        mapping = {}
        for x in range(0, 10):
            s = x * 10 + 1
            mapping[Interval(s, s + 9)] = Interval(s, s + 9)
        for pos in range(1, 101):
            assert Interval.convert_pos(mapping, pos) == pos

    def test_convert_pos_ratioed_intervals(self):
        mapping = {
            (1.0, 100): (1, 20.0),
            (101.0, 500): (21.0, 30),
            (501.0, 600): (31.0, 51),
            (601.0, 900): (52, 57.0),
            (901.0, 1100): (58.0, 100),
        }
        assert Interval.convert_ratioed_pos(mapping, 1) == Interval(1)
        assert Interval.convert_ratioed_pos(mapping, 100) == Interval(20)
        assert Interval.convert_ratioed_pos(mapping, 1100) == Interval(100, 100)

        mapping = {(1, 100): (1, 1), (101, 500): (21, 30)}
        assert Interval.convert_ratioed_pos(mapping, 1) == Interval(1, 1)
        assert Interval.convert_ratioed_pos(mapping, 100) == Interval(1, 1)

        mapping = {(1, 100.0): (20.0, 30), (100.1, 500): (1.0, 1.0)}
        assert Interval.convert_ratioed_pos(mapping, 101) == Interval(1, 1)
        assert Interval.convert_ratioed_pos(mapping, 500) == Interval(1, 1)
        assert Interval.convert_ratioed_pos(mapping, 50) == Interval(25, 25)

    def test_union(self):
        interval_list = [Interval(1, 10), Interval(5, 7), Interval(7)]
        assert Interval.union(*interval_list) == Interval(1, 10)
        m = interval_list + [Interval(11)]
        assert Interval.union(*m) == Interval(1, 11)

        with pytest.raises(AttributeError):
            Interval.union()

    def test_intersection(self):
        interval_list = [Interval(1, 10), Interval(5, 7), Interval(7)]
        assert Interval.intersection(*interval_list) == Interval(7)
        interval_list.append(Interval(11))
        assert Interval.intersection(*interval_list) is None

        with pytest.raises(AttributeError):
            Interval.intersection()

    def test_dist(self):
        x = Interval(1, 4)
        y = Interval(-1, 0)
        z = Interval(0, 3)
        assert Interval.dist(x, y) == 1
        assert Interval.dist(y, x) == -1
        assert Interval.dist(x, z) == 0
        assert Interval.dist(z, x) == 0
        assert Interval.dist(y, z) == 0
        assert Interval.dist(z, y) == 0
        assert Interval.dist((1, 4), (10, 12)) == -6

    def test_min_nonoverlapping(self):
        r = Interval.min_nonoverlapping(Interval(1, 2), Interval(4, 7), Interval(8, 9))
        assert len(r) == 3
        r = Interval.min_nonoverlapping(Interval(1, 5), Interval(4, 7), Interval(8, 9))
        assert len(r) == 2
        r = Interval.min_nonoverlapping(Interval(1, 5), Interval(4, 7), Interval(7, 9))
        assert r == [Interval(1, 9)]
        r = Interval.min_nonoverlapping((1, 2), (2, 4))
        assert r == [Interval(1, 4)]
        assert Interval.min_nonoverlapping() == []

    def test_split_overlapping_no_weight(self):
        input_intervals = [Interval(1, 10), Interval(2, 11), Interval(4, 5), Interval(4, 8)]
        exp = [
            Interval(1, 1),
            Interval(2, 3),
            Interval(4, 4),
            Interval(5, 7),
            Interval(8, 9),
            Interval(10, 11),
        ]
        print('expected', exp)
        result = Interval.split_overlap(*input_intervals)
        result = sorted(result)
        print('found', result)
        assert result == exp

    def test_split_overlapping_weighted(self):
        input_intervals = [Interval(1, 10), Interval(2, 11), Interval(4, 5), Interval(4, 8)]
        weights = {Interval(1, 10): 10, Interval(2, 11): 20, Interval(4, 5): 10, Interval(4, 8): 10}
        exp = {
            Interval(1, 1): 1,
            Interval(2, 3): 4,
            Interval(4, 4): 5,
            Interval(5, 7): 6,
            Interval(8, 9): 4,
            Interval(10, 11): 4,
        }
        result = Interval.split_overlap(*input_intervals, weight_mapping=weights)
        assert sorted(result) == sorted(exp)
        for itvl in exp:
            assert result[itvl] == exp[itvl]


class TestIntervalMapping:
    def test_convert_pos_ratioed(self):
        mapping = IntervalMapping(
            {
                (1.0, 100): (1, 20.0),
                (101.0, 500): (21.0, 30),
                (501.0, 600): (31.0, 51),
                (601.0, 900): (52, 57.0),
                (901.0, 1100): (58.0, 100),
            }
        )
        assert mapping.convert_pos(1) == 1
        assert mapping.convert_ratioed_pos(1).start == 1
        assert pytest.approx(mapping.convert_ratioed_pos(1).end) == 1.191919191919
        assert mapping.convert_pos(100) == 20
        assert mapping.convert_ratioed_pos(100).start == 20
        assert mapping.convert_pos(1100) == 100
        assert mapping.convert_ratioed_pos(1100).start == 100

        mapping = IntervalMapping({(1, 100): (1, 1.0), (101, 500): (21.0, 30)})
        assert mapping.convert_pos(1) == 1
        assert mapping.convert_pos(100) == 1

        mapping = IntervalMapping({(1, 100.0): (20.0, 30), (100.1, 500): (1.0, 1.0)})
        assert mapping.convert_pos(101) == 1
        assert mapping.convert_pos(500) == 1
        assert mapping.convert_pos(50) == 25

    def test_convert_pos(self):
        mapping = IntervalMapping({(1, 10): (101, 110), (21, 30): (201, 210), (41, 50): (301, 310)})

        assert mapping.convert_pos(5) == 105
        assert mapping.convert_pos(1) == 101
        assert mapping.convert_pos(50) == 310

        with pytest.raises(IndexError):
            mapping.convert_pos(15)

        with pytest.raises(IndexError):
            mapping.convert_pos(0)

        with pytest.raises(IndexError):
            mapping.convert_pos(80)

    def test_convert_pos_forward_to_reverse(self):
        mapping = IntervalMapping(
            {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)},
            opposing=[(41, 50), (21, 30), (1, 10)],
        )

        assert mapping.convert_pos(5) == 306
        assert mapping.convert_pos(41) == 110
        assert mapping.convert_pos(21) == 210
        assert mapping.convert_pos(1) == 310
        assert mapping.convert_pos(2) == 309

        with pytest.raises(IndexError):
            mapping.convert_pos(15)

        with pytest.raises(IndexError):
            mapping.convert_pos(51)

        with pytest.raises(IndexError):
            mapping.convert_pos(0)

        with pytest.raises(IndexError):
            mapping.convert_pos(31)

    def test_convert_pos_one_to_one(self):
        mapping = {}
        for x in range(0, 10):
            s = x * 10 + 1
            mapping[Interval(s, s + 9)] = Interval(s, s + 9)
        mapping = IntervalMapping(mapping)
        for pos in range(1, 101):
            assert mapping.convert_pos(pos) == pos
