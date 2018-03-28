import unittest
from mavis.interval import Interval, IntervalMapping


class TestInterval(unittest.TestCase):

    def test___init__error(self):
        with self.assertRaises(AttributeError):
            Interval(4, 3)
        with self.assertRaises(AttributeError):
            Interval(3, 4, 0)

    def test___contains__(self):
        self.assertTrue(Interval(1, 2) in Interval(1, 7))
        self.assertFalse(Interval(1, 7) in Interval(1, 2))
        self.assertTrue(Interval(1.0, 2) in Interval(1.0, 7))
        self.assertFalse(Interval(1, 7) in Interval(1, 2))
        self.assertTrue(1 in Interval(1, 7))
        self.assertFalse(0 in Interval(1, 7))

    def test_eq(self):
        self.assertEqual(Interval(1, 2), Interval(1, 2))
        self.assertEqual(Interval(1, 2), Interval(1, 2))

    def test_ne(self):
        self.assertNotEqual(Interval(1, 2), Interval(1, 3))
        self.assertNotEqual(Interval(1, 2), Interval(1, 3))

    def test___get_item__(self):
        temp = Interval(1, 2, 3)
        self.assertEqual(1, temp[0])
        self.assertEqual(2, temp[1])
        with self.assertRaises(IndexError):
            temp[3]
        with self.assertRaises(IndexError):
            temp[-1]
        with self.assertRaises(IndexError):
            temp['1b']

    def test___gt__(self):
        self.assertTrue(Interval(10) > Interval(1))
        self.assertFalse(Interval(1) > Interval(10))
        self.assertTrue(Interval(10) > Interval(1))
        self.assertFalse(Interval(1) > Interval(1.01))

    def test_overlaps(self):
        left = Interval(-4, 1)
        middle = Interval(0, 10)
        right = Interval(5, 12)
        self.assertFalse(Interval.overlaps(left, right))
        self.assertFalse(Interval.overlaps(right, left))
        self.assertTrue(Interval.overlaps(left, middle))
        self.assertTrue(Interval.overlaps(right, middle))
        self.assertTrue(Interval.overlaps(middle, left))
        self.assertTrue(Interval.overlaps(middle, right))
        self.assertTrue(Interval.overlaps((1, 2), (2, 5)))
        left = Interval(1148432, 1149343)
        right = Interval(1149493, 1150024)
        self.assertFalse(Interval.overlaps(left, right))

        left = Interval(-4, 0.1)
        middle = Interval(0, 10)
        right = Interval(0.11, 12)
        self.assertFalse(Interval.overlaps(left, right))
        self.assertFalse(Interval.overlaps(right, left))
        self.assertTrue(Interval.overlaps(left, middle))
        self.assertTrue(Interval.overlaps(right, middle))
        self.assertTrue(Interval.overlaps(middle, left))
        self.assertTrue(Interval.overlaps(middle, right))

    def test___len__(self):
        self.assertEqual(5, len(Interval(1, 5)))
        with self.assertRaises(TypeError):
            len(Interval(1, 5.0))
        self.assertEqual(4.0, Interval(1, 5.0).length())

    def test___lt__(self):
        self.assertTrue(Interval(1) < Interval(10))
        self.assertFalse(Interval(10) < Interval(1))

    def test___and__(self):
        self.assertEqual(None, Interval(1, 1) & Interval(2))

    def test___sub__(self):
        # x in y
        self.assertEqual([Interval(0, 4), Interval(7, 10)], Interval(0, 10) - Interval(5, 6))
        # x overlaps the start of y
        self.assertEqual([Interval(7, 10)], Interval(0, 10) - Interval(-1, 6))
        # x overlaps the end of y
        self.assertEqual([Interval(0, 4)], Interval(0, 10) - Interval(5, 11))
        # x overlaps all of y
        self.assertEqual([], Interval(0, 10) - Interval(-1, 11))
        # x does not overlap y
        self.assertEqual([Interval(0, 10)], Interval(0, 10) - Interval(11, 15))

    def test___xor__(self):
        # x in y
        self.assertEqual([], Interval(0, 10) ^ Interval(0, 10))
        # x overlaps the start of y
        self.assertEqual([Interval(7, 10), Interval(-1, -1)], Interval(0, 10) ^ Interval(-1, 6))
        # x overlaps the end of y
        self.assertEqual([Interval(0, 4), Interval(11, 11)], Interval(0, 10) ^ Interval(5, 11))
        # x overlaps all of y
        self.assertEqual([Interval(-1, -1), Interval(11, 11)], Interval(0, 10) ^ Interval(-1, 11))
        # x does not overlap y
        self.assertEqual([Interval(0, 10), Interval(11, 15)], Interval(0, 10) ^ Interval(11, 15))

    def test_center(self):
        self.assertEqual(3, Interval(1, 5).center)
        self.assertEqual(3.5, Interval(2, 5).center)

    def test_position_in_range(self):
        pos = (12, 12)
        self.assertEqual((2, False), Interval.position_in_range([(1, 2), (3, 6), (7, 15)], pos))
        self.assertEqual((3, True), Interval.position_in_range([(1, 2), (3, 6), (7, 10), (14, 16)], pos))
        self.assertEqual((3, False), Interval.position_in_range([(1, 2), (3, 6), (7, 10)], pos))
        self.assertEqual((0, True), Interval.position_in_range([(15, 16), (17, 19)], pos))

        with self.assertRaises(AttributeError):
            Interval.position_in_range([], 1)

    def test_convert_pos(self):
        mapping = {(1, 10): (101, 110), (21, 30): (201, 210), (41, 50): (301, 310)}

        self.assertEqual(105, Interval.convert_pos(mapping, 5))
        self.assertEqual(101, Interval.convert_pos(mapping, 1))
        self.assertEqual(310, Interval.convert_pos(mapping, 50))

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 15)

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 0)

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 80)

    def test_convert_pos_forward_to_reverse(self):
        mapping = {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)}

        self.assertEqual(306, Interval.convert_pos(mapping, 5))
        self.assertEqual(110, Interval.convert_pos(mapping, 41))
        self.assertEqual(210, Interval.convert_pos(mapping, 21))
        self.assertEqual(310, Interval.convert_pos(mapping, 1))
        self.assertEqual(309, Interval.convert_pos(mapping, 2))

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 15)

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 51)

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 0)

        with self.assertRaises(IndexError):
            Interval.convert_pos(mapping, 31)

    def test_convert_pos_input_errors(self):
        # test input errors
        with self.assertRaises(AttributeError):  # unequal length
            Interval.convert_pos({(1, 10): (4, 5)}, 3)

        with self.assertRaises(AttributeError):  # overlapping ranges
            Interval.convert_pos({(1, 10): (11, 20), (5, 14): (21, 30)}, 6)

        with self.assertRaises(AttributeError):  # range not increasing or decreasing
            mapping = {(1, 2): (1, 2), (3, 4): (4, 5), (5, 6): (3, 3)}
            Interval.convert_pos(mapping, 10)

        with self.assertRaises(AttributeError):  # range not increasing or decreasing
            mapping = {(1, 2): (4, 5), (3, 4): (1, 2), (5, 6): (3, 3)}
            Interval.convert_pos(mapping, 10)

    def test_convert_pos_one_to_one(self):
        mapping = {}
        for x in range(0, 10):
            s = x * 10 + 1
            mapping[Interval(s, s + 9)] = Interval(s, s + 9)
        for pos in range(1, 101):
            self.assertEqual(pos, Interval.convert_pos(mapping, pos))

    def test_convert_pos_ratioed_intervals(self):
        mapping = {
            (1.0, 100): (1, 20.0),
            (101.0, 500): (21.0, 30),
            (501.0, 600): (31.0, 51),
            (601.0, 900): (52, 57.0),
            (901.0, 1100): (58.0, 100)
        }
        self.assertEqual(Interval(1), Interval.convert_ratioed_pos(mapping, 1))
        self.assertEqual(Interval(20), Interval.convert_ratioed_pos(mapping, 100))
        self.assertEqual(Interval(100, 100), Interval.convert_ratioed_pos(mapping, 1100))

        mapping = {(1, 100): (1, 1), (101, 500): (21, 30)}
        self.assertEqual(Interval(1, 1), Interval.convert_ratioed_pos(mapping, 1))
        self.assertEqual(Interval(1, 1), Interval.convert_ratioed_pos(mapping, 100))

        mapping = {(1, 100.0): (20.0, 30), (100.1, 500): (1.0, 1.0)}
        self.assertEqual(Interval(1, 1), Interval.convert_ratioed_pos(mapping, 101))
        self.assertEqual(Interval(1, 1), Interval.convert_ratioed_pos(mapping, 500))
        self.assertEqual(Interval(25, 25), Interval.convert_ratioed_pos(mapping, 50))

    def test_union(self):
        interval_list = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(1, 10), Interval.union(*interval_list))
        m = interval_list + [Interval(11)]
        self.assertEqual(Interval(1, 11), Interval.union(*m))

        with self.assertRaises(AttributeError):
            Interval.union()

    def test_intersection(self):
        interval_list = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(7), Interval.intersection(*interval_list))
        interval_list.append(Interval(11))
        self.assertEqual(None, Interval.intersection(*interval_list))

        with self.assertRaises(AttributeError):
            Interval.intersection()

    def test_dist(self):
        x = Interval(1, 4)
        y = Interval(-1, 0)
        z = Interval(0, 3)
        self.assertEqual(1, Interval.dist(x, y))
        self.assertEqual(-1, Interval.dist(y, x))
        self.assertEqual(0, Interval.dist(x, z))
        self.assertEqual(0, Interval.dist(z, x))
        self.assertEqual(0, Interval.dist(y, z))
        self.assertEqual(0, Interval.dist(z, y))
        self.assertEqual(-6, Interval.dist((1, 4), (10, 12)))

    def test_min_nonoverlapping(self):
        r = Interval.min_nonoverlapping(Interval(1, 2), Interval(4, 7), Interval(8, 9))
        self.assertEqual(3, len(r))
        r = Interval.min_nonoverlapping(Interval(1, 5), Interval(4, 7), Interval(8, 9))
        self.assertEqual(2, len(r))
        r = Interval.min_nonoverlapping(Interval(1, 5), Interval(4, 7), Interval(7, 9))
        self.assertEqual([Interval(1, 9)], r)
        r = Interval.min_nonoverlapping((1, 2), (2, 4))
        self.assertEqual([Interval(1, 4)], r)
        self.assertEqual([], Interval.min_nonoverlapping())

    def test_split_overlapping_no_weight(self):
        input_intervals = [
            Interval(1, 10),
            Interval(2, 11),
            Interval(4, 5),
            Interval(4, 8)
        ]
        exp = [
            Interval(1, 1),
            Interval(2, 3),
            Interval(4, 4),
            Interval(5, 7),
            Interval(8, 9),
            Interval(10, 11)
        ]
        print('expected', exp)
        result = Interval.split_overlap(*input_intervals)
        result = sorted(result)
        print('found', result)
        self.assertEqual(exp, result)

    def test_split_overlapping_weighted(self):
        input_intervals = [
            Interval(1, 10),
            Interval(2, 11),
            Interval(4, 5),
            Interval(4, 8)
        ]
        weights = {
            Interval(1, 10): 10,
            Interval(2, 11): 20,
            Interval(4, 5): 10,
            Interval(4, 8): 10
        }
        exp = {
            Interval(1, 1): 1,
            Interval(2, 3): 4,
            Interval(4, 4): 5,
            Interval(5, 7): 6,
            Interval(8, 9): 4,
            Interval(10, 11): 4
        }
        result = Interval.split_overlap(*input_intervals, weight_mapping=weights)
        self.assertEqual(sorted(exp), sorted(result))
        for itvl in exp:
            self.assertEqual(exp[itvl], result[itvl])


class TestIntervalMapping(unittest.TestCase):

    def test_convert_pos_ratioed(self):
        mapping = IntervalMapping({
            (1.0, 100): (1, 20.0),
            (101.0, 500): (21.0, 30),
            (501.0, 600): (31.0, 51),
            (601.0, 900): (52, 57.0),
            (901.0, 1100): (58.0, 100)
        })
        self.assertEqual(1, mapping.convert_pos(1))
        self.assertEqual(1, mapping.convert_ratioed_pos(1).start)
        self.assertAlmostEqual(1.191919191919, mapping.convert_ratioed_pos(1).end)
        self.assertEqual(20, mapping.convert_pos(100))
        self.assertEqual(20, mapping.convert_ratioed_pos(100).start)
        self.assertEqual(100, mapping.convert_pos(1100))
        self.assertEqual(100, mapping.convert_ratioed_pos(1100).start)

        mapping = IntervalMapping({(1, 100): (1, 1.0), (101, 500): (21.0, 30)})
        self.assertEqual(1, mapping.convert_pos(1))
        self.assertEqual(1, mapping.convert_pos(100))

        mapping = IntervalMapping({(1, 100.0): (20.0, 30), (100.1, 500): (1.0, 1.0)})
        self.assertEqual(1, mapping.convert_pos(101))
        self.assertEqual(1, mapping.convert_pos(500))
        self.assertEqual(25, mapping.convert_pos(50))

    def test_convert_pos(self):
        mapping = IntervalMapping({(1, 10): (101, 110), (21, 30): (201, 210), (41, 50): (301, 310)})

        self.assertEqual(105, mapping.convert_pos(5))
        self.assertEqual(101, mapping.convert_pos(1))
        self.assertEqual(310, mapping.convert_pos(50))

        with self.assertRaises(IndexError):
            mapping.convert_pos(15)

        with self.assertRaises(IndexError):
            mapping.convert_pos(0)

        with self.assertRaises(IndexError):
            mapping.convert_pos(80)

    def test_convert_pos_forward_to_reverse(self):
        mapping = IntervalMapping(
            {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)},
            opposing=[(41, 50), (21, 30), (1, 10)]
        )

        self.assertEqual(306, mapping.convert_pos(5))
        self.assertEqual(110, mapping.convert_pos(41))
        self.assertEqual(210, mapping.convert_pos(21))
        self.assertEqual(310, mapping.convert_pos(1))
        self.assertEqual(309, mapping.convert_pos(2))

        with self.assertRaises(IndexError):
            mapping.convert_pos(15)

        with self.assertRaises(IndexError):
            mapping.convert_pos(51)

        with self.assertRaises(IndexError):
            mapping.convert_pos(0)

        with self.assertRaises(IndexError):
            mapping.convert_pos(31)

    def test_convert_pos_one_to_one(self):
        mapping = {}
        for x in range(0, 10):
            s = x * 10 + 1
            mapping[Interval(s, s + 9)] = Interval(s, s + 9)
        mapping = IntervalMapping(mapping)
        for pos in range(1, 101):
            self.assertEqual(pos, mapping.convert_pos(pos))
