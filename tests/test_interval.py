from structural_variant.interval import Interval
from structural_variant.error import DiscontiuousMappingError
import unittest

class TestInterval(unittest.TestCase):

    def test___init__error(self):
        with self.assertRaises(AttributeError):
            Interval(4, 3)
        with self.assertRaises(AttributeError):
            Interval(3, 4, 0)

    def test___contains__(self):
        self.assertTrue(Interval(1, 2) in Interval(1, 7))
        self.assertFalse(Interval(1, 7) in Interval(1, 2))

    def test_eq(self):
        self.assertEqual(Interval(1, 2), Interval(1, 2))

    def test_ne(self):
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

    def test___len__(self):
        self.assertEqual(5, len(Interval(1, 5)))

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

        with self.assertRaises(DiscontiuousMappingError) as e:
            Interval.convert_pos(mapping, 15)
        self.assertEqual(110, e.exception.before)
        self.assertEqual(201, e.exception.after)

        with self.assertRaises(DiscontiuousMappingError) as e:
            Interval.convert_pos(mapping, 0)
        self.assertEqual(101, e.exception.before)
        self.assertEqual(None, e.exception.after)

        with self.assertRaises(DiscontiuousMappingError) as e:
            Interval.convert_pos(mapping, 80)
        self.assertEqual(310, e.exception.after)
        self.assertEqual(None, e.exception.before)

        mapping = {(41, 50): (101, 110), (21, 30): (201, 210), (1, 10): (301, 310)}

        self.assertEqual(306, Interval.convert_pos(mapping, 5))
        self.assertEqual(110, Interval.convert_pos(mapping, 41))
        self.assertEqual(210, Interval.convert_pos(mapping, 21))
        self.assertEqual(310, Interval.convert_pos(mapping, 1))
        self.assertEqual(309, Interval.convert_pos(mapping, 2))

        with self.assertRaises(DiscontiuousMappingError) as e:
            Interval.convert_pos(mapping, 15)
        self.assertEqual(210, e.exception.before)
        self.assertEqual(301, e.exception.after)

        

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

    def test_convert_pos_ratioed_intervals(self):
        mapping = {(1, 100): (1, 21), (101, 500): (22, 30), (501, 600): (31, 51), (601, 900): (52, 57), (901, 1100): (58, 100)}
        self.assertEqual(11, round(Interval.convert_pos(mapping, 50), 0))
        self.assertEqual(1, round(Interval.convert_pos(mapping, 1), 0))
        self.assertEqual(100, round(Interval.convert_pos(mapping, 1100), 0))

    def test_union(self):
        l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(1, 10), Interval.union(*l))
        m = l + [Interval(11)]
        self.assertEqual(Interval(1, 11), Interval.union(*m))

        with self.assertRaises(AttributeError):
            Interval.union()

    def test_weighted_mean(self):
        m = Interval.weighted_mean((1, 2), (1, 9), (2, 10))
        self.assertEqual(Interval(1, 4), m)
        m = Interval.weighted_mean((1, 1), (10, 10))
        self.assertEqual(Interval(6), m)

    def test_weighted_mean_empty_list_error(self):
        with self.assertRaises(AttributeError):
            Interval.weighted_mean()

    def test_weighted_mean_identical_even_length(self):
        m = Interval.weighted_mean((1, 2), (1, 2), (1, 2))
        self.assertEqual(Interval(1, 2), m)

    def test_weighted_mean_identical_odd_length(self):
        m = Interval.weighted_mean((1, 3), (1, 3), (1, 3))
        self.assertEqual(Interval(1, 3), m)

    def test_intersection(self):
        l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(7), Interval.intersection(*l))
        l.append(Interval(11))
        self.assertEqual(None, Interval.intersection(*l))

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
