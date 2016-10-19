from structural_variant.interval import Interval
from structural_variant.error import DiscontiuousMappingError
import unittest


class TestInterval(unittest.TestCase):

    def test_eq(self):
        self.assertEqual(Interval(1, 2), Interval(1, 2))

    def test_ne(self):
        self.assertNotEqual(Interval(1, 2), Interval(1, 3))

    def test___get_item__(self):
        temp = Interval(1, 2, 3)
        self.assertEqual(1, temp[0])
        self.assertEqual(2, temp[1])
        self.assertEqual(3, temp[2])
        with self.assertRaises(IndexError):
            temp[4]
        with self.assertRaises(IndexError):
            temp[-1]

    def test_overlaps(self):
        left = Interval(-4, 1)
        middle = Interval(0, 10)
        right = Interval(5, 12)
        self.assertFalse(left.overlaps(right))
        self.assertFalse(right.overlaps(left))
        self.assertTrue(left.overlaps(middle))
        self.assertTrue(right.overlaps(middle))
        self.assertTrue(middle.overlaps(left))
        self.assertTrue(middle.overlaps(right))

    def test___len__(self):
        self.assertEqual(5, len(Interval(1, 5)))

    def test_center(self):
        self.assertEqual(3, Interval(1, 5).center)
        self.assertEqual(3.5, Interval(2, 5).center)

    def test_position_in_range(self):
        pos = (12, 12)
        self.assertEqual((2, False), Interval.position_in_range([(1, 2), (3, 6), (7, 15)], pos))
        self.assertEqual((3, True), Interval.position_in_range([(1, 2), (3, 6), (7, 10), (14, 16)], pos))
        self.assertEqual((3, False), Interval.position_in_range([(1, 2), (3, 6), (7, 10)], pos))
        self.assertEqual((0, True), Interval.position_in_range([(15, 16), (17, 19)], pos))

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

        with self.assertRaises(DiscontiuousMappingError) as e:
            Interval.convert_pos(mapping, 15)
        self.assertEqual(210, e.exception.before)
        self.assertEqual(301, e.exception.after)

    def test_union(self):
        l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(1, 10), Interval.union(*l))
        m = l + [Interval(11)]
        self.assertEqual(Interval(1, 11), Interval.union(*m))
        n = Interval.union(*l)

    def test_intersection(self):
        l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        self.assertEqual(Interval(7), Interval.intersection(*l))
        l.append(Interval(11))
        self.assertEqual(None, Interval.intersection(*l))

    def test___sub__(self):
        x, y, z = (Interval(1, 4), Interval(-1, 0), Interval(0, 3))
        self.assertEqual(1, x - y)
        self.assertEqual(-1, y - x)
        self.assertEqual(0, x - z)
        self.assertEqual(0, z - x)
        self.assertEqual(0, y - z)
        self.assertEqual(0, z - y)
