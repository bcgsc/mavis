from structural_variant.interval import Interval
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

    def test__len__(self):
        self.assertEqual(5, len(Interval(1, 5)))

    def test_center(self):
        self.assertEqual(3, Interval(1, 5).center)
        self.assertEqual(3.5, Interval(2, 5).center)
