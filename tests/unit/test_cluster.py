import unittest

from mavis.cluster.cluster import merge_integer_intervals
from mavis.interval import Interval


class TestMergeIntegerIntervals(unittest.TestCase):
    def test_varying_lengths(self):
        m = merge_integer_intervals((1, 2), (1, 9), (2, 10), weight_adjustment=0)
        self.assertEqual(Interval(1, 4), m)

    def test_same_length(self):
        m = merge_integer_intervals((1, 1), (10, 10))
        self.assertEqual(Interval(6), m)

    def test_empty_list_error(self):
        with self.assertRaises(AttributeError):
            merge_integer_intervals()

    def test_identical_even_length(self):
        m = merge_integer_intervals((1, 2), (1, 2), (1, 2))
        self.assertEqual(Interval(1, 2), m)

    def test_identical_odd_length(self):
        m = merge_integer_intervals((1, 3), (1, 3), (1, 3))
        self.assertEqual(Interval(1, 3), m)


if __name__ == '__main__':
    unittest.main()
