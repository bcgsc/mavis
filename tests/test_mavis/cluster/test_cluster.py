import unittest

import pytest

from mavis.cluster.cluster import merge_integer_intervals
from mavis.interval import Interval


class TestMergeIntegerIntervals:
    def test_varying_lengths(self):
        m = merge_integer_intervals((1, 2), (1, 9), (2, 10), weight_adjustment=0)
        assert m == Interval(1, 4)

    def test_same_length(self):
        m = merge_integer_intervals((1, 1), (10, 10))
        assert m == Interval(6)

    def test_empty_list_error(self):
        with pytest.raises(AttributeError):
            merge_integer_intervals()

    def test_identical_even_length(self):
        m = merge_integer_intervals((1, 2), (1, 2), (1, 2))
        assert m == Interval(1, 2)

    def test_identical_odd_length(self):
        m = merge_integer_intervals((1, 3), (1, 3), (1, 3))
        assert m == Interval(1, 3)


if __name__ == '__main__':
    unittest.main()
