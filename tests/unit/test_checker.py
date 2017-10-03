import unittest

from mavis.checker import convert_set_to_ranges


class TestConvertSetToRanges(unittest.TestCase):

    def test_single_range(self):
        self.assertEqual('1-3', convert_set_to_ranges({1, 2, 3}))

    def test_multiple_ranges(self):
        self.assertEqual('1-3, 5-8, 10, 12-14', convert_set_to_ranges({1, 2, 3, 5, 6, 7, 8, 10, 12, 13, 14}))
