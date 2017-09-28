import unittest

from mavis.constants import ORIENT
from mavis.validate.call import _call_interval_by_flanking_coverage

from .mock import Mock


class CallIntervalByFlankingCoverage(unittest.TestCase):

    def test_invalid_input_attr(self):
        pass

    def test_left(self):
        i = _call_interval_by_flanking_coverage(Mock(start=101, end=110), ORIENT.LEFT, 100, 20)
        self.assertEqual(110, i.start)
        self.assertEqual(180, i.end)

        i = _call_interval_by_flanking_coverage(Mock(start=20, end=80), ORIENT.LEFT, 230, 40)
        self.assertEqual(80, i.start)
        self.assertEqual(209, i.end)

    def test_right(self):
        i = _call_interval_by_flanking_coverage(Mock(start=101, end=110), ORIENT.RIGHT, 100, 20)
        self.assertEqual(101, i.end)
        self.assertEqual(31, i.start)

        i = _call_interval_by_flanking_coverage(Mock(start=150, end=200), ORIENT.RIGHT, 230, 40)
        self.assertEqual(11, i.start)
        self.assertEqual(150, i.end)
