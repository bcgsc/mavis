import unittest

from mavis.constants import ORIENT
from mavis.validate.call import _call_interval_by_flanking_coverage
from mavis.validate.evidence import GenomeEvidence
from mavis.validate.base import Evidence
from mavis.interval import Interval

from .mock import Mock


class CallIntervalByFlankingCoverage(unittest.TestCase):

    def test_invalid_input_attr(self):
        pass

    def test_left(self):
        i = _call_interval_by_flanking_coverage(
            Mock(start=101, end=110), ORIENT.LEFT, 100, 20,
            distance=Evidence.distance, traverse=Evidence.traverse)
        self.assertEqual(110, i.start)
        self.assertEqual(180, i.end)

        i = _call_interval_by_flanking_coverage(
            Mock(start=20, end=80), ORIENT.LEFT, 230, 40,
            distance=Evidence.distance, traverse=Evidence.traverse)
        self.assertEqual(80, i.start)
        self.assertEqual(209, i.end)

    def test_right(self):
        i = _call_interval_by_flanking_coverage(
            Mock(start=101, end=110), ORIENT.RIGHT, 100, 20,
            distance=Evidence.distance, traverse=Evidence.traverse)
        self.assertEqual(101, i.end)
        self.assertEqual(31, i.start)

        i = _call_interval_by_flanking_coverage(
            Mock(start=150, end=200), ORIENT.RIGHT, 230, 40,
            distance=Evidence.distance, traverse=Evidence.traverse)
        self.assertEqual(11, i.start)
        self.assertEqual(150, i.end)


class TestDistanceAndTraverse(unittest.TestCase):
    def test_distance(self):
        self.assertEqual(Interval(10), Evidence.distance(1, 11))

    def test_traverse_right(self):
        self.assertEqual(Interval(11), Evidence.traverse(1, 10, ORIENT.RIGHT))

    def test_traverse_left(self):
        self.assertEqual(Interval(10), Evidence.traverse(20, 10, ORIENT.LEFT))
