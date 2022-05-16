from mavis.constants import ORIENT
from mavis.interval import Interval
from mavis.validate.base import Evidence
from mavis.validate.call import _call_interval_by_flanking_coverage

from ..mock import Mock


class CallIntervalByFlankingCoverage:
    def test_invalid_input_attr(self):
        pass

    def test_left(self):
        i = _call_interval_by_flanking_coverage(
            Mock(start=101, end=110),
            ORIENT.LEFT,
            100,
            20,
            distance=Evidence.distance,
            traverse=Evidence.traverse,
        )
        assert i.start == 110
        assert i.end == 180

        i = _call_interval_by_flanking_coverage(
            Mock(start=20, end=80),
            ORIENT.LEFT,
            230,
            40,
            distance=Evidence.distance,
            traverse=Evidence.traverse,
        )
        assert i.start == 80
        assert i.end == 209

    def test_right(self):
        i = _call_interval_by_flanking_coverage(
            Mock(start=101, end=110),
            ORIENT.RIGHT,
            100,
            20,
            distance=Evidence.distance,
            traverse=Evidence.traverse,
        )
        assert i.end == 101
        assert i.start == 31

        i = _call_interval_by_flanking_coverage(
            Mock(start=150, end=200),
            ORIENT.RIGHT,
            230,
            40,
            distance=Evidence.distance,
            traverse=Evidence.traverse,
        )
        assert i.start == 11
        assert i.end == 150


class TestDistanceAndTraverse:
    def test_distance(self):
        assert Evidence.distance(1, 11) == Interval(10)

    def test_traverse_right(self):
        assert Evidence.traverse(1, 10, ORIENT.RIGHT) == Interval(11)

    def test_traverse_left(self):
        assert Evidence.traverse(20, 10, ORIENT.LEFT) == Interval(10)
