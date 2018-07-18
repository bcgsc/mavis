import unittest
from mavis.illustrate.util import generate_interval_mapping
from mavis.interval import Interval


class TestGenerateIntervalMapping(unittest.TestCase):

    def test_single_bp_window(self):
        regions = [Interval(4222347, 4222347), Interval(4221673, 4221903), Interval(2792992, 4852494)]
        target = 911.9921875
        ratio = 5
        min_width = 60
        buffer_ = None
        start = 2791992
        end = 4853494
        min_inter = 10
        mapping = generate_interval_mapping(regions, target, ratio, min_width, buffer_, start, end, min_inter)
        self.assertEqual(7, len(mapping.keys()))

    def test_no_input_intervals(self):
        target = 911.9921875
        ratio = 5
        min_width = 60
        buffer_ = None
        start = 2791992
        end = 4853494
        min_inter = 10
        mapping = generate_interval_mapping([], target, ratio, min_width, buffer_, start, end, min_inter)
        self.assertEqual(1, len(mapping.keys()))
