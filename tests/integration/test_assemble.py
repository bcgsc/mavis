import unittest

from mavis.assemble import Contig
from mavis.interval import Interval

from . import MockObject


class TestContigRemap(unittest.TestCase):
    def setUp(self):
        self.contig = Contig(' ' * 60, None)
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=10))
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=20))
        self.contig.add_mapped_sequence(MockObject(reference_start=50, reference_end=60))

    def test_depth_even_coverage(self):
        covg = self.contig.remap_depth(Interval(1, 10))
        self.assertEqual(2, covg)

    def test_depth_mixed_coverage(self):
        covg = self.contig.remap_depth(Interval(1, 20))
        self.assertEqual(1.5, covg)

    def test_depth_no_coverage(self):
        covg = self.contig.remap_depth(Interval(21, 49))
        self.assertEqual(0, covg)

    def test_depth_whole_contig_coverage(self):
        self.assertAlmostEqual(40 / 60, self.contig.remap_depth())

    def test_depth_weighted_read(self):
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=10), 5)
        self.assertAlmostEqual(42 / 60, self.contig.remap_depth())

    def test_depth_bad_query_range(self):
        with self.assertRaises(ValueError):
            self.contig.remap_depth(Interval(0, 10))
        with self.assertRaises(ValueError):
            self.contig.remap_depth(Interval(1, len(self.contig.seq) + 1))

    def test_coverage(self):
        self.assertEqual(0.5, self.contig.remap_coverage())
