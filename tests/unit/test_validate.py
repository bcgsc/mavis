import unittest
from mavis.validate.call import _call_by_flanking_pairs, _call_interval_by_flanking_coverage
from mavis.constants import STRAND, SVTYPE, PROTOCOL, ORIENT
from .mock import Mock, MockFunction
from mavis.breakpoint import Breakpoint

# evidence req
#   compute_fragment_size
#   max_expected_fragment_size
#   min_flanking_pairs_resolution
#   read_length 
#   decide_sequenced_strand
#   compute_exonic_distance
#   overlapping_transcripts
#   break1 orient chr
#   break2 orient chr
#   interchromosomal

# read req
#   reference_start
#   reference_end
#   next_reference_start

# breakpoint req
#   chr
#   orient
"""
self.ev = Mock(
    compute_fragment_size=MockFunction(1),
    max_expected_fragment_size=500,
    min_flanking_pairs_resolution=1,
    read_length=50,
    decide_sequenced_strand=MockFunction(STRAND.POS),
    overlapping_transcripts=(None, None),
    interchromosomal=False,
    flanking_pairs=[],
    event_type=SVTYPE.INS,
    protocol=PROTOCOL.GENOME,
    break1=None,
    break2=None
)
"""


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

