import unittest
from mavis.bam.read import convert_events_to_softclipping
from mavis.constants import CIGAR, ORIENT


class Mock:
    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            setattr(self, attr, val)


class TestConvertEventsToSoftclipping(unittest.TestCase):
    
    def test_left_large_deletion(self):
        read = Mock(cigar=[(CIGAR.EQ, 10), (CIGAR.D, 10), (CIGAR.EQ, 40)], query_sequence='A' * 50)
        converted = convert_events_to_softclipping(read, ORIENT.LEFT, 5, 5)
        self.assertEqual([(CIGAR.EQ, 10), (CIGAR.S, 40)], converted.cigar)
    
    def test_left_anchor_after_event(self):
        read = Mock(
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 10), (CIGAR.EQ, 40), (CIGAR.D, 10), (CIGAR.EQ, 6)], query_sequence='A' * 50)
        converted = convert_events_to_softclipping(read, ORIENT.LEFT, 5, 5)
        self.assertEqual([(CIGAR.EQ, 4), (CIGAR.D, 10), (CIGAR.EQ, 40), (CIGAR.S, 6)], converted.cigar)

    def test_left_all_mismatch_error(self):
        read = Mock(cigar=[(CIGAR.X, 10), (CIGAR.D, 10), (CIGAR.X, 40)], query_sequence='A' * 50)
        converted = convert_events_to_softclipping(read, ORIENT.LEFT, 5, 5)
        self.assertEqual(read, converted)

    def test_left_combined_small_events(self):
        read = Mock(cigar=[(CIGAR.EQ, 10), (CIGAR.D, 6), (CIGAR.I, 5), (CIGAR.EQ, 35)], query_sequence='A' * 50)
        converted = convert_events_to_softclipping(read, ORIENT.LEFT, 10, 10)
        self.assertEqual([(CIGAR.EQ, 10), (CIGAR.S, 40)], converted.cigar)

    def test_right_large_deletion(self):
        read = Mock(cigar=[(CIGAR.EQ, 10), (CIGAR.D, 10), (CIGAR.EQ, 40)], query_sequence='A' * 50, reference_start=100)
        converted = convert_events_to_softclipping(read, ORIENT.RIGHT, 5, 5)
        self.assertEqual([(CIGAR.S, 10), (CIGAR.EQ, 40)], converted.cigar)
        self.assertEqual(read.reference_start + 20, converted.reference_start)

    def test_right_anchor_after_event(self):
        read = Mock(
            cigar=[(CIGAR.EQ, 6), (CIGAR.D, 10), (CIGAR.EQ, 40), (CIGAR.D, 10), (CIGAR.EQ, 4)],
            query_sequence='A' * 50, reference_start=100)
        converted = convert_events_to_softclipping(read, ORIENT.RIGHT, 5, 5)
        self.assertEqual([(CIGAR.S, 6), (CIGAR.EQ, 40), (CIGAR.D, 10), (CIGAR.EQ, 4)], converted.cigar)
        self.assertEqual(read.reference_start + 16, converted.reference_start)

    def test_complex_alignment(self):
        cigar = [
            (CIGAR.M, 137), (CIGAR.D, 14823), (CIGAR.M, 19), (CIGAR.D, 1), (CIGAR.M, 5), (CIGAR.I, 18), (CIGAR.D, 18),
            (CIGAR.M, 16), (CIGAR.I, 1), (CIGAR.D, 120), (CIGAR.M, 22), (CIGAR.S, 147)]
        read = Mock(cigar=cigar, query_sequence='A' * 365, reference_start=88217410)
        
        converted = convert_events_to_softclipping(read, ORIENT.LEFT, 50, 50)
        self.assertEqual([(CIGAR.M, 137), (CIGAR.S, 365 - 137)], converted.cigar)
        
        converted = convert_events_to_softclipping(read, ORIENT.RIGHT, 50, 100)
        self.assertEqual(read.cigar, converted.cigar)
