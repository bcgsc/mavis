from structural_variant.constants import *
from structural_variant.align import *
import pysam
import unittest

class TestAlign(unittest.TestCase):
    """
    test class for functions in the validate namespace 
    that are not associated with a class
    """
    def test_sw_pairwise_alignment(self):
        a = sw_pairwise_alignment('ATGGACTCGGTAAA', 'CGGTAA')[0]
        self.assertEqual(a.reference_start, 7)
        self.assertEqual(a.cigar, [(CIGAR.EQ, 6)])
        self.assertEqual(a.query_sequence, 'CGGTAA')

    def test_build_string_from_reverse_path(self):
        ref = '-mxabdce'
        seq = '-abc'
        t = build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3)])
        self.assertEqual(('e', '-'), t)
        t = build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3), (5,2), (4,2)])
        self.assertEqual(('dce', '-c-'), t)
        t = build_string_from_reverse_path(ref, seq, [(6, 3), (5,2), (4,2), (3, 1), (2, 0)])
        self.assertEqual(('mxabdce', '--ab-c-'), t)
