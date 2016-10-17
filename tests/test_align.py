from structural_variant.constants import *
from structural_variant.align import *
import unittest


class TestAlign(unittest.TestCase):
    """
    test class for functions in the validate namespace
    that are not associated with a class
    """

    def test_alphabet_matching(self):
        self.assertTrue(DNA_ALPHABET.match('N', 'A'))
        self.assertTrue(DNA_ALPHABET.match('A', 'N'))


class TestCigarTools(unittest.TestCase):

    def test_longest_fuzzy_match(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(15, CigarTools.longest_fuzzy_match(c, 1))
        self.assertEqual(10, CigarTools.longest_fuzzy_match(c, 0))
        self.assertEqual(16, CigarTools.longest_fuzzy_match(c, 2))
        self.assertEqual(16, CigarTools.longest_fuzzy_match(c, 4))

    def test_score(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(22, CigarTools.score(c))

    def test_match_percent(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(0.8, CigarTools.match_percent(c))

    def test_compute(self):
        self.assertEqual(
            ([(4, 7), (7, 8)], 7),
            CigarTools.compute('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
        )
        self.assertEqual(
            ([(4, 5), (7, 8)], 7),
            CigarTools.compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------')
        )
