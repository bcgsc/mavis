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

    def test_kmers(self):
        k = kmers('ABCDEFG', 2)
        self.assertEqual(['AB', 'BC', 'CD', 'DE', 'EF', 'FG'], k)
        k = kmers('ABCDEFG', 3)
        self.assertEqual(['ABC', 'BCD', 'CDE', 'DEF', 'EFG'], k)

    def test_assemble(self):
        sequences = ['ABCD', 'BCDE', 'CDEF', 'ABCDE', 'DEFG']
        c = assemble(sequences, min_edge_weight=1)
        self.assertEqual(1, len(c))
        self.assertEqual('ABCDEFG', c[0].seq)

    def test_assemble_empty_list(self):
        self.assertEqual([], assemble([]))


class TestDeBruijnGraph(unittest.TestCase):

    def test_trim_low_weight_tails_forks(self):
        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(8, 7)
        G.add_edge(9, 8)
        G.trim_low_weight_tails(2)
        self.assertEqual([1, 2, 3, 4, 5, 6], sorted(G.nodes()))
        
        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(7, 8)
        G.add_edge(8, 7)
        G.add_edge(9, 8)
        G.trim_low_weight_tails(2)
        self.assertEqual([1, 2, 3, 4, 5, 6, 7, 8], sorted(G.nodes()))
        
        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(7, 8)
        G.add_edge(9, 8)
        G.trim_low_weight_tails(2)
        self.assertEqual([1, 2, 3, 4, 5, 6], sorted(G.nodes()))

    def test_add_edge(self):
        G = DeBruijnGraph()
        G.add_edge(1, 2)
        self.assertEqual(1, G.edge_freq[(1, 2)])
        G.add_edge(1, 2)
        self.assertEqual(2, G.edge_freq[(1, 2)])
        G.add_edge(1, 2, 5)
        self.assertEqual(7, G.edge_freq[(1, 2)])


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
