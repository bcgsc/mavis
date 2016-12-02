from structural_variant.constants import *
from structural_variant.assemble import *
import unittest


class TestModule(unittest.TestCase):
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
        self.assertEqual(5, c[0].remap_score())

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
