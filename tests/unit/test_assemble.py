from mavis.constants import *
from mavis.assemble import *
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
        c = assemble(
            sequences, assembly_min_edge_weight=0, assembly_min_nc_edge_weight=1, assembly_min_exact_match_to_remap=1)
        self.assertEqual(1, len(c))
        self.assertEqual('ABCDEFG', c[0].seq)
        self.assertEqual(5, c[0].remap_score())

    def test_assemble_empty_list(self):
        self.assertEqual([], assemble([]))

    def test_repeat_region_assembly(self):
        rep = 'ABCDEF'
        seqs = kmers(rep + rep, len(rep))
        contigs = assemble(seqs, assembly_min_edge_weight=1, assembly_min_exact_match_to_remap=1)
        self.assertEqual(0, len(contigs))


class TestDeBruijnGraph(unittest.TestCase):

    def test_trim_tails_by_freq_forks(self):
        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(8, 7)
        G.add_edge(9, 8)
        G.trim_tails_by_freq(2)
        self.assertEqual([1, 2, 3, 4, 5, 6], sorted(G.nodes()))

        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(7, 8)
        G.add_edge(8, 7)
        G.add_edge(9, 8)
        G.trim_tails_by_freq(2)
        self.assertEqual([1, 2, 3, 4, 5, 6, 7, 8], sorted(G.nodes()))

        G = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            G.add_edge(s, t)
        G.add_node(10)  # singlet
        G.add_edge(7, 6)
        G.add_edge(7, 8)
        G.add_edge(9, 8)
        G.trim_tails_by_freq(2)
        self.assertEqual([1, 2, 3, 4, 5, 6], sorted(G.nodes()))

    def test_add_edge(self):
        G = DeBruijnGraph()
        G.add_edge(1, 2)
        self.assertEqual(1, G.get_edge_freq(1, 2))
        G.add_edge(1, 2)
        self.assertEqual(2, G.get_edge_freq(1, 2))
        G.add_edge(1, 2, 5)
        self.assertEqual(7, G.get_edge_freq(1, 2))

    def test_trim_noncutting_paths_by_freq_degree_stop(self):
        g = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4], 2):
            g.add_edge(s, t, freq=4)
        for s, t in itertools.combinations([5, 6, 7, 8], 2):
            g.add_edge(s, t, freq=4)
        path1 = [5, 9, 10, 11, 12, 1]
        for s, t in zip(path1, path1[1:]):
            g.add_edge(s, t)
        for edge in g.edges():
            print(edge)
        g.trim_noncutting_paths_by_freq(3)
        self.assertEqual(list(range(1, 9)) + path1[1:-1], g.nodes())

        # add an equal weight path to force namesorting
        path2 = [5, 13, 14, 15, 16, 1]
        for s, t in zip(path2, path2[1:]):
            g.add_edge(s, t)

        g.trim_noncutting_paths_by_freq(3)
        self.assertEqual(list(range(1, 9)) + path2[1:-1], g.nodes())

        # add back the original path with a higher (but still low) weight
        for s, t in zip(path1, path1[1:]):
            g.add_edge(s, t, freq=2)

        g.trim_noncutting_paths_by_freq(3)
        self.assertEqual(list(range(1, 9)) + path1[1:-1], g.nodes())

        # add the second path with 1 high weight edge
        path2 = [5, 13, 14, 15, 16, 1]
        for s, t in zip(path2, path2[1:]):
            g.add_edge(s, t)
        g.add_edge(14, 15, freq=6)

        g.trim_noncutting_paths_by_freq(3)
        self.assertEqual(list(range(1, 9)) + path2[1:-1], g.nodes())
