from mavis.interval import Interval
import networkx as nx
import itertools
from mavis.cluster.cluster import IntervalPair, merge_integer_intervals
import unittest


class TestIntervalPair(unittest.TestCase):
    def test_sets(self):
        i = Interval(1, 3)
        h = Interval(1, 4)
        s = set([IntervalPair(i, i), IntervalPair(h, h), IntervalPair(i, i)])
        self.assertEqual(2, len(s))

    def test__lt__(self):
        self.assertTrue(IntervalPair((1, 1), (1, 1)) < IntervalPair((1, 1), (1, 10)))
        self.assertFalse(IntervalPair((1, 10), (1, 1)) < IntervalPair((1, 1), (1, 10)))

    def test_merge(self):
        pairs = [IntervalPair((1, 2), (1, 10)), IntervalPair((1, 10), (2, 11)), IntervalPair((2, 11), (1, 2))]
        m = IntervalPair.merge(*pairs)
        s = merge_integer_intervals(*[p[0] for p in pairs])
        t = merge_integer_intervals(*[p[1] for p in pairs])
        self.assertEqual(s, m[0])
        self.assertEqual(t, m[1])

        pairs = [IntervalPair((1, 2), (1, 10)), IntervalPair((1, 10), (1, 10)), IntervalPair((2, 11), (1, 10))]
        m = IntervalPair.merge(*pairs)
        s = merge_integer_intervals(*[p[0] for p in pairs])
        t = merge_integer_intervals(*[p[1] for p in pairs])
        self.assertEqual(s, m[0])
        self.assertEqual(t, m[1])

    def test_abs_dist(self):
        x = IntervalPair((1, 1), (10, 11))
        y = IntervalPair((1, 1), (40, 41))
        self.assertEqual(14.5, IntervalPair.abs_dist(x, y))

    def test__redundant_maximal_kcliques(self):
        r = 10
        a = IntervalPair(Interval(1), Interval(100), id='a')
        b = IntervalPair(Interval(1), Interval(101), id='b')
        c = IntervalPair(Interval(10), Interval(90), id='c')
        d = IntervalPair(Interval(15), Interval(85), id='d')
        e = IntervalPair(Interval(40), Interval(45), id='e')
        G = nx.Graph()
        for n in [a, b, c, d, e]:
            G.add_node(n)
        for n1, n2 in itertools.combinations([a, b, c, d, e], 2):
            if IntervalPair.abs_dist(n1, n2) <= r:
                G.add_edge(n1, n2)
        self.assertTrue(G.has_edge(a, b))
        self.assertTrue(G.has_edge(a, c))
        self.assertTrue(G.has_edge(b, c))
        self.assertEqual(4, len(G.edges()))
        cliques = IntervalPair._redundant_maximal_kcliques(G)
        cliques = sorted([sorted(list(c)) for c in cliques])
        print('expect [a, b], [c, d], [e]')
        self.assertEqual([[a, b], [c, d], [e]], cliques)
        self.assertEqual(3, len(cliques))

        c = IntervalPair(Interval(6), Interval(94), id='c')
        G = nx.Graph()
        for n in [a, b, c, d, e]:
            G.add_node(n)
        for n1, n2 in itertools.combinations([a, b, c, d, e], 2):
            if IntervalPair.abs_dist(n1, n2) <= r:
                G.add_edge(n1, n2)

        self.assertTrue(G.has_edge(a, b))
        self.assertTrue(G.has_edge(a, c))
        self.assertTrue(G.has_edge(b, c))
        self.assertEqual(4, len(G.edges()))
        cliques = IntervalPair._redundant_maximal_kcliques(G)
        cliques = sorted([sorted(list(c)) for c in cliques])
        print('expect [a, b, c], [d], [e]')
        self.assertEqual([[a, b, c], [d], [e]], cliques)
        self.assertEqual(3, len(cliques))

    def test__redundant_maximal_kcliques_equidistant(self):
        r = 5
        a = IntervalPair(Interval(0), Interval(0), id='a')
        b = IntervalPair(Interval(10), Interval(10), id='b')
        c = IntervalPair(Interval(5), Interval(5), id='c')
        G = nx.Graph()
        for n in [a, b, c]:
            G.add_node(n)
        for n1, n2 in itertools.combinations([a, b, c], 2):
            if IntervalPair.abs_dist(n1, n2) <= r:
                G.add_edge(n1, n2)
        self.assertTrue(G.has_edge(b, c))
        self.assertTrue(G.has_edge(a, c))
        self.assertEqual(2, len(G.edges()))
        cliques = IntervalPair._redundant_maximal_kcliques(G)
        cliques = sorted([sorted(list(c)) for c in cliques])
        self.assertEqual([[a, c], [c, b]], cliques)
        self.assertEqual(2, len(cliques))

    def test__redundant_ordered_hierarchical_clustering(self):
        a = IntervalPair(Interval(1), Interval(1))
        b = IntervalPair(Interval(10), Interval(10))
        c = IntervalPair(Interval(15, 20), Interval(15, 20))
        d = IntervalPair(Interval(24), Interval(24))
        e = IntervalPair(Interval(33), Interval(33))
        interval_pairs = [
            set([a, b]),
            set([c]),
            set([d, e])
        ]
        groups = IntervalPair._redundant_ordered_hierarchical_clustering(interval_pairs, r=12)
        groups = sorted([sorted(list(c)) for c in groups])
        self.assertEqual(2, len(groups))
        self.assertEqual([a, b, c], groups[0])
        self.assertEqual([c, d, e], groups[1])


class TestMergeIntegerIntervals(unittest.TestCase):
    def test_varying_lengths(self):
        m = merge_integer_intervals((1, 2), (1, 9), (2, 10))
        self.assertEqual(Interval(1, 4), m)

    def test_same_length(self):
        m = merge_integer_intervals((1, 1), (10, 10))
        self.assertEqual(Interval(6), m)

    def test_empty_list_error(self):
        with self.assertRaises(AttributeError):
            merge_integer_intervals()

    def test_identical_even_length(self):
        m = merge_integer_intervals((1, 2), (1, 2), (1, 2))
        self.assertEqual(Interval(1, 2), m)

    def test_identical_odd_length(self):
        m = merge_integer_intervals((1, 3), (1, 3), (1, 3))
        self.assertEqual(Interval(1, 3), m)


if __name__ == "__main__":
    unittest.main()