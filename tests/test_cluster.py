from structural_variant.interval import Interval
import networkx as nx
import itertools
import structural_variant.cluster
from structural_variant.cluster import IntervalPair
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.constants import STRAND

import unittest

REF_CHR='fake'

class TestIntervalPair(unittest.TestCase):
    def test_sets(self):
        i = Interval(1, 3)
        h = Interval(1, 4)
        s = set([IntervalPair(i, i), IntervalPair(h, h), IntervalPair(i, i)])
        self.assertEqual(2, len(s))

    def test__lt__(self):
        self.assertTrue(IntervalPair((1,1),(1,1)) < IntervalPair((1,1),(1,10)))
        self.assertFalse(IntervalPair((1,10),(1,1)) < IntervalPair((1,1),(1,10)))

    def test_weighted_mean(self):
        pairs = [IntervalPair((1, 2), (1, 10)), IntervalPair((1, 10), (2, 11)), IntervalPair((2, 11), (1, 2))]
        m = IntervalPair.weighted_mean(*pairs)
        self.assertEqual(IntervalPair(
            Interval.weighted_mean(*[p.start for p in pairs]), Interval.weighted_mean(*[p.end for p in pairs])), m)

        pairs = [IntervalPair((1, 2), (1, 10)), IntervalPair((1, 10), (1, 10)), IntervalPair((2, 11), (1, 10))]
        m = IntervalPair.weighted_mean(*pairs)
        self.assertEqual(IntervalPair(
            Interval.weighted_mean(*[p.start for p in pairs]), Interval.weighted_mean(*[p.end for p in pairs])), m)

    def test_dist(self):
        x = IntervalPair((1, 1), (10, 11))
        y = IntervalPair((1, 1), (40, 41))
        self.assertEqual(15, x.dist(y))

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
            if n1.dist(n2) <= r:
                G.add_edge(n1, n2)
        self.assertTrue(G.has_edge(a, b))
        self.assertTrue(G.has_edge(a, c))
        self.assertTrue(G.has_edge(b, c))
        self.assertEqual(4, len(G.edges()))
        cliques = IntervalPair._redundant_maximal_kcliques(G)
        cliques = sorted([sorted(list(c)) for c in cliques])
        self.assertEqual([[a, b], [c, d], [e]], cliques)
        self.assertEqual(3, len(cliques))

        c = IntervalPair(Interval(6), Interval(94), id='c')
        G = nx.Graph()
        for n in [a, b, c, d, e]:
            G.add_node(n)
        for n1, n2 in itertools.combinations([a, b, c, d, e], 2):
            if n1.dist(n2) <= r:
                G.add_edge(n1, n2)

        self.assertTrue(G.has_edge(a, b))
        self.assertTrue(G.has_edge(a, c))
        self.assertTrue(G.has_edge(b, c))
        self.assertEqual(4, len(G.edges()))
        cliques = IntervalPair._redundant_maximal_kcliques(G)
        cliques = sorted([sorted(list(c)) for c in cliques])
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
            if n1.dist(n2) <= r:
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

    @unittest.skip("TODO")
    def test_cluster_breakpoint_pairs(self):
        bpp1 = BreakpointPair(Breakpoint(REF_CHR, 31),
                              Breakpoint(REF_CHR, 129),
                              opposing_strands=False)
        bpp2 = BreakpointPair(Breakpoint(REF_CHR, 32),
                              Breakpoint(REF_CHR, 120),
                              opposing_strands=False)
        bpp3 = BreakpointPair(Breakpoint(REF_CHR, 40),
                              Breakpoint(REF_CHR, 139),
                              opposing_strands=False)
        bpp4 = BreakpointPair(Breakpoint(REF_CHR, 43),
                              Breakpoint(REF_CHR, 151),
                              opposing_strands=False)
        bpp4 = BreakpointPair(Breakpoint(REF_CHR, 45),
                              Breakpoint(REF_CHR, 146),
                              opposing_strands=False)
        bpp5 = BreakpointPair(Breakpoint(REF_CHR, 50),
                              Breakpoint(REF_CHR, 150),
                              opposing_strands=False)
        from pprint import pprint
        I = {bpp1, bpp2, bpp3, bpp4, bpp5}
        print(I)
        groups = structural_variant.cluster.cluster_breakpoint_pairs(I, 9, 4)
#        groups = sorted([sorted(list(c)) for c in groups])
        pprint(groups)
        groups2 = sorted(groups.keys())
        print(groups2)
        for key in groups.keys():
            print(str(key))
            print("grouping")
            for i in groups[key]:
                print(str(i))
            print("done")
        print(str(groups[list(groups.keys())[0]][0]))
        self.assertTrue(False)


if __name__ == "__main__":
    unittest.main()
