import itertools
import random

import pytest
from mavis.constants import DNA_ALPHABET
from mavis.validate.assemble import Contig, DeBruijnGraph, assemble, filter_contigs, kmers

from ..util import get_data, long_running_test


class TestModule:
    """
    test class for functions in the validate namespace
    that are not associated with a class
    """

    def test_alphabet_matching(self):
        assert DNA_ALPHABET.match('N', 'A')
        assert DNA_ALPHABET.match('A', 'N')

    def test_kmers(self):
        k = kmers('ABCDEFG', 2)
        assert k == ['AB', 'BC', 'CD', 'DE', 'EF', 'FG']
        k = kmers('ABCDEFG', 3)
        assert k == ['ABC', 'BCD', 'CDE', 'DEF', 'EFG']

    def test_assemble(self):
        sequences = ['ABCD', 'BCDE', 'CDEF', 'ABCDE', 'DEFG']
        c = assemble(sequences, 3, min_edge_trim_weight=1, remap_min_exact_match=1)
        assert len(c) == 1
        assert c[0].seq == 'ABCDEFG'
        assert c[0].remap_score() == 5

    def test_assemble_empty_list(self):
        assert assemble([], 1) == []

    def test_repeat_region_assembly(self):
        rep = 'ABCDEF'
        seqs = kmers(rep + rep, len(rep))
        contigs = assemble(seqs, len(rep) - 1, remap_min_exact_match=1)
        assert len(contigs) == 0


class TestFilterContigs:
    def test_drop_reverse_complement(self):
        c1 = Contig('atcgatcgatcgatcgatcgatcgatatagggcatcagc', 1)
        c2 = Contig('gctgatgccctatatcgatcgatcgatcgatcgatcgat', 1)
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 1
        assert result[0].seq == c1.seq

    def test_drop_alt_allele_alphabetically(self):
        c1 = Contig('atcgatcgatcgatcgatcgatcgatatagggcatcagc', 1)
        c2 = Contig('atcgatcgatcgatcgatctatcgatatagggcatcagc', 1)
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 1
        assert result[0].seq == c1.seq

    def test_drop_alt_allele_by_score(self):
        c1 = Contig('atcgatcgatcgatcgatcgatcgatatagggcatcagc', 2)
        c2 = Contig('atcgatcgatcgatcgatctatcgatatagggcatcagc', 1)
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 1
        assert result[0].seq == c1.seq

    def test_retain_disimilar(self):
        c1 = Contig('atcgatcgatcgatcgatcgatcgatatagggcatcagc', 2)
        c2 = Contig('atcgadatcgatcgatcgatctgtdstcgatatagggca', 1)
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 2

    def test_retain_disimilar_different_lengths(self):
        c1 = Contig('atcgatcgatcgatcgatcgatcgatatagggcatcagc', 2)
        c2 = Contig('atcgatcgatcgatcgatcgatcccgtgatatagggcatcagc', 1)
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 2

    def test_drop_similar_different_lengths(self):
        c1 = Contig(
            'atcgatcgatcgatcgatatcgatcgatcgatcgatatcgatcgatcgatcgatatcgatcgatcgatcgatcgatcgatatgggcatcagc',
            2,
        )
        c2 = Contig(
            'atcgatcgatcgatcgatatcgatcgatcgatcgatatcgatcgatcgatcgatatcgatcgatcgatcgatcgatcgatagggcatcagc',
            1,
        )
        result = filter_contigs([c2, c1], 0.10)
        assert len(result) == 1
        assert result[0].seq == c1.seq


class TestDeBruijnGraph:
    @long_running_test
    def test_trim_tails_by_freq_forks(self):
        g = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            g.add_edge(s, t)
            print('s => t', s, t)
        g.add_edge(6, 1)
        g.add_node(10)  # singlet
        g.add_edge(7, 6)
        g.add_edge(8, 7)
        g.add_edge(9, 8)
        g.trim_tails_by_freq(2)
        assert sorted(g.get_nodes()) == [1, 2, 3, 4, 5, 6]

        g = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            g.add_edge(s, t)
        g.add_node(10)  # singlet
        g.add_edge(6, 1)
        g.add_edge(7, 6)
        g.add_edge(7, 8)
        g.add_edge(8, 7)
        g.add_edge(9, 8)
        g.trim_tails_by_freq(2)
        assert sorted(g.get_nodes()) == [1, 2, 3, 4, 5, 6, 7, 8]

        g = DeBruijnGraph()
        for s, t in itertools.combinations([1, 2, 3, 4, 5, 6], 2):
            g.add_edge(s, t)
        g.add_node(10)  # singlet
        g.add_edge(6, 1)
        g.add_edge(7, 6)
        g.add_edge(7, 8)
        g.add_edge(9, 8)
        g.trim_tails_by_freq(2)
        assert sorted(g.get_nodes()) == [1, 2, 3, 4, 5, 6]

    def test_add_edge(self):
        g = DeBruijnGraph()
        g.add_edge(1, 2)
        assert g.get_edge_freq(1, 2) == 1
        g.add_edge(1, 2)
        assert g.get_edge_freq(1, 2) == 2
        g.add_edge(1, 2, 5)
        assert g.get_edge_freq(1, 2) == 7

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
        print('g.nodes', g.nodes)
        assert g.get_nodes() == list(range(1, 9)) + path1[1:-1]
        print('g.nodes', g.nodes)

        # add an equal weight path to force namesorting
        path2 = [5, 13, 14, 15, 16, 1]
        for s, t in zip(path2, path2[1:]):
            g.add_edge(s, t)
        print('g.nodes', g.nodes)
        g.trim_noncutting_paths_by_freq(3)
        print('g.nodes', g.nodes)
        assert g.get_nodes() == list(range(1, 9)) + path2[1:-1]

        # add back the original path with a higher (but still low) weight
        for s, t in zip(path1, path1[1:]):
            g.add_edge(s, t, freq=2)

        g.trim_noncutting_paths_by_freq(3)
        assert g.get_nodes() == list(range(1, 9)) + path1[1:-1]

        # add the second path with 1 high weight edge
        path2 = [5, 13, 14, 15, 16, 1]
        for s, t in zip(path2, path2[1:]):
            g.add_edge(s, t)
        g.add_edge(14, 15, freq=6)

        g.trim_noncutting_paths_by_freq(3)
        assert g.get_nodes() == list(range(1, 9)) + path2[1:-1]


@pytest.fixture
def assembly_sequences():
    # load the sequences
    with open(get_data('test_assembly_sequences.txt')) as fh:
        seq = [i.strip() for i in fh.readlines()]
    return seq


class TestFullAssemly:
    @long_running_test
    def test_deterministic_assembly(self, assembly_sequences):
        contig_sequences = set()
        for i in range(20):
            random.shuffle(assembly_sequences)
            contigs = assemble(
                assembly_sequences,
                111,
                min_edge_trim_weight=3,
                assembly_max_paths=8,
                assembly_min_uniq=0.1,
                min_complexity=0.1,
            )
            assert len(contigs) == 1
            contig_sequences.add(contigs[0].seq)
        assert len(contig_sequences) == 1
