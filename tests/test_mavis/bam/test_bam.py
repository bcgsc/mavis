import argparse
import warnings
from unittest import mock

import pytest
import timeout_decorator
from mavis.annotate.file_io import load_annotations, load_reference_genome
from mavis.bam import cigar as _cigar
from mavis.bam import read as _read
from mavis.bam.cache import BamCache
from mavis.bam.read import (
    breakpoint_pos,
    orientation_supports_type,
    read_pair_type,
    sequenced_strand,
)
from mavis.bam.stats import Histogram, compute_genome_bam_stats, compute_transcriptome_bam_stats
from mavis.constants import CIGAR, DNA_ALPHABET, ORIENT, READ_PAIR_TYPE, STRAND, SVTYPE
from mavis.interval import Interval

from ...util import get_data
from ..mock import MockBamFileHandle, MockRead

REFERENCE_GENOME = None


def setUpModule():
    warnings.simplefilter('ignore')
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')


class TestBamCache:
    def test___init__(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        assert b.fh == fh

    def test_add_read(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        r = mock.MagicMock(query_name='name', query_sequence='')
        b.add_read(r)
        assert len(b.cache.values()) == 1
        b.add_read(r)
        assert len(b.cache.values()) == 1
        r.reference_start = 0
        b.add_read(r)
        assert len(b.cache.values()) == 1

    @mock.patch('mavis.util.logger')
    def test_add_invalid_read(self, log_patcher):
        bad_read = mock.Mock(
            is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ'
        )
        cache = BamCache(MockBamFileHandle())
        cache.add_read(bad_read)
        assert len(cache.cache) == 0
        log_patcher.method_calls[0].assert_called_with('ignoring invalid read: BAD_READ')

    @mock.patch('mavis.util.logger')
    def test_fetch_invalid_read(self, log_patcher):
        bad_read = mock.Mock(
            is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ'
        )
        fh = mock.Mock(references=['chr'], spec=['references', 'fetch'])
        fh.configure_mock(**{'fetch.return_value': [bad_read]})
        cache = BamCache(fh)
        cache.fetch('chr', 1, 10)
        assert len(cache.cache) == 0
        log_patcher.method_calls[0].assert_called_with('ignoring invalid read: BAD_READ')

    @mock.patch('mavis.util.logger')
    def test_bin_fetch_invalid_read(self, log_patcher):
        bad_read = mock.Mock(
            is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ'
        )
        fh = mock.Mock(references=['chr'], spec=['references', 'fetch'])
        fh.configure_mock(**{'fetch.return_value': [bad_read]})
        cache = BamCache(fh)
        cache.fetch_from_bins('chr', 1, 10)
        assert len(cache.cache) == 0
        log_patcher.method_calls[0].assert_called_with('ignoring invalid read: BAD_READ')

    def test_reference_id(self):
        fh = MockBamFileHandle({'1': 0})
        b = BamCache(fh)
        assert b.reference_id('1') == 0
        with pytest.raises(KeyError):
            b.reference_id('2')

    def test_get_read_reference_name(self):
        fh = MockBamFileHandle({'1': 0})
        b = BamCache(fh)
        r = MockRead('name', 0)
        assert b.get_read_reference_name(r) == '1'

    def test_generate_fetch_bins_single(self):
        assert BamCache._generate_fetch_bins(1, 100, 1, 1) == [(1, 100)]

    def test_generate_fetch_bins_multi(self):
        assert BamCache._generate_fetch_bins(1, 100, 2, 1) == [(1, 50), (51, 100)]
        assert BamCache._generate_fetch_bins(1, 100, 5, 1) == [
            (1, 20),
            (21, 40),
            (41, 60),
            (61, 80),
            (81, 100),
        ]

    def test_generate_fetch_bins_large_min_size(self):
        assert BamCache._generate_fetch_bins(1, 100, 5, 50) == [(1, 50), (51, 100)]

    def test_fetch_single_read(self):
        b = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
        s = b.fetch_from_bins('reference3', 1382, 1383, read_limit=1, sample_bins=1)
        assert len(s) == 1
        r = list(s)[0]
        assert r.qname == 'HISEQX1_11:4:2122:14275:37717:split'
        b.close()

    def test_get_mate(self):
        # dependant on fetch working
        b = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
        s = b.fetch_from_bins('reference3', 1382, 1383, read_limit=1, sample_bins=1)
        assert len(s) == 1
        r = list(s)[0]
        assert r.qname == 'HISEQX1_11:4:2122:14275:37717:split'
        o = b.get_mate(r, allow_file_access=True)
        assert len(o) == 1
        assert o[0].qname == 'HISEQX1_11:4:2122:14275:37717:split'


class TestModule:
    """
    test class for functions in the validate namespace
    that are not associated with a class
    """

    def test_alphabet_matching(self):
        assert DNA_ALPHABET.match('N', 'A')
        assert DNA_ALPHABET.match('A', 'N')

    def test_breakpoint_pos(self):
        # ==========+++++++++>
        r = MockRead(reference_start=10, cigar=[(CIGAR.M, 10), (CIGAR.S, 10)])
        assert _read.breakpoint_pos(r) == 19

        with pytest.raises(AttributeError):
            breakpoint_pos(r, ORIENT.RIGHT)

        assert _read.breakpoint_pos(r, ORIENT.LEFT) == 19

        # ++++++++++=========>
        r = MockRead(reference_start=10, cigar=[(CIGAR.S, 10), (CIGAR.M, 10)])
        assert _read.breakpoint_pos(r) == 10

        with pytest.raises(AttributeError):
            breakpoint_pos(r, ORIENT.LEFT)

        assert _read.breakpoint_pos(r, ORIENT.RIGHT) == 10

        with pytest.raises(AttributeError):
            r = MockRead(reference_start=10, cigar=[(CIGAR.X, 10), (CIGAR.M, 10)])
            _read.breakpoint_pos(r, ORIENT.LEFT)

    def test_nsb_align(self):
        ref = (
            'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG'
            'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        )
        seq = (
            'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC'
            'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        )
        _read.nsb_align(ref, seq)
        # GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG


class TestNsbAlign:
    def test_length_seq_le_ref(self):
        ref = (
            'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG'
            'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        )
        seq = (
            'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC'
            'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        )
        alignment = _read.nsb_align(ref, seq)
        assert len(alignment) == 1
        alignment = _read.nsb_align(ref, seq, min_consecutive_match=20)
        assert len(alignment) == 0

    def test_length_ref_le_seq(self):
        pass

    def test_length_ref_eq_seq(self):
        pass

    @timeout_decorator.timeout(5)
    def test_long_ref_seq(self):
        ref = str(REFERENCE_GENOME['test_bam_long_ref'].seq)
        seq = (
            'TGAGGTCAGGAGTTTGAGACCAGCCTGGACAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAAATTAGCCAGGCATGGTGGTGGATGCCTGTAAT'
            'CGCAGCTACTCAGGAGATCGGAAG'
        )
        alignment = _read.nsb_align(ref, seq, min_consecutive_match=6)
        assert len(alignment) == 1

    def test_left_softclipping(self):
        ref = 'TAAGCTTCTTCCTTTTTCTATGCCACCTACATAGGCATTTTGCATGGTCAGATTGGAATTTACATAATGCATACATGCAAAGAAATATATAGAAGCCAGATATATAAGGTAGTACATTGGCAGGCTTCATATATATAGACTCCCCCATATTGTCTATATGCTAAAAAAGTATTTTAAATCCTTAAATTTTATTTTTGTTCTCTGCATTTGAAATCTTTATCAACTAGGTCATGAAAATAGCCAGTCGGTTCTCCTTTTGGTCTATTAGAATAAAATCTGGACTGCAACTGAGAAGCAGAAGGTAATGTCAGAATGTAT'
        seq = 'GCTAAAAAAGTATTTTAAATCCTTAAATGTTATTTTTGTTCTC'
        alignments = _read.nsb_align(ref, seq, min_consecutive_match=6)
        assert len(alignments) == 1
        print(alignments)
        seq = 'CTTATAAAGCTGGAGTATCTGCTGAGAGCATCAGGAATTGACATCTAGGATAATGAGAGAAGGCTGATCATGGACAACATATAGCCTTTCTAGTAGATGCAGCTGAGGCTAAAAAAGTATTTTAAATCCTTAAATGTTATTTTTGTTCTC'
        alignments = _read.nsb_align(ref, seq, min_consecutive_match=6, min_overlap_percent=0.5)
        assert len(alignments) == 1

    def test_min_overlap(self):
        ref = 'ATTACATTAAAGATTCAAACTCCTAGAGTTTTTTTGATTTTTAGTATGATCTTTAGATAAAAAAAAAGGAAGAAAAAGAAAAAAAAACAGAGTCTATTAAGGCATCTTCTATGGTCAGATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGAT'
        seq = 'GATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTGTATATATATATATATATATATATATATATATACCCAGTTTCAAGCAG'
        alignments = _read.nsb_align(
            ref,
            seq,
            min_consecutive_match=15,
            min_match=0.95,
            min_overlap_percent=(len(seq) - 15) / len(seq),
        )
        print(alignments)
        assert len(alignments) == 0


@pytest.fixture
def stranded_reads():
    n = argparse.Namespace()
    n.read1_pos_neg = MockRead(is_reverse=False, is_read1=True, mate_is_reverse=True)
    assert not n.read1_pos_neg.is_read2
    n.read1_neg_pos = MockRead(is_reverse=True, is_read1=True, mate_is_reverse=False)
    n.read1_pos_pos = MockRead(is_reverse=False, is_read1=True, mate_is_reverse=False)
    n.read1_neg_neg = MockRead(is_reverse=True, is_read1=True, mate_is_reverse=True)

    n.read2_pos_neg = MockRead(is_reverse=True, is_read1=False, mate_is_reverse=True)
    assert n.read2_pos_neg.is_read2
    n.read2_neg_pos = MockRead(is_reverse=False, is_read1=False, mate_is_reverse=False)
    n.read2_pos_pos = MockRead(is_reverse=False, is_read1=False, mate_is_reverse=False)
    n.read2_neg_neg = MockRead(is_reverse=True, is_read1=False, mate_is_reverse=True)

    n.unpaired_pos = MockRead(is_reverse=False, is_paired=False)
    n.unpaired_neg = MockRead(is_reverse=True, is_paired=False)
    return n


class TestReadPairStrand:
    def test_read_pair_strand_det1_read1(self, stranded_reads):
        assert (
            sequenced_strand(stranded_reads.read1_pos_neg, strand_determining_read=1) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read1_neg_pos, strand_determining_read=1) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read1_pos_pos, strand_determining_read=1) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read1_neg_neg, strand_determining_read=1) == STRAND.NEG
        )

    def test_read_pair_strand_det1_read2(self, stranded_reads):
        assert (
            sequenced_strand(stranded_reads.read2_pos_neg, strand_determining_read=1) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read2_neg_pos, strand_determining_read=1) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read2_pos_pos, strand_determining_read=1) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read2_neg_neg, strand_determining_read=1) == STRAND.POS
        )

    def test_read_pair_strand_det2_read2(self, stranded_reads):
        assert (
            sequenced_strand(stranded_reads.read2_pos_neg, strand_determining_read=2) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read2_neg_pos, strand_determining_read=2) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read2_pos_pos, strand_determining_read=2) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read2_neg_neg, strand_determining_read=2) == STRAND.NEG
        )

    def test_read_pair_strand_det2_read1(self, stranded_reads):
        assert (
            sequenced_strand(stranded_reads.read1_pos_neg, strand_determining_read=2) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read1_neg_pos, strand_determining_read=2) == STRAND.POS
        )
        assert (
            sequenced_strand(stranded_reads.read1_pos_pos, strand_determining_read=2) == STRAND.NEG
        )
        assert (
            sequenced_strand(stranded_reads.read1_neg_neg, strand_determining_read=2) == STRAND.POS
        )

    def test_read_pair_strand_unpaired(self, stranded_reads):
        with pytest.raises(ValueError):
            sequenced_strand(stranded_reads.unpaired_pos)
        with pytest.raises(ValueError):
            sequenced_strand(stranded_reads.unpaired_neg)

    def test_read_pair_strand_det_error(self, stranded_reads):
        with pytest.raises(ValueError):
            sequenced_strand(stranded_reads.read1_pos_neg, strand_determining_read=3)


@pytest.fixture
def read_pairs():
    n = argparse.Namespace()
    n.LR = MockRead(
        reference_id=0,
        next_reference_id=0,
        reference_start=1,
        next_reference_start=2,
        is_reverse=False,
        mate_is_reverse=True,
    )
    n.LL = MockRead(
        reference_id=0,
        next_reference_id=0,
        reference_start=1,
        next_reference_start=2,
        is_reverse=False,
        mate_is_reverse=False,
    )
    n.RR = MockRead(
        reference_id=0,
        next_reference_id=0,
        reference_start=1,
        next_reference_start=2,
        is_reverse=True,
        mate_is_reverse=True,
    )
    n.RL = MockRead(
        reference_id=0,
        next_reference_id=0,
        reference_start=1,
        next_reference_start=2,
        is_reverse=True,
        mate_is_reverse=False,
    )
    return n


class TestReadPairType:
    def test_read_pair_type_LR(self, read_pairs):
        assert read_pair_type(read_pairs.LR) == READ_PAIR_TYPE.LR

    def test_read_pair_type_LL(self, read_pairs):
        assert read_pair_type(read_pairs.LL) == READ_PAIR_TYPE.LL

    def test_read_pair_type_RR(self, read_pairs):
        assert read_pair_type(read_pairs.RR) == READ_PAIR_TYPE.RR

    def test_read_pair_type_RL(self, read_pairs):
        assert read_pair_type(read_pairs.RL) == READ_PAIR_TYPE.RL

    def test_orientation_supports_type_deletion(self, read_pairs):
        assert orientation_supports_type(read_pairs.LR, SVTYPE.DEL)
        assert not orientation_supports_type(read_pairs.RL, SVTYPE.DEL)
        assert not orientation_supports_type(read_pairs.LL, SVTYPE.DEL)
        assert not orientation_supports_type(read_pairs.RR, SVTYPE.DEL)

    def test_orientation_supports_type_insertion(self, read_pairs):
        assert orientation_supports_type(read_pairs.LR, SVTYPE.INS)
        assert not orientation_supports_type(read_pairs.RL, SVTYPE.INS)
        assert not orientation_supports_type(read_pairs.LL, SVTYPE.INS)
        assert not orientation_supports_type(read_pairs.RR, SVTYPE.INS)

    def test_orientation_supports_type_inversion(self, read_pairs):
        assert not orientation_supports_type(read_pairs.LR, SVTYPE.INV)
        assert not orientation_supports_type(read_pairs.RL, SVTYPE.INV)
        assert orientation_supports_type(read_pairs.LL, SVTYPE.INV)
        assert orientation_supports_type(read_pairs.RR, SVTYPE.INV)

    def test_orientation_supports_type_translocation_inversion(self, read_pairs):
        assert not orientation_supports_type(read_pairs.LR, SVTYPE.ITRANS)
        assert not orientation_supports_type(read_pairs.RL, SVTYPE.ITRANS)
        assert orientation_supports_type(read_pairs.LL, SVTYPE.ITRANS)
        assert orientation_supports_type(read_pairs.RR, SVTYPE.ITRANS)

    def test_orientation_supports_type_trans_duplication(self, read_pairs):
        assert not orientation_supports_type(read_pairs.LR, SVTYPE.DUP)
        assert orientation_supports_type(read_pairs.RL, SVTYPE.DUP)
        assert not orientation_supports_type(read_pairs.LL, SVTYPE.DUP)
        assert not orientation_supports_type(read_pairs.RR, SVTYPE.DUP)

    def test_orientation_supports_type_translocation(self, read_pairs):
        assert orientation_supports_type(read_pairs.LR, SVTYPE.TRANS)
        assert orientation_supports_type(read_pairs.RL, SVTYPE.TRANS)
        assert not orientation_supports_type(read_pairs.LL, SVTYPE.TRANS)
        assert not orientation_supports_type(read_pairs.RR, SVTYPE.TRANS)


class TestHistogram:
    def test_add(self):
        h = Histogram()
        h.add(1)
        h.add(1)
        assert h[1] == 2
        h.add(1, 4)
        assert h[1] == 6

    def test_median(self):
        h = Histogram()
        for i in range(1, 11):
            h.add(i)
        assert h.median() == 5.5
        h.add(11)
        assert h.median() == 6

    def test_distib_stderr(self):
        h = Histogram()
        for i in range(0, 11):
            h.add(i)
        for i in range(4, 8):
            h.add(i)
        m = h.median()
        assert m == 5
        err = h.distribution_stderr(m, 1)
        assert err == 116 / 15

    def test_add_operator(self):
        x = Histogram()
        y = Histogram()
        x.add(1)
        y.add(1, 4)
        z = x + y
        assert x[1] == 1
        assert y[1] == 4
        assert z[1] == 5


class TestBamStats:
    def test_genome_bam_stats(self):
        bamfh = BamCache(get_data('mock_reads_for_events.sorted.bam'))
        stats = compute_genome_bam_stats(
            bamfh, 1000, 100, min_mapping_quality=1, sample_cap=10000, distribution_fraction=0.99
        )
        assert 50 >= abs(stats.median_fragment_size - 420)
        assert stats.read_length == 150
        bamfh.close()

    def test_trans_bam_stats(self):
        bamfh = BamCache(get_data('mock_trans_reads_for_events.sorted.bam'))
        annotations = load_annotations(get_data('mock_annotations.json'))
        stats = compute_transcriptome_bam_stats(
            bamfh,
            annotations,
            100,
            min_mapping_quality=1,
            stranded=True,
            sample_cap=10000,
            distribution_fraction=0.99,
        )
        assert abs(stats.median_fragment_size - 185) < 5
        assert stats.read_length == 75
        assert stats.stdev_fragment_size < 50
        bamfh.close()


@pytest.fixture
def contig_read():
    return MockRead(
        cigar=_cigar.convert_string_to_cigar('275M18I12041D278M'),
        reference_start=89700025,
        reference_name='10',
    )


class TestMapRefRangeToQueryRange:
    def test_full_aligned_portion(self, contig_read):
        ref_range = Interval(89700026, 89712619)
        qrange = _read.map_ref_range_to_query_range(contig_read, ref_range)
        assert len(qrange) == 571
        assert qrange.start == 1
        assert qrange.end == 571

    def test_multiple_events(self, contig_read):
        ref_range = Interval(89700067, 89712347)
        qrange = _read.map_ref_range_to_query_range(contig_read, ref_range)
        assert len(qrange) == len(ref_range) - 12041 + 18

    def test_no_events(self, contig_read):
        ref_range = Interval(89700031, 89700040)
        qrange = _read.map_ref_range_to_query_range(contig_read, ref_range)
        assert len(qrange) == 10
        assert qrange.start == 6
        assert qrange.end == 15
