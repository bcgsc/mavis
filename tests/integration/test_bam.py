import logging
import os
import unittest
from unittest import mock
import warnings

from mavis.annotate.file_io import load_reference_genes, load_reference_genome
from mavis.bam import cigar as _cigar
from mavis.bam import read as _read
from mavis.bam.cache import BamCache
from mavis.bam.read import breakpoint_pos, orientation_supports_type, read_pair_type, sequenced_strand
from mavis.bam.stats import compute_genome_bam_stats, compute_transcriptome_bam_stats, Histogram
from mavis.constants import CIGAR, DNA_ALPHABET, ORIENT, READ_PAIR_TYPE, STRAND, SVTYPE, NA_MAPPING_QUALITY
from mavis.interval import Interval
import timeout_decorator

from . import MockRead, MockBamFileHandle
from ..util import get_data


REFERENCE_GENOME = None


def setUpModule():
    warnings.simplefilter('ignore')
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestBamCache(unittest.TestCase):

    def test___init__(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        self.assertEqual(fh, b.fh)

    def test_add_read(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        r = mock.MagicMock(query_name='name', query_sequence='')
        b.add_read(r)
        self.assertEqual(1, len(b.cache.values()))
        b.add_read(r)
        self.assertEqual(1, len(b.cache.values()))
        r.reference_start = 0
        b.add_read(r)
        self.assertEqual(1, len(b.cache.values()))

    @mock.patch('mavis.util.LOG')
    def test_add_invalid_read(self, log_patcher):
        bad_read = mock.Mock(is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ')
        cache = BamCache(MockBamFileHandle())
        cache.add_read(bad_read)
        self.assertEqual(0, len(cache.cache))
        log_patcher.assert_called_with('ignoring invalid read', 'BAD_READ', level=logging.DEBUG)

    @mock.patch('mavis.util.LOG')
    def test_fetch_invalid_read(self, log_patcher):
        bad_read = mock.Mock(is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ')
        fh = mock.Mock(references=['chr'], spec=['references', 'fetch'])
        fh.configure_mock(**{'fetch.return_value': [bad_read]})
        cache = BamCache(fh)
        cache.fetch('chr', 1, 10)
        self.assertEqual(0, len(cache.cache))
        log_patcher.assert_called_with('ignoring invalid read', 'BAD_READ', level=logging.DEBUG)

    @mock.patch('mavis.util.LOG')
    def test_bin_fetch_invalid_read(self, log_patcher):
        bad_read = mock.Mock(is_unmapped=False, reference_start=0, reference_end=0, query_name='BAD_READ')
        fh = mock.Mock(references=['chr'], spec=['references', 'fetch'])
        fh.configure_mock(**{'fetch.return_value': [bad_read]})
        cache = BamCache(fh)
        cache.fetch_from_bins('chr', 1, 10)
        self.assertEqual(0, len(cache.cache))
        log_patcher.assert_called_with('ignoring invalid read', 'BAD_READ', level=logging.DEBUG)

    def test_reference_id(self):
        fh = MockBamFileHandle({'1': 0})
        b = BamCache(fh)
        self.assertEqual(0, b.reference_id('1'))
        with self.assertRaises(KeyError):
            b.reference_id('2')

    def test_get_read_reference_name(self):
        fh = MockBamFileHandle({'1': 0})
        b = BamCache(fh)
        r = MockRead('name', 0)
        self.assertEqual('1', b.get_read_reference_name(r))

    def test_generate_fetch_bins_single(self):
        self.assertEqual([(1, 100)], BamCache._generate_fetch_bins(1, 100, 1, 1))

    def test_generate_fetch_bins_multi(self):
        self.assertEqual([(1, 50), (51, 100)], BamCache._generate_fetch_bins(1, 100, 2, 1))
        self.assertEqual(
            [(1, 20), (21, 40), (41, 60), (61, 80), (81, 100)], BamCache._generate_fetch_bins(1, 100, 5, 1))

    def test_generate_fetch_bins_large_min_size(self):
        self.assertEqual([(1, 50), (51, 100)], BamCache._generate_fetch_bins(1, 100, 5, 50))

    def test_fetch_single_read(self):
        b = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
        s = b.fetch_from_bins('reference3', 1382, 1383, read_limit=1, sample_bins=1)
        self.assertEqual(1, len(s))
        r = list(s)[0]
        self.assertEqual('HISEQX1_11:4:2122:14275:37717:split', r.qname)
        b.close()

    def test_get_mate(self):
        # dependant on fetch working
        b = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
        s = b.fetch_from_bins('reference3', 1382, 1383, read_limit=1, sample_bins=1)
        self.assertEqual(1, len(s))
        r = list(s)[0]
        self.assertEqual('HISEQX1_11:4:2122:14275:37717:split', r.qname)
        o = b.get_mate(r, allow_file_access=True)
        self.assertEqual(1, len(o))
        self.assertEqual('HISEQX1_11:4:2122:14275:37717:split', o[0].qname)


class TestModule(unittest.TestCase):
    """
    test class for functions in the validate namespace
    that are not associated with a class
    """

    def test_alphabet_matching(self):
        self.assertTrue(DNA_ALPHABET.match('N', 'A'))
        self.assertTrue(DNA_ALPHABET.match('A', 'N'))

    def test_breakpoint_pos(self):
        # ==========+++++++++>
        r = MockRead(reference_start=10, cigar=[(CIGAR.M, 10), (CIGAR.S, 10)])
        self.assertEqual(19, _read.breakpoint_pos(r))

        with self.assertRaises(AttributeError):
            breakpoint_pos(r, ORIENT.RIGHT)

        self.assertEqual(19, _read.breakpoint_pos(r, ORIENT.LEFT))

        # ++++++++++=========>
        r = MockRead(reference_start=10, cigar=[(CIGAR.S, 10), (CIGAR.M, 10)])
        self.assertEqual(10, _read.breakpoint_pos(r))

        with self.assertRaises(AttributeError):
            breakpoint_pos(r, ORIENT.LEFT)

        self.assertEqual(10, _read.breakpoint_pos(r, ORIENT.RIGHT))

        with self.assertRaises(AttributeError):
            r = MockRead(reference_start=10, cigar=[(CIGAR.X, 10), (CIGAR.M, 10)])
            _read.breakpoint_pos(r, ORIENT.LEFT)

    def test_nsb_align(self):
        ref = 'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG' \
              'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        seq = 'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC' \
              'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        _read.nsb_align(ref, seq)
        # GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG


class TestNsbAlign(unittest.TestCase):

    def test_length_seq_le_ref(self):
        ref = 'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG' \
              'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        seq = 'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC' \
              'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        alignment = _read.nsb_align(ref, seq)
        self.assertEqual(1, len(alignment))
        alignment = _read.nsb_align(ref, seq, min_consecutive_match=20)
        self.assertEqual(0, len(alignment))

    def test_length_ref_le_seq(self):
        pass

    def test_length_ref_eq_seq(self):
        pass

    @timeout_decorator.timeout(5)
    def test_long_ref_seq(self):
        ref = str(REFERENCE_GENOME['test_bam_long_ref'].seq)
        seq = 'TGAGGTCAGGAGTTTGAGACCAGCCTGGACAACATGGTGAAACCCCATCTCTACTAAAAATACAAAAAAATTAGCCAGGCATGGTGGTGGATGCCTGTAAT' \
            'CGCAGCTACTCAGGAGATCGGAAG'
        alignment = _read.nsb_align(ref, seq, min_consecutive_match=6)
        self.assertEqual(1, len(alignment))

    def test_left_softclipping(self):
        ref = 'TAAGCTTCTTCCTTTTTCTATGCCACCTACATAGGCATTTTGCATGGTCAGATTGGAATTTACATAATGCATACATGCAAAGAAATATATAGAAGCCAGATATATAAGGTAGTACATTGGCAGGCTTCATATATATAGACTCCCCCATATTGTCTATATGCTAAAAAAGTATTTTAAATCCTTAAATTTTATTTTTGTTCTCTGCATTTGAAATCTTTATCAACTAGGTCATGAAAATAGCCAGTCGGTTCTCCTTTTGGTCTATTAGAATAAAATCTGGACTGCAACTGAGAAGCAGAAGGTAATGTCAGAATGTAT'
        seq = 'GCTAAAAAAGTATTTTAAATCCTTAAATGTTATTTTTGTTCTC'
        alignments = _read.nsb_align(ref, seq, min_consecutive_match=6)
        self.assertEqual(1, len(alignments))
        print(alignments)
        seq = 'CTTATAAAGCTGGAGTATCTGCTGAGAGCATCAGGAATTGACATCTAGGATAATGAGAGAAGGCTGATCATGGACAACATATAGCCTTTCTAGTAGATGCAGCTGAGGCTAAAAAAGTATTTTAAATCCTTAAATGTTATTTTTGTTCTC'
        alignments = _read.nsb_align(ref, seq, min_consecutive_match=6, min_overlap_percent=0.5)
        self.assertEqual(1, len(alignments))

    def test_min_overlap(self):
        ref = 'ATTACATTAAAGATTCAAACTCCTAGAGTTTTTTTGATTTTTAGTATGATCTTTAGATAAAAAAAAAGGAAGAAAAAGAAAAAAAAACAGAGTCTATTAAGGCATCTTCTATGGTCAGATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGAT'
        seq = 'GATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTGTATATATATATATATATATATATATATATATACCCAGTTTCAAGCAG'
        alignments = _read.nsb_align(
            ref, seq,
            min_consecutive_match=15,
            min_match=0.95,
            min_overlap_percent=(len(seq) - 15) / len(seq)
        )
        print(alignments)
        self.assertEqual(0, len(alignments))


class TestReadPairStrand(unittest.TestCase):
    def setUp(self):
        self.read1_pos_neg = MockRead(is_reverse=False, is_read1=True, mate_is_reverse=True)
        assert(not self.read1_pos_neg.is_read2)
        self.read1_neg_pos = MockRead(is_reverse=True, is_read1=True, mate_is_reverse=False)
        self.read1_pos_pos = MockRead(is_reverse=False, is_read1=True, mate_is_reverse=False)
        self.read1_neg_neg = MockRead(is_reverse=True, is_read1=True, mate_is_reverse=True)

        self.read2_pos_neg = MockRead(is_reverse=True, is_read1=False, mate_is_reverse=True)
        assert(self.read2_pos_neg.is_read2)
        self.read2_neg_pos = MockRead(is_reverse=False, is_read1=False, mate_is_reverse=False)
        self.read2_pos_pos = MockRead(is_reverse=False, is_read1=False, mate_is_reverse=False)
        self.read2_neg_neg = MockRead(is_reverse=True, is_read1=False, mate_is_reverse=True)

        self.unpaired_pos = MockRead(is_reverse=False, is_paired=False)
        self.unpaired_neg = MockRead(is_reverse=True, is_paired=False)

    def test_read_pair_strand_det1_read1(self):
        self.assertEqual(STRAND.POS, sequenced_strand(self.read1_pos_neg, strand_determining_read=1))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read1_neg_pos, strand_determining_read=1))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read1_pos_pos, strand_determining_read=1))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read1_neg_neg, strand_determining_read=1))

    def test_read_pair_strand_det1_read2(self):
        self.assertEqual(STRAND.POS, sequenced_strand(self.read2_pos_neg, strand_determining_read=1))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read2_neg_pos, strand_determining_read=1))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read2_pos_pos, strand_determining_read=1))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read2_neg_neg, strand_determining_read=1))

    def test_read_pair_strand_det2_read2(self):
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read2_pos_neg, strand_determining_read=2))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read2_neg_pos, strand_determining_read=2))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read2_pos_pos, strand_determining_read=2))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read2_neg_neg, strand_determining_read=2))

    def test_read_pair_strand_det2_read1(self):
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read1_pos_neg, strand_determining_read=2))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read1_neg_pos, strand_determining_read=2))
        self.assertEqual(STRAND.NEG, sequenced_strand(self.read1_pos_pos, strand_determining_read=2))
        self.assertEqual(STRAND.POS, sequenced_strand(self.read1_neg_neg, strand_determining_read=2))

    def test_read_pair_strand_unpaired(self):
        with self.assertRaises(ValueError):
            sequenced_strand(self.unpaired_pos)
        with self.assertRaises(ValueError):
            sequenced_strand(self.unpaired_neg)

    def test_read_pair_strand_det_error(self):
        with self.assertRaises(ValueError):
            sequenced_strand(self.read1_pos_neg, strand_determining_read=3)


class TestReadPairType(unittest.TestCase):
    def setUp(self):
        self.LR = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=True
        )
        self.LL = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=False
        )
        self.RR = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=True
        )
        self.RL = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=False
        )

    def test_read_pair_type_LR(self):
        self.assertEqual(READ_PAIR_TYPE.LR, read_pair_type(self.LR))

    def test_read_pair_type_LL(self):
        self.assertEqual(READ_PAIR_TYPE.LL, read_pair_type(self.LL))

    def test_read_pair_type_RR(self):
        self.assertEqual(READ_PAIR_TYPE.RR, read_pair_type(self.RR))

    def test_read_pair_type_RL(self):
        self.assertEqual(READ_PAIR_TYPE.RL, read_pair_type(self.RL))

    def test_orientation_supports_type_deletion(self):
        self.assertTrue(orientation_supports_type(self.LR, SVTYPE.DEL))
        self.assertFalse(orientation_supports_type(self.RL, SVTYPE.DEL))
        self.assertFalse(orientation_supports_type(self.LL, SVTYPE.DEL))
        self.assertFalse(orientation_supports_type(self.RR, SVTYPE.DEL))

    def test_orientation_supports_type_insertion(self):
        self.assertTrue(orientation_supports_type(self.LR, SVTYPE.INS))
        self.assertFalse(orientation_supports_type(self.RL, SVTYPE.INS))
        self.assertFalse(orientation_supports_type(self.LL, SVTYPE.INS))
        self.assertFalse(orientation_supports_type(self.RR, SVTYPE.INS))

    def test_orientation_supports_type_inversion(self):
        self.assertFalse(orientation_supports_type(self.LR, SVTYPE.INV))
        self.assertFalse(orientation_supports_type(self.RL, SVTYPE.INV))
        self.assertTrue(orientation_supports_type(self.LL, SVTYPE.INV))
        self.assertTrue(orientation_supports_type(self.RR, SVTYPE.INV))

    def test_orientation_supports_type_translocation_inversion(self):
        self.assertFalse(orientation_supports_type(self.LR, SVTYPE.ITRANS))
        self.assertFalse(orientation_supports_type(self.RL, SVTYPE.ITRANS))
        self.assertTrue(orientation_supports_type(self.LL, SVTYPE.ITRANS))
        self.assertTrue(orientation_supports_type(self.RR, SVTYPE.ITRANS))

    def test_orientation_supports_type_trans_duplication(self):
        self.assertFalse(orientation_supports_type(self.LR, SVTYPE.DUP))
        self.assertTrue(orientation_supports_type(self.RL, SVTYPE.DUP))
        self.assertFalse(orientation_supports_type(self.LL, SVTYPE.DUP))
        self.assertFalse(orientation_supports_type(self.RR, SVTYPE.DUP))

    def test_orientation_supports_type_translocation(self):
        self.assertTrue(orientation_supports_type(self.LR, SVTYPE.TRANS))
        self.assertTrue(orientation_supports_type(self.RL, SVTYPE.TRANS))
        self.assertFalse(orientation_supports_type(self.LL, SVTYPE.TRANS))
        self.assertFalse(orientation_supports_type(self.RR, SVTYPE.TRANS))


class TestHistogram(unittest.TestCase):
    def test_add(self):
        h = Histogram()
        h.add(1)
        h.add(1)
        self.assertEqual(2, h[1])
        h.add(1, 4)
        self.assertEqual(6, h[1])

    def test_median(self):
        h = Histogram()
        for i in range(1, 11):
            h.add(i)
        self.assertEqual(5.5, h.median())
        h.add(11)
        self.assertEqual(6, h.median())

    def test_distib_stderr(self):
        h = Histogram()
        for i in range(0, 11):
            h.add(i)
        for i in range(4, 8):
            h.add(i)
        m = h.median()
        self.assertEqual(5, m)
        err = h.distribution_stderr(m, 1)
        self.assertEqual(116 / 15, err)

    def test_add_operator(self):
        x = Histogram()
        y = Histogram()
        x.add(1)
        y.add(1, 4)
        z = x + y
        self.assertEqual(1, x[1])
        self.assertEqual(4, y[1])
        self.assertEqual(5, z[1])


class TestBamStats(unittest.TestCase):
    def test_genome_bam_stats(self):
        bamfh = BamCache(get_data('mock_reads_for_events.sorted.bam'))
        stats = compute_genome_bam_stats(
            bamfh,
            1000,
            100,
            min_mapping_quality=1,
            sample_cap=10000,
            distribution_fraction=0.99
        )
        self.assertGreaterEqual(50, abs(stats.median_fragment_size - 420))
        self.assertEqual(150, stats.read_length)
        bamfh.close()

    def test_trans_bam_stats(self):
        bamfh = BamCache(get_data('mock_trans_reads_for_events.sorted.bam'))
        annotations = load_reference_genes(get_data('mock_annotations.json'))
        stats = compute_transcriptome_bam_stats(
            bamfh,
            annotations,
            100,
            min_mapping_quality=1,
            stranded=True,
            sample_cap=10000,
            distribution_fraction=0.99
        )
        self.assertTrue(abs(stats.median_fragment_size - 185) < 5)
        self.assertEqual(75, stats.read_length)
        self.assertTrue(stats.stdev_fragment_size < 50)
        bamfh.close()


class TestMapRefRangeToQueryRange(unittest.TestCase):
    def setUp(self):
        self.contig_read = MockRead(
            cigar=_cigar.convert_string_to_cigar('275M18I12041D278M'),
            reference_start=89700025,
            reference_name='10'
        )

    def test_full_aligned_portion(self):
        ref_range = Interval(89700026, 89712619)
        qrange = _read.map_ref_range_to_query_range(self.contig_read, ref_range)
        self.assertEqual(571, len(qrange))
        self.assertEqual(1, qrange.start)
        self.assertEqual(571, qrange.end)

    def test_multiple_events(self):
        ref_range = Interval(89700067, 89712347)
        qrange = _read.map_ref_range_to_query_range(self.contig_read, ref_range)
        self.assertEqual(len(ref_range) - 12041 + 18, len(qrange))

    def test_no_events(self):
        ref_range = Interval(89700031, 89700040)
        qrange = _read.map_ref_range_to_query_range(self.contig_read, ref_range)
        self.assertEqual(10, len(qrange))
        self.assertEqual(6, qrange.start)
        self.assertEqual(15, qrange.end)
