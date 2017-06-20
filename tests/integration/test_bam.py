from mavis.constants import *
from mavis.bam import read as read_tools
from mavis.bam import cigar as cigar_tools
from mavis.bam.read import sequenced_strand, read_pair_type, breakpoint_pos, orientation_supports_type, get_samtools_version
from mavis.bam.cache import BamCache
from mavis.bam.stats import Histogram, compute_transcriptome_bam_stats, compute_genome_bam_stats
from mavis.annotate import load_reference_genome, load_reference_genes
import pysam
import unittest
import warnings
import os
from . import MockRead, MockBamFileHandle
from . import REFERENCE_GENOME_FILE, TRANSCRIPTOME_BAM_INPUT, FULL_REFERENCE_ANNOTATIONS_FILE_JSON
from . import BAM_INPUT, FULL_BAM_INPUT
from .config import samtools_versions
import timeout_decorator


REFERENCE_GENOME = None


def setUpModule():
    warnings.simplefilter("ignore")
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestGetSamtoolsVersion(unittest.TestCase):
    
    def test_get_samtools_version(self):
        env = os.environ
        for version, path in samtools_versions.items():
            env['PATH'] = os.path.dirname(path) + ':' + env['PATH']
            self.assertEqual(version, get_samtools_version())


class TestBamCache(unittest.TestCase):

    def test___init__(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        self.assertEqual(fh, b.fh)

    def test_add_read(self):
        fh = MockBamFileHandle()
        b = BamCache(fh)
        r = MockRead('name')
        b.add_read(r)
        self.assertEqual(1, len(b.cache.values()))
        self.assertEqual(set([r]), b.cache['name'])

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
        b = BamCache(BAM_INPUT)
        s = b.fetch_from_bins('reference3', 1382, 1383, read_limit=1, sample_bins=1)
        self.assertEqual(1, len(s))
        r = list(s)[0]
        self.assertEqual('HISEQX1_11:4:2122:14275:37717:split', r.qname)
        b.close()

    def test_get_mate(self):
        #dependant on fetch working
        b = BamCache(BAM_INPUT)
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
        self.assertEqual(19, read_tools.breakpoint_pos(r))

        with self.assertRaises(AttributeError):
            breakpoint_pos(r, ORIENT.RIGHT)

        self.assertEqual(19, read_tools.breakpoint_pos(r, ORIENT.LEFT))

        # ++++++++++=========>
        r = MockRead(reference_start=10, cigar=[(CIGAR.S, 10), (CIGAR.M, 10)])
        self.assertEqual(10, read_tools.breakpoint_pos(r))

        with self.assertRaises(AttributeError):
            breakpoint_pos(r, ORIENT.LEFT)

        self.assertEqual(10, read_tools.breakpoint_pos(r, ORIENT.RIGHT))

        with self.assertRaises(AttributeError):
            r = MockRead(reference_start=10, cigar=[(CIGAR.X, 10), (CIGAR.M, 10)])
            read_tools.breakpoint_pos(r, ORIENT.LEFT)

    def test_nsb_align(self):
        ref = 'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG' \
              'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        seq = 'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC' \
              'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        alignment = read_tools.nsb_align(ref, seq)
        # GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG


class TestNsbAlign(unittest.TestCase):

    def test_length_seq_le_ref(self):
        ref = 'GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAG' \
              'TTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG'
        seq = 'TGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGAC' \
              'TGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAAC'
        alignment = read_tools.nsb_align(ref, seq)
        self.assertEqual(1, len(alignment))
        alignment = read_tools.nsb_align(ref, seq, min_consecutive_match=20)
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
        alignment = read_tools.nsb_align(ref, seq, min_consecutive_match=6)
        self.assertEqual(1, len(alignment))



class TestCigarTools(unittest.TestCase):

    def test_recompute_cigar_mismatch(self):
        r = MockRead(
            reference_start=1456,
            query_sequence='CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.M, 45)]
        )
        self.assertEqual([(CIGAR.EQ, 45)], cigar_tools.recompute_cigar_mismatch(r, REFERENCE_GENOME['fake']))

        r = MockRead(
            reference_start=1456,
            query_sequence='TATA'
                           'CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.S, 4), (CIGAR.M, 10), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.M, 25)]
        )
        self.assertEqual(
            [(CIGAR.S, 4), (CIGAR.EQ, 10), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.EQ, 25)],
            cigar_tools.recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )
        r = MockRead(
            reference_start=1452,
            query_sequence='CAGC'
                           'CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.X, 4), (CIGAR.M, 10), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.M, 25)]
        )
        self.assertEqual(
            [(CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.EQ, 25)],
            cigar_tools.recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )
        r = MockRead(
            reference_start=1452,
            query_sequence='CAGC'
                           'CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.M, 14), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.M, 25)]
        )
        self.assertEqual(
            [(CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.EQ, 25)],
            cigar_tools.recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )

    def test_recompute_cigar_mismatch_invalid_cigar_value(self):
        r = MockRead(
            reference_start=1452,
            query_sequence='CAGC'
                           'CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(50, 14), (CIGAR.D, 10), (CIGAR.I, 10), (CIGAR.M, 25)]
        )
        with self.assertRaises(NotImplementedError):
            cigar_tools.recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])

    def test_longest_fuzzy_match(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(15, cigar_tools.longest_fuzzy_match(c, 1))
        self.assertEqual(10, cigar_tools.longest_fuzzy_match(c, 0))
        self.assertEqual(16, cigar_tools.longest_fuzzy_match(c, 2))
        self.assertEqual(16, cigar_tools.longest_fuzzy_match(c, 4))

    def test_score(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(22, cigar_tools.score(c))

    def test_score_error(self):
        with self.assertRaises(AssertionError):
            c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (99, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
            cigar_tools.score(c)

    def test_match_percent(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(0.8, cigar_tools.match_percent(c))
        with self.assertRaises(AttributeError):
            cigar_tools.match_percent([(CIGAR.M, 100)])
        with self.assertRaises(AttributeError):
            cigar_tools.match_percent([(CIGAR.S, 100)])

    def test_compute(self):
        # GTGAGTAAATTCAACATCGTTTTT
        # aacttagAATTCAAC---------
        self.assertEqual(
            ([(CIGAR.S, 7), (CIGAR.EQ, 8)], 7),
            cigar_tools.compute('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
        )
        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 7),
            cigar_tools.compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------')
        )
        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 7),
            cigar_tools.compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 5), (CIGAR.I, 2), (CIGAR.EQ, 1)], 7),
            cigar_tools.compute('GTGAGTAAATTC--CATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

        with self.assertRaises(AttributeError):
            cigar_tools.compute('CCTG', 'CCG')

        self.assertEqual(
            ([(CIGAR.EQ, 2), (CIGAR.X, 2)], 0),
            cigar_tools.compute('CCTG', 'CCGT', min_exact_to_stop_softclipping=10)
        )

        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 5),
            cigar_tools.compute('--GAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

    def test_convert_for_igv(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 30)], cigar_tools.convert_for_igv(c))

    def test_extend_softclipping(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.M, 10)], 0),
            cigar_tools.extend_softclipping([(CIGAR.S, 10), (CIGAR.M, 10)], 1)
        )

    def test_extend_softclipping_deletions(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.M, 10)], 1),
            cigar_tools.extend_softclipping([(CIGAR.I, 10), (CIGAR.D, 1), (CIGAR.M, 10)], 1)
        )

    def test_extend_softclipping_mismatch(self):
        with self.assertRaises(AttributeError):
            cigar_tools.extend_softclipping([(CIGAR.X, 10), (CIGAR.M, 20), (CIGAR.X, 10)], 30)

    def test_extend_softclipping_insert(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.S, 2), (CIGAR.S, 5), (CIGAR.M, 10), (CIGAR.S, 5)], 2),
            cigar_tools.extend_softclipping([(CIGAR.S, 10), (CIGAR.M, 2), (CIGAR.I, 5), (CIGAR.M, 10), (CIGAR.I, 5)], 5)
        )

    def test_alignment_matches(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        self.assertEqual(30, cigar_tools.alignment_matches(c))

    def test_join(self):
        c = [(CIGAR.M, 10), (CIGAR.X, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 10), (CIGAR.X, 20)], cigar_tools.join(c))
        k = [(CIGAR.X, 10), (CIGAR.M, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 10), (CIGAR.X, 30), (CIGAR.M, 10), (CIGAR.X, 10)], cigar_tools.join(c, k))
        k = [(4, 1), (4, 2), (7, 5), (8, 7), (7, 2), (8, 5), (7, 28), (8, 1), (7, 99)]
        self.assertEqual([(4, 3), (7, 5), (8, 7), (7, 2), (8, 5), (7, 28), (8, 1), (7, 99)], cigar_tools.join(k))


class TestHgvsStandardizeCigars(unittest.TestCase):
    def no_change_aligned(self):
        ref = 'AAATTTGGGCCCAATT'
        read = MockRead('name', '1', 1, cigar=[(CIGAR.M, 10)], query_sequence='AAATTTGGGC')
        self.assertEqual([(CIGAR.M, 10)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def no_change_proper_indel(self):
        ref = 'ATAGGC' 'ATCTACGAG' 'ATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCTAC' 'CCC' 'ATCG',
            cigar=[(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def ins_after_deletion(self):
        ref = 'ATAGGC' 'ATCTACGAG' 'ATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCTAC' 'CCC' 'ATCG',
            cigar=[(CIGAR.EQ, 6), (CIGAR.D, 3), (CIGAR.I, 3), (CIGAR.EQ, 4)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def test_insertion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.I, 4), (CIGAR.EQ, 8)]
        )
        self.assertEqual([(CIGAR.EQ, 8), (CIGAR.I, 4), (CIGAR.EQ, 4)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def test_deletion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 4), (CIGAR.EQ, 12)]
        )
        self.assertEqual([(CIGAR.EQ, 12), (CIGAR.D, 4), (CIGAR.EQ, 4)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def test_bubble_sort_indel_sections(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'TTTTT' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 2), (CIGAR.I, 3), (CIGAR.D, 2), (CIGAR.I, 2), (CIGAR.EQ, 12)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 4), (CIGAR.EQ, 12)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def test_bubble_sort_indel_sections_drop_mismatch(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGAC' 'TTTTT' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.X, 1), (CIGAR.D, 2), (CIGAR.I, 3), (CIGAR.D, 2), (CIGAR.I, 2), (CIGAR.EQ, 12)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 4), (CIGAR.I, 6), (CIGAR.D, 5), (CIGAR.EQ, 12)], cigar_tools.hgvs_standardize_cigar(read, ref))

    def test_homopolymer_even_odd(self):
        ref = 'ATCGAGAT' + 'A' * 15 + 'TCGAGAT'
        read = MockRead(
            'name',
            1,
            1,
            query_sequence='ATCGAGATA' + 'A' * 12 + 'TCGAGAT',
            cigar=[(CIGAR.EQ, 8), (CIGAR.D, 2), (CIGAR.EQ, 20)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 9 + 12), (CIGAR.D, 2), (CIGAR.EQ, 7)], cigar_tools.hgvs_standardize_cigar(read, ref))
        ref = 'CCCCGGCTCATGTCTGGTTTTGTTTTCCGGGGGCGGGGGGGCTCCCTGGGGATGATGGTGATTTTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAA' \
              'TGAGGCAGAGACAATGTGGGGAGCGAGAGAGGGGAAAAGGACGGGGGAGG'

        read = MockRead(
            'name', '1', 0, 149,
            query_sequence=(
                'CCCCGGCTCATGTCTGGTTTTGTTTTCCGGGGGCGGGGGGGCTCCCTGGGGATGATGGTGATTTTTTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAA'
                'TGAGGCAGAGACAATGTGGGGAGCGAGAGAGGGGAAAAGGACGGGGGAGG'),
            cigar=[(CIGAR.EQ, 61), (CIGAR.I, 2), (CIGAR.EQ, 87)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 61 + 15), (CIGAR.I, 2), (CIGAR.EQ, 87 - 15)], cigar_tools.hgvs_standardize_cigar(read, ref))

        ref = 'CCTCCTCGGTCGGGCAGATCTTTCAGAAGCAGGAGCCCAGGATCATGTCTGGTTTTGTTTTCCGAGGGCGAGGGGGCTCCCTGAGGATGATGGTGATTT' \
            'TTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAATGAGGCAGAGACA'

        read = MockRead(
            'name', '1', 0, 149,
            query_sequence=(
                'CCCCTCCTCGGTCGGGCAGATCTTTCAGAAGCAGGAGCCCAGGATCATGTCTGGTTTTGTTTTCCGAGGGCGAGGGGGCTCCCTGAGGATGATGGTGATTTT'
                'TTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAATGAGGCAGAGACA'),
            cigar=[(CIGAR.S, 2), (CIGAR.EQ, 96), (CIGAR.I, 2), (CIGAR.EQ, 50)]
        )
        self.assertEqual(
            [(CIGAR.S, 2), (CIGAR.EQ, 96 + 15), (CIGAR.I, 2), (CIGAR.EQ, 50 - 15)],
            cigar_tools.hgvs_standardize_cigar(read, ref))


    def test_smallest_nonoverlapping_repeat(self):
        s = 'ATATATATAA'
        self.assertEqual(s, cigar_tools.smallest_nonoverlapping_repeat(s))
        s = 'ATATATATAT'
        self.assertEqual('AT', cigar_tools.smallest_nonoverlapping_repeat(s))
        s = 'AAAAAA'
        self.assertEqual('A', cigar_tools.smallest_nonoverlapping_repeat(s))
        s = 'ATGGCATGGCATGGC'
        self.assertEqual('ATGGC', cigar_tools.smallest_nonoverlapping_repeat(s))


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
        bamfh = None
        try:
            bamfh = pysam.AlignmentFile(FULL_BAM_INPUT, 'rb')
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
        finally:
            try:
                bamfh.close()
            except AttributeError:
                pass

    def test_trans_bam_stats(self):
        bamfh = None
        try:
            bamfh = pysam.AlignmentFile(TRANSCRIPTOME_BAM_INPUT, 'rb')
            annotations = load_reference_genes(FULL_REFERENCE_ANNOTATIONS_FILE_JSON)
            stats = compute_transcriptome_bam_stats(
                bamfh,
                annotations,
                100,
                best_transcripts_only=False,
                min_mapping_quality=1,
                stranded=True,
                sample_cap=10000,
                distribution_fraction=0.99
            )
            self.assertTrue(abs(stats.median_fragment_size - 185) < 5)
            self.assertEqual(75, stats.read_length)
            self.assertTrue(stats.stdev_fragment_size < 50)
        finally:
            try:
                bamfh.close()
            except AttributeError:
                pass


