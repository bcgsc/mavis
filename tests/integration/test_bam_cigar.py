from mavis.constants import CIGAR
from mavis.bam.cigar import recompute_cigar_mismatch, hgvs_standardize_cigar, extend_softclipping
from mavis.bam.cigar import alignment_matches, longest_fuzzy_match, join, score, match_percent, compute
from mavis.bam.cigar import smallest_nonoverlapping_repeat, convert_for_igv
from mavis.annotate import load_reference_genome
import unittest
import warnings
from . import REFERENCE_GENOME_FILE
from . import MockRead


REFERENCE_GENOME = None


def setUpModule():
    warnings.simplefilter("ignore")
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestRecomputeCigarMismatch(unittest.TestCase):

    def test_simple(self):
        r = MockRead(
            reference_start=1456,
            query_sequence='CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.M, 45)]
        )
        self.assertEqual([(CIGAR.EQ, 45)], recompute_cigar_mismatch(r, REFERENCE_GENOME['fake']))

    def test_hardclipping(self):
        r = MockRead(
            reference_start=1456,
            query_sequence='CCCAAACAAC'
                           'TATAAATTTT'
                           'GTAATACCTA'
                           'GAACAATATA'
                           'AATAT',
            cigar=[(CIGAR.H, 20), (CIGAR.M, 45)]
        )
        self.assertEqual([(CIGAR.H, 20), (CIGAR.EQ, 45)], recompute_cigar_mismatch(r, REFERENCE_GENOME['fake']))

    def test_with_events(self):
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
            recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )

    def test_mismatch_to_mismatch(self):
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
            recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )

    def test_m_to_mismatch(self):
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
            recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])
        )

    def test_invalid_cigar_value(self):
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
            recompute_cigar_mismatch(r, REFERENCE_GENOME['fake'])


class TestExtendSoftclipping(unittest.TestCase):

    def test_simple(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.M, 10)], 0),
            extend_softclipping([(CIGAR.S, 10), (CIGAR.M, 10)], 1)
        )

    def test_deletions(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.M, 10)], 1),
            extend_softclipping([(CIGAR.I, 10), (CIGAR.D, 1), (CIGAR.M, 10)], 1)
        )

    def test_mismatch(self):
        with self.assertRaises(AttributeError):
            extend_softclipping([(CIGAR.X, 10), (CIGAR.M, 20), (CIGAR.X, 10)], 30)

    def test_insert(self):
        self.assertEqual(
            ([(CIGAR.S, 10), (CIGAR.S, 2), (CIGAR.S, 5), (CIGAR.M, 10), (CIGAR.S, 5)], 2),
            extend_softclipping([(CIGAR.S, 10), (CIGAR.M, 2), (CIGAR.I, 5), (CIGAR.M, 10), (CIGAR.I, 5)], 5)
        )

    def test_hardclipping(self):
        c = [(CIGAR.H, 10), (CIGAR.EQ, 10)]
        cnew, prefix = extend_softclipping(c, 1)
        self.assertEqual(0, prefix)
        self.assertEqual(c, cnew)


class TestCigarTools(unittest.TestCase):

    def test_alignment_matches(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        self.assertEqual(30, alignment_matches(c))

    def test_join(self):
        c = [(CIGAR.M, 10), (CIGAR.X, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 10), (CIGAR.X, 20)], join(c))
        k = [(CIGAR.X, 10), (CIGAR.M, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 10), (CIGAR.X, 30), (CIGAR.M, 10), (CIGAR.X, 10)], join(c, k))
        k = [(4, 1), (4, 2), (7, 5), (8, 7), (7, 2), (8, 5), (7, 28), (8, 1), (7, 99)]
        self.assertEqual([(4, 3), (7, 5), (8, 7), (7, 2), (8, 5), (7, 28), (8, 1), (7, 99)], join(k))

    def test_join_hardclipping(self):
        c = [(CIGAR.H, 10), (CIGAR.M, 10), (CIGAR.X, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.H, 10), (CIGAR.M, 10), (CIGAR.X, 20)], join(c))

    def test_longest_fuzzy_match(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(15, longest_fuzzy_match(c, 1))
        self.assertEqual(10, longest_fuzzy_match(c, 0))
        self.assertEqual(16, longest_fuzzy_match(c, 2))
        self.assertEqual(16, longest_fuzzy_match(c, 4))

    def test_score(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(22, score(c))

    def test_score_error(self):
        with self.assertRaises(AssertionError):
            c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (99, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
            score(c)

    def test_match_percent(self):
        c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        self.assertEqual(0.8, match_percent(c))
        with self.assertRaises(AttributeError):
            match_percent([(CIGAR.M, 100)])
        with self.assertRaises(AttributeError):
            match_percent([(CIGAR.S, 100)])

    def test_compute(self):
        # GTGAGTAAATTCAACATCGTTTTT
        # aacttagAATTCAAC---------
        self.assertEqual(
            ([(CIGAR.S, 7), (CIGAR.EQ, 8)], 7),
            compute('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
        )
        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 7),
            compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------')
        )
        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 7),
            compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 5), (CIGAR.I, 2), (CIGAR.EQ, 1)], 7),
            compute('GTGAGTAAATTC--CATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

        with self.assertRaises(AttributeError):
            compute('CCTG', 'CCG')

        self.assertEqual(
            ([(CIGAR.EQ, 2), (CIGAR.X, 2)], 0),
            compute('CCTG', 'CCGT', min_exact_to_stop_softclipping=10)
        )

        self.assertEqual(
            ([(CIGAR.S, 5), (CIGAR.EQ, 8)], 5),
            compute('--GAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------', False)
        )

    def test_convert_for_igv(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        self.assertEqual([(CIGAR.M, 30)], convert_for_igv(c))


class TestHgvsStandardizeCigars(unittest.TestCase):
    def no_change_aligned(self):
        ref = 'AAATTTGGGCCCAATT'
        read = MockRead('name', '1', 1, cigar=[(CIGAR.M, 10)], query_sequence='AAATTTGGGC')
        self.assertEqual([(CIGAR.M, 10)], hgvs_standardize_cigar(read, ref))

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
            [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)], hgvs_standardize_cigar(read, ref))

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
            [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)], hgvs_standardize_cigar(read, ref))

    def test_insertion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.I, 4), (CIGAR.EQ, 8)]
        )
        self.assertEqual([(CIGAR.EQ, 8), (CIGAR.I, 4), (CIGAR.EQ, 4)], hgvs_standardize_cigar(read, ref))

    def test_deletion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 4), (CIGAR.EQ, 12)]
        )
        self.assertEqual([(CIGAR.EQ, 12), (CIGAR.D, 4), (CIGAR.EQ, 4)], hgvs_standardize_cigar(read, ref))

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
            [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 4), (CIGAR.EQ, 12)], hgvs_standardize_cigar(read, ref))

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
            [(CIGAR.EQ, 4), (CIGAR.I, 6), (CIGAR.D, 5), (CIGAR.EQ, 12)], hgvs_standardize_cigar(read, ref))

    def test_bubble_sort_indel_sections_drop_mismatch_with_hardclipping(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGAC' 'TTTTT' 'ACGA' 'GATC',
            cigar=[(CIGAR.H, 10), (CIGAR.EQ, 4), (CIGAR.X, 1), (CIGAR.D, 2), (CIGAR.I, 3), (CIGAR.D, 2), (CIGAR.I, 2), (CIGAR.EQ, 12)]
        )
        self.assertEqual(
            [(CIGAR.H, 10), (CIGAR.EQ, 4), (CIGAR.I, 6), (CIGAR.D, 5), (CIGAR.EQ, 12)], hgvs_standardize_cigar(read, ref))

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
            [(CIGAR.EQ, 9 + 12), (CIGAR.D, 2), (CIGAR.EQ, 7)], hgvs_standardize_cigar(read, ref))
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
            [(CIGAR.EQ, 61 + 15), (CIGAR.I, 2), (CIGAR.EQ, 87 - 15)], hgvs_standardize_cigar(read, ref))

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
            hgvs_standardize_cigar(read, ref))

    def test_smallest_nonoverlapping_repeat(self):
        s = 'ATATATATAA'
        self.assertEqual(s, smallest_nonoverlapping_repeat(s))
        s = 'ATATATATAT'
        self.assertEqual('AT', smallest_nonoverlapping_repeat(s))
        s = 'AAAAAA'
        self.assertEqual('A', smallest_nonoverlapping_repeat(s))
        s = 'ATGGCATGGCATGGC'
        self.assertEqual('ATGGC', smallest_nonoverlapping_repeat(s))
