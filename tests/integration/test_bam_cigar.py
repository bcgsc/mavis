import unittest
import warnings

from mavis.annotate.file_io import load_reference_genome
from mavis.bam.cigar import alignment_matches, compute, convert_for_igv, convert_string_to_cigar, extend_softclipping, hgvs_standardize_cigar, join, longest_fuzzy_match, match_percent, merge_internal_events, recompute_cigar_mismatch, score
from mavis.constants import CIGAR
from mavis.bam.read import SamRead
from mavis.bam import read as _read
import timeout_decorator

from . import MockRead, MockObject
from ..util import get_data


REFERENCE_GENOME = None


def setUpModule():
    warnings.simplefilter('ignore')
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
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


class TestExtendSoftclipping(unittest.TestCase):

    def test_softclipped_right(self):
        c = convert_string_to_cigar('70=2X1=8X4=1X1=4X1=6X1=4X1=4X2=5X3=3X1=4X1=3X1=14X1=1X2=1S')
        cnew, prefix = extend_softclipping(c, 6)
        self.assertEqual(0, prefix)
        self.assertEqual(convert_string_to_cigar('70=80S'), cnew)


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
        rseq = 'ATAGGC' 'ATCT'        'GG'    'GA' 'GCGA' 'GATCGCTACG'
        qseq = 'ATCT'  'TTT'     'TT'     'GCGA' 'GATC'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence=qseq,
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 2), (CIGAR.I, 3), (CIGAR.D, 2), (CIGAR.I, 2), (CIGAR.EQ, 8)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 4), (CIGAR.EQ, 8)], hgvs_standardize_cigar(read, rseq))

    def test_bubble_sort_indel_sections_drop_mismatch(self):
        rseq = 'ATAGGC' 'ATCT' 'A' 'CGA'   'AGCAT'     'ACGA' 'GATCGCTACG'
        #                ATCT   CTTTT                 TACGA
        qseq = 'ATCT' 'C'      'TT'     'TTT' 'ACGA' 'GATC'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence=qseq,
            cigar=[(CIGAR.EQ, 4), (CIGAR.X, 1), (CIGAR.D, 3), (CIGAR.I, 2), (CIGAR.D, 5), (CIGAR.I, 3), (CIGAR.EQ, 8)]
        )
        self.assertEqual(
            [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 8), (CIGAR.EQ, 9)], hgvs_standardize_cigar(read, rseq))

    def test_bubble_sort_indel_sections_drop_mismatch_with_hardclipping(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        # original
        # ATAGGCATCTACG   AA  CGAACGAGATCGCTACG
        #       ATCTC  TTT  TTCGAACG
        # expected
        # ATAGGCATCT      ACGAACGAACGAGATCGCTACG
        #       ATCTCTTTTT     CGAACG
        read = MockRead(
            'name',
            1,
            6,
            reference_name='1',
            query_sequence='ATCTCTTTTTCGAACG',
            cigar=[(CIGAR.H, 10), (CIGAR.EQ, 4), (CIGAR.X, 1), (CIGAR.D, 2), (CIGAR.I, 3), (CIGAR.D, 2), (CIGAR.I, 2), (CIGAR.EQ, 6)]
        )
        print(SamRead.deletion_sequences(read, {'1': MockObject(seq=ref)}))
        print(SamRead.insertion_sequences(read))
        print(read.query_sequence, len(read.query_sequence))
        self.assertEqual(
            [(CIGAR.H, 10), (CIGAR.EQ, 4), (CIGAR.I, 6), (CIGAR.D, 5), (CIGAR.EQ, 6)], hgvs_standardize_cigar(read, ref))

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

    def test_even_deletion_in_repeat(self):
        rseq = 'AAAGAAAAAAAAAAAAT' 'ATATATATATA' 'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        qseq = 'TTTTAAAAAAAAAAAAT' 'ATATATATATA'   'AATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        print(len(qseq) - 28)
        read = MockRead(
            'name', reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=2D64='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('4S24=2D53=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        print(SamRead.deletion_sequences(read, reference_genome))
        read.cigar = new_cigar
        print(SamRead.deletion_sequences(read, reference_genome))
        self.assertEqual(exp, new_cigar)

    def test_odd_deletion_in_repeat(self):
        rseq = 'AAAGAAAAAAAAAAAAT' 'ATATATATATA' 'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        qseq = 'TTTTAAAAAAAAAAAAT' 'ATATATATATA'    'ATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        print(len(qseq) - 28)
        read = MockRead(
            'name', reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=3D63='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('4S24=3D52=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        print(SamRead.deletion_sequences(read, reference_genome))
        read.cigar = new_cigar
        print(SamRead.deletion_sequences(read, reference_genome))
        self.assertEqual(exp, new_cigar)

    def test_unecessary_indel(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfghjklzxcvbnm'
        read = MockRead(
            'name', reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('13=1I1D12='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('26=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        self.assertEqual(exp, new_cigar)

    def test_unecessary_indel2(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfkghjklzxcvbnm'
        read = MockRead(
            'name', reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('13=2I1D12='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('14=1I12=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        self.assertEqual(exp, new_cigar)

    def test_unecessary_indel_end_match(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfkmkghjklzxcvbnm'
        read = MockRead(
            'name', reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('14=5I2D10='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('14=3I12=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        self.assertEqual(exp, new_cigar)

    def test_unecessary_indel_end_match2(self):
        rseq = 'GGGTGCAGTGGCTTACACCT' 'GTAATCCAAACACCTTGGGAGCCGCCCCCTGAG' 'CCTCCAGGCCCGGGACAGA'
        qseq = 'GGGTGCAGTGGCTTACACCT' 'CCAGG'                             'CCTCCAGGCCCGGGACAGA'
        read = MockRead(
            'name', reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('20=5I33D19='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('20=4I32D20=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        self.assertEqual(exp, new_cigar)

    def test_even_insertion_in_repeat(self):
        rseq = 'AAAGAAAAAAAAAAAAT' 'ATATATATATATA'   'AATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        qseq = 'TTTTAAAAAAAAAAAAT' 'ATATATATATATA' 'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        print(len(qseq) - 13 - 4)
        read = MockRead(
            'name', reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=2I66='),
            query_sequence=qseq
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('4S26=2I53=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        read.cigar = new_cigar
        self.assertEqual(exp, new_cigar)

    def test_deletion_repeat(self):
        qseq = (
            'GAGT'
            'GAGACTCTGT'
            'GAA'
            'AAAGAAAAAAAAAA'
            'A'
            'ATATATATATATATAAATATA'
            'C'
            'ATATTATGTATCAAATATATAT'
            'TATGTGTAATATACATCATGTATCAAATATATATTATGTATAATATACATCATATATCAAATATATATTATGTG'
        )
        # deleted reference: TATGTGTAATATACATCATGTATCAAA
        print(qseq[:76], qseq[76:])
        read = MockRead(
            'name', reference_name='11_86018001-86018500',
            reference_start=28,
            cigar=[
                (CIGAR.S, 4), (CIGAR.EQ, 10), (CIGAR.X, 3), (CIGAR.EQ, 14), (CIGAR.X, 1),
                (CIGAR.EQ, 21), (CIGAR.X, 1), (CIGAR.EQ, 22), (CIGAR.D, 27), (CIGAR.EQ, 74)
            ],
            query_sequence=qseq
        )
        expected_cigar = [
            (CIGAR.S, 4), (CIGAR.EQ, 10), (CIGAR.X, 3), (CIGAR.EQ, 14), (CIGAR.X, 1),
            (CIGAR.EQ, 21), (CIGAR.X, 1), (CIGAR.EQ, 22 + 30), (CIGAR.D, 27), (CIGAR.EQ, 74 - 30)
        ]
        std_cigar = hgvs_standardize_cigar(read, REFERENCE_GENOME[read.reference_name].seq)
        print(SamRead.deletion_sequences(read, REFERENCE_GENOME))
        read.cigar = std_cigar
        print(SamRead.deletion_sequences(read, REFERENCE_GENOME))
        self.assertEqual(expected_cigar, std_cigar)

    @timeout_decorator.timeout(1)
    def test_complex(self):
        qseq = (
            'TATTTGGAAATATTTGTAAGATAGATGTCTCTG' 'C'
            'CTCCTTCTGTTTCTGTCTCTGTCTCTTGCACTCTCTCTCTCCCTCTCTT'
            'TCTCTCTCTCTCTCTCTCTCTCTCTC'
            'TCTATATATATATATATA'
            'T' 'A' 'T' 'C' 'T'
            'ACACACACACACACACAC')
        rseq = (
            'TATTTGGAAATATTTGTAAGATAGATGTCTCTG' 'T'
            'CTCCTTCTGTTTCTGTCTCTGTCTCTTGCACTCTCTCTCTCCCTCTCTT'
            'TCTATATATATATATATA'
            'C' 'A' 'C'
            'ACACACACACACACACAC')
        read = MockRead(
            'name', reference_name='mock', reference_start=0, query_sequence=qseq,
            cigar=[
                (CIGAR.EQ, 33), (CIGAR.X, 1), (CIGAR.EQ, 49), (CIGAR.I, 26),
                (CIGAR.EQ, 18), (CIGAR.X, 1), (CIGAR.EQ, 1), (CIGAR.I, 1),
                (CIGAR.EQ, 1), (CIGAR.I, 1), (CIGAR.EQ, 18)]
        )
        print(rseq)
        print(read.query_sequence[:83], read.query_sequence[83 + 26: 83 + 26 + 20], read.query_sequence[83 + 26 + 22:])
        print(read.query_sequence)
        print(SamRead.insertion_sequences(read))
        new_cigar = [
            (CIGAR.EQ, 33), (CIGAR.X, 1), (CIGAR.EQ, 52), (CIGAR.I, 26),
            (CIGAR.EQ, 15), (CIGAR.X, 1), (CIGAR.EQ, 1), (CIGAR.I, 1),
            (CIGAR.EQ, 1), (CIGAR.I, 1), (CIGAR.EQ, 18)]
        std_cigar = hgvs_standardize_cigar(read, rseq)
        print(new_cigar)
        print(std_cigar)
        self.assertEqual(new_cigar, std_cigar)

    def test_deletion_partial_repeat(self):
        qseq = ('ATCTTAGCCAGGT'          'AGTTACATACATATC')
        rseq = ('ATCTTAGCCAGGT' 'AGCTAT' 'AGTTACATACATATC')
        read = MockRead(
            'name', reference_name='mock', reference_start=0, query_sequence=qseq,
            cigar=convert_string_to_cigar('13=6D15=')
        )
        self.assertEqual(convert_string_to_cigar('15=6D13='), hgvs_standardize_cigar(read, rseq))

    def test_indel_repeat(self):
        qseq = ('ATCTTAGCCAGGT' 'C'      'AGTTACATACATATC')
        rseq = ('ATCTTAGCCAGGT' 'AGCTAT' 'AGTTACATACATATC')
        print(qseq)
        print(rseq)
        read = MockRead(
            'name', reference_name='mock', reference_start=0, query_sequence=qseq,
            cigar=convert_string_to_cigar('13=1I6D15=')
        )
        self.assertEqual(convert_string_to_cigar('13=1I6D15='), hgvs_standardize_cigar(read, rseq))

    def test_shift_complex_indel(self):
        refseq = 'ATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAATGCTTTCTTTTATTAGTTGCCCTGTTTCAGATTCAGCTTTGTATCTATATCACCTGTTAATATGTGTGGACTCACAGAAATGATCATTGAGGGAATGCACCCTGTTTGGGTGTAAGTAGCTCAGGGAAAAAATCCTAG'
        read = MockRead(
            'name', reference_name='18',
            reference_start=40237946 - 40237890,
            query_sequence='AGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTGTATATATATATATATATATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAGTGCTTTCTTTTATTAGTGGCCCTG',
            cigar=convert_string_to_cigar('44M18I88M')
        )
        print(_read.convert_cigar_to_string(read.cigar))
        read.cigar = recompute_cigar_mismatch(read, refseq)
        self.assertEqual(convert_string_to_cigar('44=18I63=1X17=1X6='), read.cigar)
        print(_read.convert_cigar_to_string(read.cigar))
        read.cigar = hgvs_standardize_cigar(read, refseq)
        print(_read.convert_cigar_to_string(read.cigar))
        self.assertEqual(convert_string_to_cigar('45=18I62=1X17=1X6='), read.cigar)


class TestMergeInternalEvents(unittest.TestCase):

    def test_small_exact_match(self):
        cigar = convert_string_to_cigar('283M17506D5M21275D596M17506D5M21275D313M')
        # [(0, 283), (2, 17506), (0, 5), (2, 21275), (0, 596), (2, 17506), (0, 5), (2, 21275), (0, 313)]
        new_cigar = merge_internal_events(cigar, 20, 15)
        exp = [
            (CIGAR.M, 283),
            (CIGAR.I, 5),
            (CIGAR.D, 17506 + 21275 + 5),
            (CIGAR.M, 596),
            (CIGAR.I, 5),
            (CIGAR.D, 17506 + 21275 + 5),
            (CIGAR.M, 313)
        ]
        self.assertEqual(exp, new_cigar)


class TestConvertStringToCigar(unittest.TestCase):
    def test(self):
        string = '283M' '17506D' '5M' '21275D' '596M' '17506D' '5M' '21275D' '313M'
        exp = [
            (CIGAR.M, 283),
            (CIGAR.D, 17506),
            (CIGAR.M, 5),
            (CIGAR.D, 21275),
            (CIGAR.M, 596),
            (CIGAR.D, 17506),
            (CIGAR.M, 5),
            (CIGAR.D, 21275),
            (CIGAR.M, 313)
        ]
        self.assertEqual(exp, convert_string_to_cigar(string))


class TestGetSequences(unittest.TestCase):

    def setUp(self):
        self.reference_genome = {'1': MockObject(seq='abcdefghijklmnopqrstuvwxyz')}

    def test_deletions(self):
        exp = ['cde', 'nopq']
        read = MockRead(
            reference_start=0, reference_name='1', query_sequence='',
            cigar=convert_string_to_cigar('2=3D8=4D9=')
        )
        self.assertEqual(exp, SamRead.deletion_sequences(read, self.reference_genome))

    def test_insertions(self):
        exp = ['kkk', 'kkkk']
        read = MockRead(
            reference_start=0, reference_name='1', query_sequence='abcdekkkfghijklmnopqkkkkrstuvwxyz',
            cigar=convert_string_to_cigar('5=3I12=4I9=')
        )
        self.assertEqual(exp, SamRead.insertion_sequences(read))
