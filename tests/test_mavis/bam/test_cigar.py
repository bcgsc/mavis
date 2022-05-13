import warnings

import pytest
import timeout_decorator
from mavis.annotate.file_io import load_reference_genome
from mavis.bam import read as _read
from mavis.bam.cigar import (
    alignment_matches,
    compute,
    convert_for_igv,
    convert_string_to_cigar,
    extend_softclipping,
    hgvs_standardize_cigar,
    join,
    longest_fuzzy_match,
    match_percent,
    merge_internal_events,
    recompute_cigar_mismatch,
    score,
)
from mavis.bam.read import SamRead
from mavis.constants import CIGAR

from ...util import get_data
from ..mock import MockObject, MockRead


@pytest.fixture(scope='module')
def mock_reference_genome():
    warnings.simplefilter('ignore')
    ref = load_reference_genome(get_data('mock_reference_genome.fa'))
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != ref['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    return ref


@pytest.mark.parametrize(
    'reference_start,query_sequence,cigar_in,cigar_out',
    [
        [1456, 'CCCAAACAAC' 'TATAAATTTT' 'GTAATACCTA' 'GAACAATATA' 'AATAT', '45M', '45='],
        [1456, 'CCCAAACAAC' 'TATAAATTTT' 'GTAATACCTA' 'GAACAATATA' 'AATAT', '20H45M', '20H45='],
        [
            1456,
            'TATA' 'CCCAAACAAC' 'TATAAATTTT' 'GTAATACCTA' 'GAACAATATA' 'AATAT',
            '4S10M10D10I25M',
            '4S10=10D10I25=',
        ],
        [
            1452,
            'CAGC' 'CCCAAACAAC' 'TATAAATTTT' 'GTAATACCTA' 'GAACAATATA' 'AATAT',
            '4X10M10D10I25M',
            '4X10=10D10I25=',
        ],
        [
            1452,
            'CAGC' 'CCCAAACAAC' 'TATAAATTTT' 'GTAATACCTA' 'GAACAATATA' 'AATAT',
            '14M10D10I25M',
            '4X10=10D10I25=',
        ],
    ],
)
def test_recompute_cigar_mismatch(
    reference_start, query_sequence, cigar_in, cigar_out, mock_reference_genome
):
    r = MockRead(
        reference_start=reference_start,
        query_sequence=query_sequence,
        cigar=convert_string_to_cigar(cigar_in),
    )
    assert recompute_cigar_mismatch(r, mock_reference_genome['fake']) == convert_string_to_cigar(
        cigar_out
    )


class TestExtendSoftclipping:
    def test_softclipped_right(self):
        c = convert_string_to_cigar('70=2X1=8X4=1X1=4X1=6X1=4X1=4X2=5X3=3X1=4X1=3X1=14X1=1X2=1S')
        cnew, prefix = extend_softclipping(c, 6)
        assert prefix == 0
        assert cnew == convert_string_to_cigar('70=80S')


@pytest.mark.parametrize(
    'input_cigars,output_cigar',
    [
        [['10M10X10X'], '10M20X'],
        [['10X10M10X'], '10X10M10X'],
        [['10M10X10X', '10X10M10X'], '10M30X10M10X'],
        [['1S2S5=7X2=5X28=1X99='], '3S5=7X2=5X28=1X99='],
        [['10H10M10X10X'], '10H10M20X'],
    ],
)
def test_join(input_cigars, output_cigar):
    assert join(*[convert_string_to_cigar(c) for c in input_cigars]) == convert_string_to_cigar(
        output_cigar
    )


@pytest.mark.parametrize(
    'seq1,seq2,opt,cigar_string,ref_start',
    [
        ['GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------', {}, '7S8=', 7],
        ['GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------', {}, '5S8=', 7],
        [
            'GTGAGTAAATTCAACATCGTTTTT',
            '--CTTAGAATTCAAC---------',
            dict(force_softclipping=False),
            '5S8=',
            7,
        ],
        [
            'GTGAGTAAATTC--CATCGTTTTT',
            '--CTTAGAATTCAAC---------',
            dict(force_softclipping=False),
            '5S5=2I1=',
            7,
        ],
        ['CCTG', 'CCGT', dict(min_exact_to_stop_softclipping=10), '2=2X', 0],
        [
            '--GAGTAAATTCAACATCGTTTTT',
            '--CTTAGAATTCAAC---------',
            dict(force_softclipping=False),
            '5S8=',
            5,
        ],
    ],
)
def test_compute(seq1, seq2, opt, cigar_string, ref_start):
    assert compute(seq1, seq2, **opt) == (convert_string_to_cigar(cigar_string), ref_start)


def test_compute_error():
    with pytest.raises(AttributeError):
        compute('CCTG', 'CCG')


class TestCigarTools:
    def test_alignment_matches(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        assert alignment_matches(c) == 30

    def test_longest_fuzzy_match(self):
        c = [
            (CIGAR.S, 10),
            (CIGAR.EQ, 1),
            (CIGAR.X, 4),
            (CIGAR.EQ, 10),
            (CIGAR.I, 3),
            (CIGAR.EQ, 5),
        ]
        assert longest_fuzzy_match(c, 1) == 15
        assert longest_fuzzy_match(c, 0) == 10
        assert longest_fuzzy_match(c, 2) == 16
        assert longest_fuzzy_match(c, 4) == 16

    def test_score(self):
        c = [
            (CIGAR.S, 10),
            (CIGAR.EQ, 1),
            (CIGAR.X, 4),
            (CIGAR.EQ, 10),
            (CIGAR.I, 3),
            (CIGAR.EQ, 5),
        ]
        assert score(c) == 22

    def test_score_error(self):
        with pytest.raises(AssertionError):
            c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (99, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
            score(c)

    def test_match_percent(self):
        c = [
            (CIGAR.S, 10),
            (CIGAR.EQ, 1),
            (CIGAR.X, 4),
            (CIGAR.EQ, 10),
            (CIGAR.I, 3),
            (CIGAR.EQ, 5),
        ]
        assert match_percent(c) == 0.8
        with pytest.raises(AttributeError):
            match_percent([(CIGAR.M, 100)])
        with pytest.raises(AttributeError):
            match_percent([(CIGAR.S, 100)])

    def test_convert_for_igv(self):
        c = [(CIGAR.M, 10), (CIGAR.EQ, 10), (CIGAR.X, 10)]
        assert convert_for_igv(c) == [(CIGAR.M, 30)]


class TestHgvsStandardizeCigars:
    def no_change_aligned(self):
        ref = 'AAATTTGGGCCCAATT'
        read = MockRead('name', '1', 1, cigar=[(CIGAR.M, 10)], query_sequence='AAATTTGGGC')
        assert hgvs_standardize_cigar(read, ref) == [(CIGAR.M, 10)]

    def no_change_proper_indel(self):
        ref = 'ATAGGC' 'ATCTACGAG' 'ATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCTAC' 'CCC' 'ATCG',
            cigar=[(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)],
        )
        assert [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)] == hgvs_standardize_cigar(
            read, ref
        )

    def ins_after_deletion(self):
        ref = 'ATAGGC' 'ATCTACGAG' 'ATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCTAC' 'CCC' 'ATCG',
            cigar=[(CIGAR.EQ, 6), (CIGAR.D, 3), (CIGAR.I, 3), (CIGAR.EQ, 4)],
        )
        assert [(CIGAR.EQ, 6), (CIGAR.I, 3), (CIGAR.D, 3), (CIGAR.EQ, 4)] == hgvs_standardize_cigar(
            read, ref
        )

    def test_insertion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.I, 4), (CIGAR.EQ, 8)],
        )
        assert [(CIGAR.EQ, 8), (CIGAR.I, 4), (CIGAR.EQ, 4)] == hgvs_standardize_cigar(read, ref)

    def test_deletion_in_repeat(self):
        ref = 'ATAGGC' 'ATCT' 'ACGA' 'ACGA' 'ACGA' 'GATCGCTACG'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence='ATCT' 'ACGA' 'ACGA' 'GATC',
            cigar=[(CIGAR.EQ, 4), (CIGAR.D, 4), (CIGAR.EQ, 12)],
        )
        assert [(CIGAR.EQ, 12), (CIGAR.D, 4), (CIGAR.EQ, 4)] == hgvs_standardize_cigar(read, ref)

    def test_bubble_sort_indel_sections(self):
        rseq = 'ATAGGC' 'ATCT' 'GG' 'GA' 'GCGA' 'GATCGCTACG'
        qseq = 'ATCT' 'TTT' 'TT' 'GCGA' 'GATC'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence=qseq,
            cigar=[
                (CIGAR.EQ, 4),
                (CIGAR.D, 2),
                (CIGAR.I, 3),
                (CIGAR.D, 2),
                (CIGAR.I, 2),
                (CIGAR.EQ, 8),
            ],
        )
        assert [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 4), (CIGAR.EQ, 8)] == hgvs_standardize_cigar(
            read, rseq
        )

    def test_bubble_sort_indel_sections_drop_mismatch(self):
        rseq = 'ATAGGC' 'ATCT' 'A' 'CGA' 'AGCAT' 'ACGA' 'GATCGCTACG'
        #                ATCT   CTTTT                 TACGA
        qseq = 'ATCT' 'C' 'TT' 'TTT' 'ACGA' 'GATC'
        read = MockRead(
            'name',
            1,
            6,
            query_sequence=qseq,
            cigar=[
                (CIGAR.EQ, 4),
                (CIGAR.X, 1),
                (CIGAR.D, 3),
                (CIGAR.I, 2),
                (CIGAR.D, 5),
                (CIGAR.I, 3),
                (CIGAR.EQ, 8),
            ],
        )
        assert [(CIGAR.EQ, 4), (CIGAR.I, 5), (CIGAR.D, 8), (CIGAR.EQ, 9)] == hgvs_standardize_cigar(
            read, rseq
        )

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
            cigar=[
                (CIGAR.H, 10),
                (CIGAR.EQ, 4),
                (CIGAR.X, 1),
                (CIGAR.D, 2),
                (CIGAR.I, 3),
                (CIGAR.D, 2),
                (CIGAR.I, 2),
                (CIGAR.EQ, 6),
            ],
        )
        print(SamRead.deletion_sequences(read, {'1': MockObject(seq=ref)}))
        print(SamRead.insertion_sequences(read))
        print(read.query_sequence, len(read.query_sequence))
        assert [
            (CIGAR.H, 10),
            (CIGAR.EQ, 4),
            (CIGAR.I, 6),
            (CIGAR.D, 5),
            (CIGAR.EQ, 6),
        ] == hgvs_standardize_cigar(read, ref)

    def test_homopolymer_even_odd(self):
        ref = 'ATCGAGAT' + 'A' * 15 + 'TCGAGAT'
        read = MockRead(
            'name',
            1,
            1,
            query_sequence='ATCGAGATA' + 'A' * 12 + 'TCGAGAT',
            cigar=[(CIGAR.EQ, 8), (CIGAR.D, 2), (CIGAR.EQ, 20)],
        )
        assert [(CIGAR.EQ, 9 + 12), (CIGAR.D, 2), (CIGAR.EQ, 7)] == hgvs_standardize_cigar(
            read, ref
        )
        ref = (
            'CCCCGGCTCATGTCTGGTTTTGTTTTCCGGGGGCGGGGGGGCTCCCTGGGGATGATGGTGATTTTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAA'
            'TGAGGCAGAGACAATGTGGGGAGCGAGAGAGGGGAAAAGGACGGGGGAGG'
        )

        read = MockRead(
            'name',
            '1',
            0,
            149,
            query_sequence=(
                'CCCCGGCTCATGTCTGGTTTTGTTTTCCGGGGGCGGGGGGGCTCCCTGGGGATGATGGTGATTTTTTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAA'
                'TGAGGCAGAGACAATGTGGGGAGCGAGAGAGGGGAAAAGGACGGGGGAGG'
            ),
            cigar=[(CIGAR.EQ, 61), (CIGAR.I, 2), (CIGAR.EQ, 87)],
        )
        assert [(CIGAR.EQ, 61 + 15), (CIGAR.I, 2), (CIGAR.EQ, 87 - 15)] == hgvs_standardize_cigar(
            read, ref
        )

        ref = (
            'CCTCCTCGGTCGGGCAGATCTTTCAGAAGCAGGAGCCCAGGATCATGTCTGGTTTTGTTTTCCGAGGGCGAGGGGGCTCCCTGAGGATGATGGTGATTT'
            'TTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAATGAGGCAGAGACA'
        )

        read = MockRead(
            'name',
            '1',
            0,
            149,
            query_sequence=(
                'CCCCTCCTCGGTCGGGCAGATCTTTCAGAAGCAGGAGCCCAGGATCATGTCTGGTTTTGTTTTCCGAGGGCGAGGGGGCTCCCTGAGGATGATGGTGATTTT'
                'TTTTTTTTTTTTTAATCCTCAACTAGGAGAGAAAATGAGGCAGAGACA'
            ),
            cigar=[(CIGAR.S, 2), (CIGAR.EQ, 96), (CIGAR.I, 2), (CIGAR.EQ, 50)],
        )
        assert [
            (CIGAR.S, 2),
            (CIGAR.EQ, 96 + 15),
            (CIGAR.I, 2),
            (CIGAR.EQ, 50 - 15),
        ] == hgvs_standardize_cigar(read, ref)

    def test_even_deletion_in_repeat(self):
        rseq = (
            'AAAGAAAAAAAAAAAAT'
            'ATATATATATA'
            'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        qseq = (
            'TTTTAAAAAAAAAAAAT'
            'ATATATATATA'
            'AATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        print(len(qseq) - 28)
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=2D64='),
            query_sequence=qseq,
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('4S24=2D53=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        print(SamRead.deletion_sequences(read, reference_genome))
        read.cigar = new_cigar
        print(SamRead.deletion_sequences(read, reference_genome))
        assert new_cigar == exp

    def test_odd_deletion_in_repeat(self):
        rseq = (
            'AAAGAAAAAAAAAAAAT'
            'ATATATATATA'
            'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        qseq = (
            'TTTTAAAAAAAAAAAAT' 'ATATATATATA' 'ATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        print(len(qseq) - 28)
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=3D63='),
            query_sequence=qseq,
        )
        reference_genome = {'1': MockObject(seq=rseq)}
        exp = convert_string_to_cigar('4S24=3D52=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        print(SamRead.deletion_sequences(read, reference_genome))
        read.cigar = new_cigar
        print(SamRead.deletion_sequences(read, reference_genome))
        assert new_cigar == exp

    def test_unecessary_indel(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfghjklzxcvbnm'
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('13=1I1D12='),
            query_sequence=qseq,
        )
        exp = convert_string_to_cigar('26=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        assert new_cigar == exp

    def test_unecessary_indel2(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfkghjklzxcvbnm'
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('13=2I1D12='),
            query_sequence=qseq,
        )
        exp = convert_string_to_cigar('14=1I12=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        assert new_cigar == exp

    def test_unecessary_indel_end_match(self):
        rseq = 'qwertyuiopasdfghjklzxcvbnm'
        qseq = 'qwertyuiopasdfkmkghjklzxcvbnm'
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('14=5I2D10='),
            query_sequence=qseq,
        )
        exp = convert_string_to_cigar('14=3I12=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        assert new_cigar == exp

    def test_unecessary_indel_end_match2(self):
        rseq = 'GGGTGCAGTGGCTTACACCT' 'GTAATCCAAACACCTTGGGAGCCGCCCCCTGAG' 'CCTCCAGGCCCGGGACAGA'
        qseq = 'GGGTGCAGTGGCTTACACCT' 'CCAGG' 'CCTCCAGGCCCGGGACAGA'
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=0,
            cigar=convert_string_to_cigar('20=5I33D19='),
            query_sequence=qseq,
        )
        exp = convert_string_to_cigar('20=4I32D20=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        assert new_cigar == exp

    def test_even_insertion_in_repeat(self):
        rseq = (
            'AAAGAAAAAAAAAAAAT'
            'ATATATATATATA'
            'AATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        qseq = (
            'TTTTAAAAAAAAAAAAT'
            'ATATATATATATA'
            'TAAATATACATATTATGTATCAAATATATATTATGTGTAATATACATCATGTATC'
        )
        print(len(qseq) - 13 - 4)
        read = MockRead(
            'name',
            reference_name='1',
            reference_start=4,
            cigar=convert_string_to_cigar('4S13=2I66='),
            query_sequence=qseq,
        )
        exp = convert_string_to_cigar('4S26=2I53=')
        new_cigar = hgvs_standardize_cigar(read, rseq)
        read.cigar = new_cigar
        assert new_cigar == exp

    def test_deletion_repeat(self, mock_reference_genome):
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
            'name',
            reference_name='11_86018001-86018500',
            reference_start=28,
            cigar=[
                (CIGAR.S, 4),
                (CIGAR.EQ, 10),
                (CIGAR.X, 3),
                (CIGAR.EQ, 14),
                (CIGAR.X, 1),
                (CIGAR.EQ, 21),
                (CIGAR.X, 1),
                (CIGAR.EQ, 22),
                (CIGAR.D, 27),
                (CIGAR.EQ, 74),
            ],
            query_sequence=qseq,
        )
        expected_cigar = [
            (CIGAR.S, 4),
            (CIGAR.EQ, 10),
            (CIGAR.X, 3),
            (CIGAR.EQ, 14),
            (CIGAR.X, 1),
            (CIGAR.EQ, 21),
            (CIGAR.X, 1),
            (CIGAR.EQ, 22 + 30),
            (CIGAR.D, 27),
            (CIGAR.EQ, 74 - 30),
        ]
        std_cigar = hgvs_standardize_cigar(read, mock_reference_genome[read.reference_name].seq)
        print(SamRead.deletion_sequences(read, mock_reference_genome))
        read.cigar = std_cigar
        print(SamRead.deletion_sequences(read, mock_reference_genome))
        assert std_cigar == expected_cigar

    @timeout_decorator.timeout(1)
    def test_complex(self):
        qseq = (
            'TATTTGGAAATATTTGTAAGATAGATGTCTCTG'
            'C'
            'CTCCTTCTGTTTCTGTCTCTGTCTCTTGCACTCTCTCTCTCCCTCTCTT'
            'TCTCTCTCTCTCTCTCTCTCTCTCTC'
            'TCTATATATATATATATA'
            'T'
            'A'
            'T'
            'C'
            'T'
            'ACACACACACACACACAC'
        )
        rseq = (
            'TATTTGGAAATATTTGTAAGATAGATGTCTCTG'
            'T'
            'CTCCTTCTGTTTCTGTCTCTGTCTCTTGCACTCTCTCTCTCCCTCTCTT'
            'TCTATATATATATATATA'
            'C'
            'A'
            'C'
            'ACACACACACACACACAC'
        )
        read = MockRead(
            'name',
            reference_name='mock',
            reference_start=0,
            query_sequence=qseq,
            cigar=[
                (CIGAR.EQ, 33),
                (CIGAR.X, 1),
                (CIGAR.EQ, 49),
                (CIGAR.I, 26),
                (CIGAR.EQ, 18),
                (CIGAR.X, 1),
                (CIGAR.EQ, 1),
                (CIGAR.I, 1),
                (CIGAR.EQ, 1),
                (CIGAR.I, 1),
                (CIGAR.EQ, 18),
            ],
        )
        print(rseq)
        print(
            read.query_sequence[:83],
            read.query_sequence[83 + 26 : 83 + 26 + 20],
            read.query_sequence[83 + 26 + 22 :],
        )
        print(read.query_sequence)
        print(SamRead.insertion_sequences(read))
        new_cigar = [
            (CIGAR.EQ, 33),
            (CIGAR.X, 1),
            (CIGAR.EQ, 52),
            (CIGAR.I, 26),
            (CIGAR.EQ, 15),
            (CIGAR.X, 1),
            (CIGAR.EQ, 1),
            (CIGAR.I, 1),
            (CIGAR.EQ, 1),
            (CIGAR.I, 1),
            (CIGAR.EQ, 18),
        ]
        std_cigar = hgvs_standardize_cigar(read, rseq)
        print(new_cigar)
        print(std_cigar)
        assert std_cigar == new_cigar

    def test_deletion_partial_repeat(self):
        qseq = 'ATCTTAGCCAGGT' 'AGTTACATACATATC'
        rseq = 'ATCTTAGCCAGGT' 'AGCTAT' 'AGTTACATACATATC'
        read = MockRead(
            'name',
            reference_name='mock',
            reference_start=0,
            query_sequence=qseq,
            cigar=convert_string_to_cigar('13=6D15='),
        )
        assert convert_string_to_cigar('15=6D13=') == hgvs_standardize_cigar(read, rseq)

    def test_indel_repeat(self):
        qseq = 'ATCTTAGCCAGGT' 'C' 'AGTTACATACATATC'
        rseq = 'ATCTTAGCCAGGT' 'AGCTAT' 'AGTTACATACATATC'
        print(qseq)
        print(rseq)
        read = MockRead(
            'name',
            reference_name='mock',
            reference_start=0,
            query_sequence=qseq,
            cigar=convert_string_to_cigar('13=1I6D15='),
        )
        assert convert_string_to_cigar('13=1I6D15=') == hgvs_standardize_cigar(read, rseq)

    def test_shift_complex_indel(self):
        refseq = 'ATATATCTATTTTTTTCTTTCTTTTTTTTACTTTCATTAAGTGCCACTAAAAAATTAGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAATGCTTTCTTTTATTAGTTGCCCTGTTTCAGATTCAGCTTTGTATCTATATCACCTGTTAATATGTGTGGACTCACAGAAATGATCATTGAGGGAATGCACCCTGTTTGGGTGTAAGTAGCTCAGGGAAAAAATCCTAG'
        read = MockRead(
            'name',
            reference_name='18',
            reference_start=40237946 - 40237890,
            query_sequence='AGGTTCAATTAAACTTTATTAATCTCTTCTGAGTTTTGATTGAGTGTATATATATATATATATATATATATATATATACCCAGTTTCAAGCAGGTATCTGCCTTTAAAGATAAGAGACCTCCTAAGTGCTTTCTTTTATTAGTGGCCCTG',
            cigar=convert_string_to_cigar('44M18I88M'),
        )
        print(_read.convert_cigar_to_string(read.cigar))
        read.cigar = recompute_cigar_mismatch(read, refseq)
        assert read.cigar == convert_string_to_cigar('44=18I63=1X17=1X6=')
        print(_read.convert_cigar_to_string(read.cigar))
        read.cigar = hgvs_standardize_cigar(read, refseq)
        print(_read.convert_cigar_to_string(read.cigar))
        assert read.cigar == convert_string_to_cigar('45=18I62=1X17=1X6=')


class TestMergeInternalEvents:
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
            (CIGAR.M, 313),
        ]
        assert new_cigar == exp


@pytest.mark.parametrize(
    'string,cigartuples',
    [
        [
            '283M' '17506D' '5M' '21275D' '596M' '17506D' '5M' '21275D' '313M',
            [
                (CIGAR.M, 283),
                (CIGAR.D, 17506),
                (CIGAR.M, 5),
                (CIGAR.D, 21275),
                (CIGAR.M, 596),
                (CIGAR.D, 17506),
                (CIGAR.M, 5),
                (CIGAR.D, 21275),
                (CIGAR.M, 313),
            ],
        ],
        [
            '42H25M',
            [
                (CIGAR.H, 42),
                (CIGAR.M, 25),
            ],
        ],
    ],
)
def test_convert_string_to_cigar(string, cigartuples):
    assert convert_string_to_cigar(string) == cigartuples


class TestGetSequences:
    def test_deletions(self):
        exp = ['cde', 'nopq']
        read = MockRead(
            reference_start=0,
            reference_name='1',
            query_sequence='',
            cigar=convert_string_to_cigar('2=3D8=4D9='),
        )
        assert (
            SamRead.deletion_sequences(read, {'1': MockObject(seq='abcdefghijklmnopqrstuvwxyz')})
            == exp
        )

    def test_insertions(self):
        exp = ['kkk', 'kkkk']
        read = MockRead(
            reference_start=0,
            reference_name='1',
            query_sequence='abcdekkkfghijklmnopqkkkkrstuvwxyz',
            cigar=convert_string_to_cigar('5=3I12=4I9='),
        )
        assert SamRead.insertion_sequences(read) == exp
