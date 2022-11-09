import pytest
from Bio import SeqIO

import mavis.bam.cigar as _cigar
from mavis.annotate.file_io import load_reference_genome
from mavis.bam.cache import BamCache
from mavis.constants import CIGAR, reverse_complement
from mavis.interval import Interval
from mavis.validate.align import query_coverage_interval
from mavis.validate.blat import Blat

from ..util import get_data
from .mock import MockBamFileHandle, MockLongString, MockObject

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))


@pytest.fixture
def cache():
    return BamCache(MockBamFileHandle({'Y': 23, 'fake': 0, 'reference3': 3, '14': 13}))


class TestBlat:
    def test_read_pslx(self):
        mapping = {}
        for record in SeqIO.parse(get_data('blat_input.fa'), 'fasta'):
            mapping[record.id] = record.seq
        header, rows = Blat.read_pslx(get_data('blat_output.pslx'), mapping)
        assert len(rows) == 11067
        expect_pslx_header = [
            'match',
            'mismatch',
            'repmatch',
            'ncount',
            'qgap_count',
            'qgap_bases',
            'tgap_count',
            'tgap_bases',
            'strand',
            'qname',
            'qsize',
            'qstart',
            'qend',
            'tname',
            'tsize',
            'tstart',
            'tend',
            'block_count',
            'block_sizes',
            'qstarts',
            'tstarts',
            'qseqs',
            'tseqs',
        ]
        assert header == expect_pslx_header

    def test_pslx_row_to_pysam_single_block(self, cache):
        pslx_row = {
            'score': 20,
            'tseqs': ['AATACCAAATACATGATATA'],
            'tstarts': [3432307],
            'tstart': 3432307,
            'block_sizes': [20],
            'qname': 'seq1',
            'tname': 'Y',
            'qstarts': [93],
            'strand': '+',
            'qseqs': ['AATACCATACATGATATA'],
            'percent_ident': 100.0,
            'qseq_full': 'AGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTT'
            'AGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCGCCCGCCTCGGC',
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, None)
        assert read.reference_id == 23
        assert query_coverage_interval(read) == Interval(93, 112)

    def test_pslx_row_to_pysam_full_reverse(self, cache):
        pslx_row = {
            'match': 128,
            'mismatch': 0,
            'repmatch': 0,
            'ncount': 0,
            'qgap_count': 0,
            'qgap_bases': 0,
            'tgap_count': 0,
            'tgap_bases': 0,
            'strand': '-',
            'qname': 'seq1',
            'tname': 'reference3',
            'tsize': 3711,
            'block_sizes': [128],
            'qstarts': [117],
            'tstarts': [2187],
            'qseqs': [
                'TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'
            ],
            'tseqs': [
                'TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'
            ],
            '_index': 1,
            'score': 128,
            'percent_ident': 100.0,
            'qseq_full': 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT',
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, None)
        assert read.reference_id == 3
        assert read.cigar == [(CIGAR.S, 117), (CIGAR.M, 128)]
        assert read.reference_start == 2187
        assert query_coverage_interval(read) == Interval(117, 244)

    def test_pslx_row_to_pysam_simple(self, cache):
        pslx_row = {
            'tstarts': [950],
            'block_sizes': [53],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0],
            'strand': '+',
            'score': 0,
            'qseq_full': 'ATCTAATAACTTGATCAATA' 'TCTGTGATTATATTTTCATT' 'GCCTTCCAATTTT',
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, None)
        assert read.reference_id == 0
        assert query_coverage_interval(read) == Interval(0, 52)
        assert read.reference_start == 950
        assert read.reference_end == 1003
        assert read.cigar == [(CIGAR.M, 53)]

    def test_pslx_row_to_pysam_simple_with_reference(self, cache):
        pslx_row = {
            'tstarts': [950],
            'block_sizes': [53],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0],
            'strand': '+',
            'score': 0,
            'qseq_full': 'ATCTAATAACTTGATCAATA' 'TCTGTGATTATATTTTCATT' 'GCCTTCCAATTTT',
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, REFERENCE_GENOME)
        assert read.reference_id == 0
        assert query_coverage_interval(read) == Interval(0, 52)
        assert read.reference_start == 950
        assert read.reference_end == 1003
        assert read.cigar == [(CIGAR.EQ, 53)]

    def test_pslx_row_to_pysam_gapped_alignment(self, cache):
        pslx_row = {
            'block_count': 1,
            'tstarts': [950, 7233],
            'block_sizes': [47, 100],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0, 47],
            'strand': '+',
            'qseq_full': 'ATCTAATAACTTGATCAATA'
            'TCTGTGATTATATTTTCATT'
            'GCCTTCC'
            'AATTTTGCAGATTATAAGAT'
            'CAATAGATATTTATTGTAAA'
            'ATGCACAAATAGTGCAACAT'
            'TTCTTAAAGTAGACCGTGAA'
            'ATACTTCATGTTGCCATGTT',
            'score': 1,
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, None)
        assert read.reference_id == 0
        assert query_coverage_interval(read) == Interval(0, 146)
        assert read.reference_start == 950
        assert read.cigar == [(CIGAR.M, 47), (CIGAR.D, 6236), (CIGAR.M, 100)]

    def test_pslx_row_to_pysam_gapped_alignment_with_reference(self, cache):
        pslx_row = {
            'block_count': 1,
            'tstarts': [950, 7233],
            'block_sizes': [47, 100],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0, 47],
            'strand': '+',
            'qseq_full': 'ATCTAATAACTTGATCAATA'
            'TCTGTGATTATATTTTCATT'
            'GCCTTCC'
            'AATTTTGCAGATTATAAGAT'
            'CAATAGATATTTATTGTAAA'
            'ATGCACAAATAGTGCAACAT'
            'TTCTTAAAGTAGACCGTGAA'
            'ATACTTCATGTTGCCATGTT',
            'score': 1,
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, REFERENCE_GENOME)
        assert read.reference_id == 0
        assert query_coverage_interval(read) == Interval(0, 146)
        assert read.reference_start == 950
        assert read.cigar == [(CIGAR.EQ, 53), (CIGAR.D, 6236), (CIGAR.EQ, 94)]

    def test_pslx_row_to_pysam_revcomp_deletion(self, cache):
        pslx_row = {
            'block_count': 2,
            'tstarts': [2205, 2281],
            'block_sizes': [50, 34],
            'qname': 'seq1',
            'tname': 'reference3',
            'qstarts': [0, 50],
            'strand': '-',
            'qseq_full': 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTA',
            'score': 1,
            'qseqs': [
                'TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA',
                'CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG',
            ],
            'tseqs': [
                'TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA',
                'CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG',
            ],
        }
        read = Blat.pslx_row_to_pysam(pslx_row, cache, REFERENCE_GENOME)
        assert read.reference_id == 3
        assert query_coverage_interval(read) == Interval(0, 83)
        assert read.reference_start == 2205
        assert read.cigar == [(CIGAR.EQ, 51), (CIGAR.D, 26), (CIGAR.EQ, 33)]
        assert read.query_sequence[0:50] == 'TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA'
        assert read.query_sequence[50:] == 'CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'

    def test_pslx_row_to_pysam_inversion(self, cache):
        s = 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'
        # first part of the inversion
        pslx_row = {
            'block_count': 1,
            'tstarts': [1114],
            'block_sizes': [120],
            'qname': 'seq1',
            'tname': 'reference3',
            'qstarts': [125],
            'strand': '+',
            'qseq_full': s,
            'score': 1,
            'qseqs': [
                'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGG'
                'TTTTCATTTCTGTATGTTAAT'
            ],
            'tseqs': [
                'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGG'
                'TTTTCATTTCTGTATGTTAAT'
            ],
        }
        read1 = Blat.pslx_row_to_pysam(pslx_row, cache, REFERENCE_GENOME)
        assert read1.reference_id == 3
        assert query_coverage_interval(read1) == Interval(125, 244)
        assert read1.reference_start == 1114
        assert read1.cigar == [(CIGAR.S, 125), (CIGAR.EQ, 120)]

        # second part of the inversion
        pslx_row = {
            'block_count': 1,
            'tstarts': [2187],
            'block_sizes': [128],
            'qname': 'seq1',
            'tname': 'reference3',
            'qstarts': [117],
            'strand': '-',
            'qseq_full': s,
            'score': 1,
            'qseqs': [
                'TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAAT'
                'TCTGTGTTTACAGGGCTTTCATGCTCAG'
            ],
            'tseqs': [
                'TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAAT'
                'TCTGTGTTTACAGGGCTTTCATGCTCAG'
            ],
        }
        read2 = Blat.pslx_row_to_pysam(pslx_row, cache, REFERENCE_GENOME)
        assert read2.reference_id == 3
        assert read2.reference_start == 2187
        assert read2.cigar == [(CIGAR.S, 117), (CIGAR.EQ, 128)]
        assert query_coverage_interval(read2) == Interval(117, 244)
        assert reverse_complement(read2.query_sequence) == read1.query_sequence
        # test that this is selected for duplication or insertion evidence

    def test_pslx_row_to_pysam_duplication(self, cache):
        reference = {
            '14': MockObject(
                seq=MockLongString(
                    'TTCTTCCATGCCCCCTAATCATGGCCACATTGTATCAGCCTGAGCATGAGCAACAGCACCATGGCCACATACGGGAATGGGCCTCATTGGTGTAATATTTGGCAGATTCTCTCCACACCCCCCGTGGCGGTCTGGCTTACTGTTAAGAAGGGTAACCTTAAAAAATACATTTCCCACTCCAGAAAATACTCATATGTGGCCTGTTAGCAGCACAAGAAGGGTGAAAGCAATGCCCATTCCTGCCTCCCTCCCCCTGCTCACCTCCACGTCCCTGTTTGCCCCTTTGTAGGTGAAGTGAGTATATTCAGCGTCTTCATGGCAGGGGAGAGGGTGTATTAATCCGTCTATGTCCGCTGGAAAGGCAGTCTCTGAGCGGGCCACAAGGGTTCAGCCATGGCCCATCCAATAACCTTTTTGATGACTTGGATGAAGAGACAAACATTCCAACCACATTCAAAGATCCAGACCTCCAAAGTGTGGCTCATTTGGTAGATAATGGAATTATATTTGGAAAGCATTTCCCGCAGCTGGGATGATGGGTCAAAAACAGATAGCATTTTACCAGATCATATTTGTGTGTGTGTGTGTGCGCGCGTGTGTGTGTGTGTGTGTGTGTGTTTTAAATTCAGTTTCCCAACTACAGGATG',
                    offset=73014463,
                )
            )
        }
        pslx_row = {
            'block_count': 2,
            'tstarts': [73014606, 73014747],
            'block_sizes': [141, 30],
            'qname': '',
            'tname': '14',
            'qstarts': [0, 239],
            'strand': '+',
            'qseq_full': 'AAGAAGGGTAACCTTAAAAAATACATTTCCCACTCCAGAAAATACTCATATGTGGCCTGTTAGCAGCACAAGAAGGGTGAAAGCAATGCCCATTCCTGCCTCCCTCCCCCTGCTCACCTCCACGTCCCTGTTTGCCCCTTTACTCATATGTGGCCTGTTAGCAGCACAAGAAGGGTGAAAGCAATGCCCATTCCTGCCTCCCTCCCCCTGCTCACCTCCACGTCCCTGTTTGCCCCTTTGTAGGTGAAGTGAGTATATTCAGCGTCTTC',
            'score': 1,
        }
        read2 = Blat.pslx_row_to_pysam(pslx_row, cache, reference)
        assert read2.reference_id == 13
        assert read2.reference_start == 73014606
        assert _cigar.convert_for_igv(read2.cigar) == [(CIGAR.M, 141), (CIGAR.I, 98), (CIGAR.M, 30)]
        assert query_coverage_interval(read2) == Interval(0, len(pslx_row['qseq_full']) - 1)
