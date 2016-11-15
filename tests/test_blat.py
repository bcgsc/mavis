from structural_variant.blat import Blat
from structural_variant.interval import Interval
from structural_variant.annotate import load_reference_genome
from structural_variant.constants import CIGAR
import unittest
from Bio import SeqIO


REFERENCE_GENOME = None


class MockBamCache:
    def __init__(self, **kwargs):
        self.chr_to_referenceid = kwargs

    def reference_id(self, chrom):
        return self.chr_to_referenceid[chrom]


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome('/home/creisle/svn/sv_compile/trunk/mock_reference_genome.fa')
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestBlat(unittest.TestCase):
    def setUp(self):
        self.pslx_file = 'tests/blat_output.pslx'
        self.fasta_file = 'tests/blat_input.fa'
        self.cache = MockBamCache(Y=23, fake=0)

    def test_read_pslx(self):
        mapping = {}
        for record in SeqIO.parse(self.fasta_file, 'fasta'):
            mapping[record.id] = record.seq
        header, rows = Blat.read_pslx(self.pslx_file, mapping)
        self.assertEqual(11067, len(rows))
        expect_pslx_header = [
            'match', 'mismatch', 'repmatch', 'ncount',
            'qgap_count', 'qgap_bases',
            'tgap_count', 'tgap_bases',
            'strand',
            'qname', 'qsize', 'qstart', 'qend',
            'tname', 'tsize', 'tstart', 'tend',
            'block_count', 'block_sizes',
            'qstarts', 'tstarts',
            'qseqs', 'tseqs'
        ]
        self.assertEqual(expect_pslx_header, header)

    def test_pslx_row_to_pysam_single_block(self):
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
                         'AGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCGCCCGCCTCGGC'
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, None)
        self.assertEqual(23, read.reference_id)
        self.assertEqual(Interval(93, 112), read.query_coverage_interval())
    
    def test_pslx_row_to_pysam_simple(self):
        pslx_row = {
            'tstarts': [950],
            'block_sizes': [53],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0],
            'strand': '+',
            'score': 0,
            'qseq_full':
                'ATCTAATAACTTGATCAATA'
                'TCTGTGATTATATTTTCATT'
                'GCCTTCCAATTTT'
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, None)
        self.assertEqual(0, read.reference_id)
        self.assertEqual(Interval(0, 52), read.query_coverage_interval())
        self.assertEqual(950, read.reference_start)
        self.assertEqual(1003, read.reference_end)
        self.assertEqual([(CIGAR.M, 53)], read.cigar)

    def test_pslx_row_to_pysam_simple_with_reference(self):
        pslx_row = {
            'tstarts': [950],
            'block_sizes': [53],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0],
            'strand': '+',
            'score': 0,
            'qseq_full':
                'ATCTAATAACTTGATCAATA'
                'TCTGTGATTATATTTTCATT'
                'GCCTTCCAATTTT'
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, REFERENCE_GENOME)
        self.assertEqual(0, read.reference_id)
        self.assertEqual(Interval(0, 52), read.query_coverage_interval())
        self.assertEqual(950, read.reference_start)
        self.assertEqual(1003, read.reference_end)
        self.assertEqual([(CIGAR.EQ, 53)], read.cigar)

    def test_pslx_row_to_pysam_gapped_alignment(self):
        pslx_row = {
            'block_count': 1,
            'tstarts': [950, 7233],
            'block_sizes': [47, 100],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0, 47],
            'strand': '+',
            'qseq_full':
                'ATCTAATAACTTGATCAATA'
                'TCTGTGATTATATTTTCATT'
                'GCCTTCC'
                'AATTTTGCAGATTATAAGAT'
                'CAATAGATATTTATTGTAAA'
                'ATGCACAAATAGTGCAACAT'
                'TTCTTAAAGTAGACCGTGAA'
                'ATACTTCATGTTGCCATGTT',
            'score': 1
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, None)
        self.assertEqual(0, read.reference_id)
        self.assertEqual(Interval(0, 146), read.query_coverage_interval())
        self.assertEqual(950, read.reference_start)
        self.assertEqual([(CIGAR.M, 47), (CIGAR.D, 6236), (CIGAR.M, 100)], read.cigar)

    def test_pslx_row_to_pysam_gapped_alignment_with_reference(self):
        pslx_row = {
            'block_count': 1,
            'tstarts': [950, 7233],
            'block_sizes': [47, 100],
            'qname': 'seq1',
            'tname': 'fake',
            'qstarts': [0, 47],
            'strand': '+',
            'qseq_full':
                'ATCTAATAACTTGATCAATA'
                'TCTGTGATTATATTTTCATT'
                'GCCTTCC'
                'AATTTTGCAGATTATAAGAT'
                'CAATAGATATTTATTGTAAA'
                'ATGCACAAATAGTGCAACAT'
                'TTCTTAAAGTAGACCGTGAA'
                'ATACTTCATGTTGCCATGTT',
            'score': 1
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, REFERENCE_GENOME)
        self.assertEqual(0, read.reference_id)
        self.assertEqual(Interval(0, 146), read.query_coverage_interval())
        self.assertEqual(950, read.reference_start)
        self.assertEqual([(CIGAR.EQ, 53), (CIGAR.D, 6236), (CIGAR.EQ, 94)], read.cigar)
