from mavis.blat import *
from mavis.interval import Interval
from mavis.annotate import load_reference_genome
from mavis.constants import ORIENT, CIGAR
from mavis.bam.cache import BamCache
from mavis.assemble import Contig
from mavis.breakpoint import Breakpoint
from mavis.validate.evidence import GenomeEvidence
from mavis.align import query_coverage_interval

import unittest
import shutil
from Bio import SeqIO
from . import REFERENCE_GENOME_FILE, REFERENCE_GENOME_FILE_2BIT
from . import BLAT_INPUT, BLAT_OUTPUT, BAM_INPUT, MockBamFileHandle


REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(BAM_INPUT)


class TestBlat(unittest.TestCase):
    def setUp(self):
        self.cache = BamCache(MockBamFileHandle({'Y': 23, 'fake': 0, 'reference3': 3}))

    def test_read_pslx(self):
        mapping = {}
        for record in SeqIO.parse(BLAT_INPUT, 'fasta'):
            mapping[record.id] = record.seq
        header, rows = Blat.read_pslx(BLAT_OUTPUT, mapping)
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
        self.assertEqual(Interval(93, 112), query_coverage_interval(read))

    def test_pslx_row_to_pysam_full_reverse(self):
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
            'qseqs': ['TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'],
            'tseqs': ['TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'],
            '_index': 1,
            'score': 128,
            'percent_ident': 100.0,
            'qseq_full': 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, None)
        self.assertEqual(3, read.reference_id)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.M, 128)], read.cigar)
        self.assertEqual(2187, read.reference_start)
        self.assertEqual(Interval(117, 244), query_coverage_interval(read))

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
        self.assertEqual(Interval(0, 52), query_coverage_interval(read))
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
        self.assertEqual(Interval(0, 52), query_coverage_interval(read))
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
        self.assertEqual(Interval(0, 146), query_coverage_interval(read))
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
        self.assertEqual(Interval(0, 146), query_coverage_interval(read))
        self.assertEqual(950, read.reference_start)
        self.assertEqual([(CIGAR.EQ, 53), (CIGAR.D, 6236), (CIGAR.EQ, 94)], read.cigar)

    def test_pslx_row_to_pysam_revcomp_deletion(self):
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
            'qseqs': ['TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA', 'CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'],
            'tseqs': ['TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA', 'CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG']
        }
        read = Blat.pslx_row_to_pysam(pslx_row, self.cache, REFERENCE_GENOME)
        self.assertEqual(3, read.reference_id)
        self.assertEqual(Interval(0, 83), query_coverage_interval(read))
        self.assertEqual(2205, read.reference_start)
        self.assertEqual([(CIGAR.EQ, 51), (CIGAR.D, 26), (CIGAR.EQ, 33)], read.cigar)
        self.assertEqual('TAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCA', read.query_sequence[0:50])
        self.assertEqual('CCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG', read.query_sequence[50:])

    def test_pslx_row_to_pysam_inversion(self):
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
                'TTTTCATTTCTGTATGTTAAT'],
            'tseqs': [
                'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGG'
                'TTTTCATTTCTGTATGTTAAT']
        }
        read1 = Blat.pslx_row_to_pysam(pslx_row, self.cache, REFERENCE_GENOME)
        self.assertEqual(3, read1.reference_id)
        self.assertEqual(Interval(125, 244), query_coverage_interval(read1))
        self.assertEqual(1114, read1.reference_start)
        self.assertEqual([(CIGAR.S, 125), (CIGAR.EQ, 120)], read1.cigar)

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
                'TCTGTGTTTACAGGGCTTTCATGCTCAG'],
            'tseqs': [
                'TGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAAT'
                'TCTGTGTTTACAGGGCTTTCATGCTCAG']
        }
        read2 = Blat.pslx_row_to_pysam(pslx_row, self.cache, REFERENCE_GENOME)
        self.assertEqual(3, read2.reference_id)
        self.assertEqual(2187, read2.reference_start)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.EQ, 128)], read2.cigar)
        self.assertEqual(Interval(117, 244), query_coverage_interval(read2))
        self.assertEqual(read1.query_sequence, reverse_complement(read2.query_sequence))

    @unittest.skipIf(not shutil.which('blat'), "missing the blat command")
    def test_pslx_row_to_pysam_duplication(self):
        raise unittest.SkipTest('TODO')
        """
        s = 'CTCCCACCAGGAGCTCGTCCTCACCACGTCCTGCACCAGCACCTCCAGCTCCCGCAGCAGCGCCTCGCCCCCACGGTGCGC' \
            'GCTCCGCGCCGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGT' \
            'AGGTTCCATCGGCTCCGTGGGTTCCATGGACTCTGTGGGCTCGGGCCCGACGCGCACGGAGGACTGGAGGACTGGGGCGTG' \
            'TGTCTGCGGTGCAGGCGAGGCGGGGCGGGCCGGGGC'

        pslx_row = {
            'block_count': 2,
            'tstarts': [2187],
            'block_sizes': [126,10],
            'qname': 'seq1',
            'tname': 'reference17',
            'qstarts': [117],
            'strand': '+',
            'qseq_full': s,
            'score': 1
        }
        # 136	0	0	0	1	19	1	1	+	contig2	175	0	155	reference17	4000	1882	2019	2	126,10,	0,145,	1882,2009,	ctcccaccaggagctcgtcctcaccacgtcctgcaccagcacctccagctcccgcagcagcgcctcgcccccacggtgcgcgctccgcgccggttccatgggctccgtaggttccatgggctccgt,ggttccatgg,	ctcccaccaggagctcgtcctcaccacgtcctgcaccagcacctccagctcccgcagcagcgcctcgcccccacggtgcgcgctccgcgccggttccatgggctccgtaggttccatgggctccgt,ggttccatgg,
        """

    def test_ranking_alignments(self):
        pass
