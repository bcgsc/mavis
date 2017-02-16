from structural_variant.blat import *
from structural_variant.interval import Interval
from structural_variant.annotate import load_reference_genome
from structural_variant.constants import ORIENT, CIGAR
from structural_variant.read_tools import BamCache
from structural_variant.assemble import Contig
from structural_variant.validate import Evidence
from structural_variant.breakpoint import Breakpoint, BreakpointPair

import unittest
import shutil
from Bio import SeqIO
from tests import REFERENCE_GENOME_FILE, REFERENCE_GENOME_FILE_2BIT
from tests import BLAT_INPUT, BLAT_OUTPUT, BAM_INPUT, MockBamFileHandle


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
        self.cache = BamCache(MockBamFileHandle({'Y': 23, 'fake': 0,'reference3':3}))

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
        self.assertEqual(Interval(93, 112), read.query_coverage_interval())

    def test_pslx_row_to_pysam_reverse(self):
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
            'qsize': 245,
            'qstart': 0,
            'qend': 128,
            'tname': 'reference3',
            'tsize': 3711,
            'tstart': 2187,
            'tend': 2315,
            'block_count': 1,
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
        print(read.cigar)
        raise(UserWarning)

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

    @unittest.skipIf(not shutil.which('blat'), "missing the blat command")
    def test_blat_contigs(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
                Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
                opposing_strands=True
                ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_reads_resolution=1
            )
        c = Contig("CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT",0)
        ev.contigs=[c]
        blat_contigs([ev], BAM_CACHE, REFERENCE_GENOME, ref_2bit=REFERENCE_GENOME_FILE_2BIT)
        c2 = ev.contigs[0]
        print(c2.alignments)
        raise(UserWarning)
