from structural_variant.blat import Blat
from structural_variant.interval import Interval
import unittest
from Bio import SeqIO


class MockBamCache:
    def __init__(self, **kwargs):
        self.chr_to_referenceid = kwargs

    def reference_id(self, chrom):
        return self.chr_to_referenceid[chrom]


class TestBlat(unittest.TestCase):
    def setUp(self):
        self.pslx_file = 'tests/blat_output.pslx'
        self.fasta_file = 'tests/blat_input.fa'

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

    def test_pslx_row_to_pysam(self):
        pslx_row = {
            'block_count': 1,
            '_index': 0,
            'qgap_count': 0,
            'ncount': 0,
            'match': 20,
            'score': 20,
            'qstart': 93,
            'tend': 3432327,
            'qend': 113,
            'tseqs': ['AATACCAAATACATGATATA'],
            'tstarts': [3432307],
            'tgap_bases': 0,
            'tstart': 3432307,
            'mismatch': 0,
            'qgap_bases': 0,
            'block_sizes': [20],
            'qname': 'seq1',
            'tgap_count': 0,
            'tname': 'Y',
            'repmatch': 0,
            'qsize': 667,
            'qstarts': [93],
            'strand': '+',
            'qseqs': ['AATACCAAATACATGATATA'],
            'tsize': 59373566,
            'percent_ident': 100.0,
            'qseq_full': 'AGCCTCCCAAGTAGCTGGGACTACAGGCGCCCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTT'
                         'AGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCGCCCGCCTCGGC'
        }
        cache = MockBamCache(Y=23)
        read = Blat.pslx_row_to_pysam(pslx_row, cache)
        self.assertEqual(23, read.reference_id)
        self.assertEqual(Interval(93, 112), read.query_coverage_interval())
