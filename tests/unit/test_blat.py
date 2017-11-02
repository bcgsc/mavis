import unittest

from mavis.blat import Blat
from mavis.constants import CIGAR, reverse_complement

from .mock import Mock, MockFunction, MockLongString


class TestConvertPslxToPysam(unittest.TestCase):

    def test_simple(self):
        row = {
            'match': 142,
            'mismatch': 0,
            'repmatch': 0,
            'ncount': 0,
            'qgap_count': 0,
            'qgap_bases': 0,
            'tgap_count': 0,
            'tgap_bases': 0,
            'strand': '-',
            'qname': 'seq1',
            'qsize': 204,
            'qstart': 0,
            'qend': 142,
            'tname': '17',
            'tsize': 81195210,
            'tstart': 32673408,
            'tend': 32673550,
            'block_count': 1,
            'block_sizes': [142],
            'qstarts': [62],
            'tstarts': [32673408],
            '_index': 880,
            'score': 142,
            'percent_ident': 100.0,
            'qseq_full': (
                'ACATGTGCACAACGTGCAGGTTTGTTACATATGTATACATGTGCCATGTTGGTTTGCTGCACCCATTAACTCGTCCTAGTTTATTACTAGTCTTCAGACATC'
                'CAGAAAATAGAGTAAGATACTAGGTAGACATAACACCTAGATACATCCGTAAGGCATTTGTTTCCTATCACATGGCCCATTCTAGCTTAACACCCACCAACT'
            )}
        refseq = {'17': Mock(seq=MockLongString(
            'ACTAGGTGTTATGTCTACCTAGTATCTTACTCTATTTTCTGGATGTCTGAAGACTAGTAATAAACTAGGACGAGTTAATGGGTGCAGCAAACCAACATGGCACATG'
            'TATACATATGTAACAAACCTGCACGTTGTGCACATGTACCCTAAAACTTAAAGTATAAAAAAAAATTTCACTGAGCATAAGACTTCAGACACAAAAGAGTGCATGC'
            'CATATAATTCCATTTATGTGAATTTCAAGAACAATCAGTGATGACAGAAGTCAAAGTAGTGGTCACCTCTGGAAGGTGGGACATTGACC',
            32673407))}
        cache = Mock(reference_id=MockFunction(16))
        read = Blat.pslx_row_to_pysam(row, cache, refseq)
        self.assertEqual(16, read.reference_id)
        self.assertEqual('17', read.reference_name)
        self.assertEqual(row['qseq_full'], reverse_complement(read.query_sequence))
        self.assertEqual([(CIGAR.S, 62), (CIGAR.EQ, 142)], read.cigar)

    def test_overlapping_blat_blocks_error(self):
        row = {
            'strand': '+',
            'qname': 'seq23',
            'tname': '7',
            'block_sizes': [54, 53, 36, 80, 29],
            'qstarts': [0, 55, 108, 143, 223],
            'tstarts': [61279112, 61279166, 61397315, 61990208, 62366144],
            'score': 207,
            'percent_ident': 91.3,
            'qseq_full': (
                'CAAAAGGAAATACCTTCACATAAATTCTAGACGGAAGCAATCTGAGAAACTTTTATTGTGATTTGTGCATTCACTTCACAGAGTTAAAACTTTCTTTTGATT'
                'GAGCAGTTTGAAACTCTGTTTTTGTAGAATCTGCAAGTGGACATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAATATCTTCACAGGAAAACTAGATA'
                'GAAGTATTCTGAGAAACTTCTTTGTGATGTATGCAGTCATATCTCAGA')
        }
        cache = Mock(reference_id=MockFunction(6))
        with self.assertRaises(AssertionError):
            Blat.pslx_row_to_pysam(row, cache, None)
