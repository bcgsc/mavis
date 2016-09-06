from structural_variant.breakpoint import BreakpointPair, Breakpoint
from structural_variant.validate import *
from structural_variant.annotate import load_reference_genome
from structural_variant.constants import *
import pysam
import unittest


class MockRead(pysam.AlignedSegment):
    """
    utility class to clean up code for generating a mock read 
    for testing validate functions
    """
    def __init__(self, **kwargs):
        pysam.AlignedSegment.__init__(self)
        self.reference_start = kwargs.pop('reference_start')
        self.reference_id = kwargs.pop('reference_id')
        self.query_name = kwargs.pop('query_name', 'read-{0}:{1}'.format(self.reference_id, self.reference_start))
        self.query_sequence = kwargs.pop('query_sequence')
        self.cigar = kwargs.pop('cigar')

def setUpModule():
    load_reference_genome('chr11_chr22.fa')

class TestEvidence(unittest.TestCase):
    def setUp(self):
        self.first_read = MockRead(
            query_sequence = 'TTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGGTAG' \
                    'CTGGAAGGTATTTAAAAGAAACTTGTTTCACCAAATTC',
            reference_id = 10,
            cigar = [(CIGAR.S, 6), (CIGAR.M, 119)],
            reference_start = 128664209
            )
        self.second_read = MockRead(
            query_sequence = 'CTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTT' \
                    'AGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGG',
            reference_id = 21,
            cigar = [(CIGAR.M, 49), (CIGAR.S, 76)],
            reference_start = 29684317,
            )
        
        self.evidence = Evidence(
                BreakpointPair(
                    Breakpoint('11', 128664209, 128664209, orient = ORIENT.RIGHT),
                    Breakpoint('22', 29684365, 29684365, orient = ORIENT.LEFT),
                    stranded = False, opposing_strands = False
                ), 
                None, # don't need bamfile b/c mocking reads 
                convert_chr_to_index={'11': 10, '22': 21}
                )

    def test_init_conversion_dict(self):
        self.assertEqual(self.evidence.convert_chr_to_index['11'], 10)
        self.assertEqual(self.evidence.convert_chr_to_index['22'], 21)
        self.assertEqual(self.evidence.convert_index_to_chr[10], '11')
        self.assertEqual(self.evidence.convert_index_to_chr[21], '22')

    def test_recompute_cigar(self):
        c = self.evidence.recompute_cigar(self.first_read)
        self.assertEqual(c, [(CIGAR.S, 6), (CIGAR.EQ, 119)])
        c = self.evidence.recompute_cigar(self.second_read)
        self.assertEqual(c, [(CIGAR.EQ, 49), (CIGAR.S, 76)])

    def test_add_split_read_second_break_shifted(self):
        self.evidence.add_split_read(self.second_read, False)
        # since this read would have been shifted, make sure the positions are correct
        # at the first breakpoint (the novel read)
        self.assertEqual(1, len(self.evidence.split_reads[self.evidence.break1]))
        temp = list(self.evidence.split_reads[self.evidence.break1])[0]
        self.assertEqual(128664209, temp.reference_start)
        self.assertEqual(10, temp.reference_id)
        self.assertEqual([(CIGAR.S, 47), (CIGAR.EQ, 78)], temp.cigar)
        self.assertEqual(78, temp.reference_end - temp.reference_start)
        # at the second breakpoint (the copied/shifted input read)
        self.assertEqual(1, len(self.evidence.split_reads[self.evidence.break2]))
        temp = list(self.evidence.split_reads[self.evidence.break2])[0]
        self.assertEqual(self.second_read.reference_start, temp.reference_start)
        self.assertEqual(self.second_read.reference_end - 2, temp.reference_end)
        self.assertEqual(21, temp.reference_id, 21)
        self.assertEqual([(CIGAR.EQ, 47), (CIGAR.S, 78)], temp.cigar)
        self.assertEqual(47, temp.reference_end - temp.reference_start)

    def test_add_split_read_first_break(self):
        self.evidence.add_split_read(self.first_read)

        # at the first breakpoint (the copied input read)
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break1]), 1)
        temp = list(self.evidence.split_reads[self.evidence.break1])[0] 
        self.assertEqual(self.first_read.reference_start, temp.reference_start)
        self.assertEqual(self.first_read.reference_end, temp.reference_end)
        self.assertEqual([(CIGAR.S, 6), (CIGAR.EQ, 119)], temp.cigar)
        self.assertEqual(119, temp.reference_end - temp.reference_start)
        # at the second breakpoint (the novel read) 
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break2]), 0) # multi-map alignment

class TestValidate(unittest.TestCase):
    """
    test class for functions in the validate namespace 
    that are not associated with a class
    """
    def test_sw_pairwise_alignment(self):
        a = sw_pairwise_alignment('ATGGACTCGGTAAA', 'CGGTAA')[0]
        self.assertEqual(a.reference_start, 7)
        self.assertEqual(a.cigar, [(CIGAR.EQ, 6)])
        self.assertEqual(a.query_sequence, 'CGGTAA')

    def test_build_string_from_reverse_path(self):
        ref = '-mxabdce'
        seq = '-abc'
        t = build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3)])
        self.assertEqual(('e', '-'), t)
        t = build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3), (5,2), (4,2)])
        self.assertEqual(('dce', '-c-'), t)
        t = build_string_from_reverse_path(ref, seq, [(6, 3), (5,2), (4,2), (3, 1), (2, 0)])
        self.assertEqual(('mxabdce', '--ab-c-'), t)

if __name__ == "__main__":
    unittest.main()
