from sv import BreakpointPair, Breakpoint
from validate import Evidence
import validate
from constants import *
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

class TestEvidence(unittest.TestCase):
    def setUp(self):
        self.first_read = MockRead(
            query_sequence = 'TTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGGTAGCTGGAAGGTATTTAAAAGAAACTTGTTTCACCAAATTC',
            reference_id = 10,
            cigar = [(CIGAR.S, 6), (CIGAR.M, 119)],
            reference_start = 128664209
            )
        self.second_read = MockRead(
            query_sequence = 'CTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGG',
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
                None, 
                convert_chr_to_index={'11': 10, '22': 21}
                )
        validate.load_reference('chr11_chr22.fa')

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

    def test_add_split_read_second_break(self):
        self.evidence.add_split_read(self.second_read, False)
        # since this read would have been shifted, make sure the positions are correct
        # at the first breakpoint (the novel read)
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break1]), 1)
        temp = list(self.evidence.split_reads[self.evidence.break1])[0]
        self.assertEqual(temp.reference_start, 128664209)
        self.assertEqual(temp.reference_id, 10)
        self.assertEqual(temp.cigar, [(CIGAR.S, 47), (CIGAR.EQ, 78)])
        # at the second breakpoint (the copied/shifted input read)
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break2]), 1)
        temp = list(self.evidence.split_reads[self.evidence.break2])[0]
        self.assertEqual(temp.reference_start, 29684319)
        self.assertEqual(temp.reference_id, 21)
        self.assertEqual(temp.cigar, [(CIGAR.S, 49), (CIGAR.EQ, 76)])


    
    def test_add_split_read_first_break(self):
        self.evidence.add_split_read(self.first_read)
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break1]), 1)
        self.assertEqual(len(self.evidence.split_reads[self.evidence.break2]), 1)

class TestValidate(unittest.TestCase):
    """
    test class for functions in the validate namespace 
    that are not associated with a class
    """
    def test_sw_pairwise_alignment(self):
        a = validate.sw_pairwise_alignment('ATGGACTCGGTAAA', 'CGGTAA')
        self.assertEqual(a.reference_start, 7)
        self.assertEqual(a.cigar, [(CIGAR.EQ, 6)])
        self.assertEqual(a.query_sequence, 'CGGTAA')

    def test__build_string_from_reverse_path(self):
        ref = '-mxabdce'
        seq = '-abc'
        t = validate._build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3)])
        self.assertEqual(('e', '-'), t)
        t = validate._build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3), (5,2), (4,2)])
        self.assertEqual(('dce', '-c-'), t)
        t = validate._build_string_from_reverse_path(ref, seq, [(6, 3), (5,2), (4,2), (3, 1), (2, 0)])
        self.assertEqual(('mxabdce', '--ab-c-'), t)
#f = '/projects/seqref/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'
#f = 'chr11_chr22.fa'
#print('loading the human reference genome', f)
#validate.load_reference(f)
#print('finished loading:', f)
#
#bp = BreakpointPair(
#            Breakpoint('11', 128664209, 128664209, orient = ORIENT.RIGHT),
#            Breakpoint('22', 29684365, 29684365, orient = ORIENT.LEFT),
#            stranded = False, opposing_strands = False)
#
#bf = '/projects/analysis/analysis14/P00159/merge_bwa/125nt/hg19a/P00159_5_lanes_dupsFlagged.bam'
#
#split_read = MockRead(
#        query_sequence = 'TTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGGTAGCTGG' \
#                'AAGGTATTTAAAAGAAACTTGTTTCACCAAATTC',
#        reference_id = 10,
#        cigar = [(CIGAR.S, 6), (CIGAR.M, 119)],
#        reference_start = 128664209
#        )
#print(dir(split_read))
#e = Evidence(bp, bf, convert_chr_to_index={'11': 10, '22': 21})
#
#print(e.convert_chr_to_index, e.convert_index_to_chr)
#
#print(split_read.cigar)
#print(e.recompute_cigar(split_read))
#e.add_split_read(split_read)
#
#split_read = MockRead(
#        query_sequence = 'CTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTA' \
#            'GAAATTCCACCACATGGAGGTTATTTGAGGCGATTGG',
#        reference_id = 21,
#        cigar = [(CIGAR.M, 49), (CIGAR.S, 76)],
#        reference_start = 29684317,
#        )
#
#e.add_split_read(split_read, False)
#
#exit()
#
##e.load_evidence()
##e.resolve_breakpoint(e.break1)
#for breakpoint in e.split_reads:
#    print(breakpoint)
#    for read in e.split_reads[breakpoint]:
#        print('split read evidence')
#        print(read)
#        print(validate.str_cigar(read))
#        print(read.cigar)
#        ref = validate.HUMAN_REFERENCE_GENOME[breakpoint.chr].seq[read.reference_start:read.reference_end]
#        print(ref)
#        print(read.query_sequence[read.query_alignment_start:read.query_alignment_end])
#        print()

if __name__ == "__main__":
    unittest.main()
