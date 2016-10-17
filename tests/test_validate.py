from structural_variant.breakpoint import BreakpointPair, Breakpoint
from structural_variant.validate import *
from structural_variant.annotate import load_reference_genome
from structural_variant.constants import *
import structural_variant.align as align
import pysam
import unittest

HUMAN_REFERENCE_GENOME = None


class MockRead(pysam.AlignedSegment):
    """
    utility class to clean up code for generating a mock read
    for testing validate functions
    """

    def __init__(self, **kwargs):
        pysam.AlignedSegment.__init__(self)
        self.reference_start = kwargs.pop('reference_start')
        self.reference_id = kwargs.pop('reference_id')
        self.query_name = kwargs.pop(
            'query_name', 'read-{0}:{1}'.format(self.reference_id, self.reference_start))
        self.query_sequence = kwargs.pop('query_sequence')
        self.cigar = kwargs.pop('cigar')


def setUpModule():
    global HUMAN_REFERENCE_GENOME
    HUMAN_REFERENCE_GENOME = load_reference_genome('chr11_chr22.fa')


class TestEvidence(unittest.TestCase):

    def setUp(self):
        self.first_read = MockRead(
            query_sequence='TTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGGTAG'
            'CTGGAAGGTATTTAAAAGAAACTTGTTTCACCAAATTC',
            reference_id=10,
            cigar=[(CIGAR.S, 6), (CIGAR.M, 119)],
            reference_start=128664209
        )
        self.second_read = MockRead(
            query_sequence='CTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTT'
            'AGAAATTCCACCACATGGAGGTTATTTGAGGCGATTGG',
            reference_id=21,
            cigar=[(CIGAR.M, 49), (CIGAR.S, 76)],
            reference_start=29684317,
        )

        self.evidence = Evidence(
            BreakpointPair(
                Breakpoint('11', 128664209, 128664209, orient=ORIENT.RIGHT),
                Breakpoint('22', 29684365, 29684365, orient=ORIENT.LEFT),
                stranded=False, opposing_strands=False
            ),
            SVTYPE.TRANS,
            None,  # don't need bamfile b/c mocking reads
            convert_chr_to_index={'11': 10, '22': 21}
        )


if __name__ == "__main__":
    unittest.main()
