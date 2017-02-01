import unittest
from sv_pair import compare_genome_events
from structural_variant.constants import SVTYPE, COLUMNS, CALL_METHOD, ORIENT, STRAND
from structural_variant.breakpoint import BreakpointPair, Breakpoint


class TestPairing(unittest.TestCase):
    def setUp(self):
        self.g1 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.break1_call_method: CALL_METHOD.CONTIG,
                COLUMNS.break2_call_method: CALL_METHOD.CONTIG
            }
        )
        self.g2 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.break1_call_method: CALL_METHOD.CONTIG,
                COLUMNS.break2_call_method: CALL_METHOD.CONTIG
            }
        )

    def test_genome_genome_diff_chrom(self):
        self.g2.break1.chr = '2'
        self.assertFalse(compare_genome_events(self.g1, self.g2, 0, 10, 10))
        raise unittest.SkipTest('TODO')

    def test_genome_genome_diff_orient(self):
        self.g2.break1.orient = ORIENT.LEFT
        self.g1.break1.orient = ORIENT.RIGHT
        self.assertFalse(compare_genome_events(self.g1, self.g2, 0, 0, 0))

    def test_genome_genome_diff_strand(self):
        self.g2.break1.strand = STRAND.POS
        self.g1.break1.strand = STRAND.NEG
        self.assertFalse(compare_genome_events(self.g1, self.g2, 0, 0, 0))

    def test_genome_genome_diff_event_type(self):
        self.g2.data[COLUMNS.event_type] = SVTYPE.DEL
        self.g1.data[COLUMNS.event_type] = SVTYPE.INS
        self.assertFalse(compare_genome_events(self.g1, self.g2, 0, 0, 0))

    def test_genome_genome_ns_orient(self):
        self.g2.break1.orient = ORIENT.LEFT
        self.g1.break2.orient = ORIENT.RIGHT
        self.assertTrue(compare_genome_events(self.g1, self.g2, 0, 0, 0))

    def test_genome_genome_by_contig(self):
        raise unittest.SkipTest('TODO')

    def test_genome_genome_by_split(self):
        raise unittest.SkipTest('TODO')

    def test_genome_genome_by_flanking(self):
        raise unittest.SkipTest('TODO')

    def test_genome_genome_by_mixed(self):
        raise unittest.SkipTest('TODO')

    def test_genome_transcriptome(self):
        raise unittest.SkipTest('TODO')

    def test_transcriptome_transcriptome(self):
        raise unittest.SkipTest('TODO')
