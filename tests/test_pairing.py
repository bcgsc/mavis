import unittest
from structural_variant.pairing import *
from structural_variant.constants import SVTYPE, COLUMNS, CALL_METHOD, ORIENT, STRAND, PROTOCOL
from structural_variant.breakpoint import BreakpointPair, Breakpoint
from structural_variant.annotate.genomic import usTranscript, Exon


class TestPairing(unittest.TestCase):
    def setUp(self):
        self.gev1 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.break1_call_method: CALL_METHOD.CONTIG,
                COLUMNS.break2_call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME
            }
        )
        self.gev2 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.break1_call_method: CALL_METHOD.CONTIG,
                COLUMNS.break2_call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME
            }
        )

        self.ust1 = usTranscript(
            exons=[(1, 100), (301, 400), (501, 600)],
            strand=STRAND.POS,
            name='t1'
        )
        self.ust2 = usTranscript(
            exons=[(1001, 1100), (1301, 1400), (1501, 1600)],
            strand=STRAND.POS,
            name='t2'
        )
        self.DISTANCES = {CALL_METHOD.CONTIG: 0, CALL_METHOD.FLANK: 0, CALL_METHOD.SPLIT: 10}
        self.TRANSCRIPTS = {
            self.ust1.name: self.ust1,
            self.ust2.name: self.ust2
        }

    def test_genome_protocol_diff_chrom(self):
        self.gev2.break1.chr = '2'
        self.assertFalse(equivalent_events(self.gev1, self.gev2, self.TRANSCRIPTS))
        raise unittest.SkipTest('TODO')

    def test_genome_protocol_diff_orient(self):
        self.gev2.break1.orient = ORIENT.LEFT
        self.gev1.break1.orient = ORIENT.RIGHT
        self.assertFalse(equivalent_events(self.gev1, self.gev2, self.TRANSCRIPTS))

    def test_genome_protocol_diff_strand(self):
        self.gev2.break1.strand = STRAND.POS
        self.gev1.break1.strand = STRAND.NEG
        self.assertFalse(equivalent_events(self.gev1, self.gev2, self.TRANSCRIPTS))

    def test_genome_protocol_diff_event_type(self):
        self.gev2.data[COLUMNS.event_type] = SVTYPE.DEL
        self.gev1.data[COLUMNS.event_type] = SVTYPE.INS
        self.assertFalse(equivalent_events(self.gev1, self.gev2, self.TRANSCRIPTS))

    def test_genome_protocol_ns_orient(self):
        self.gev2.break1.orient = ORIENT.LEFT
        self.gev1.break2.orient = ORIENT.RIGHT
        self.assertTrue(equivalent_events(self.gev1, self.gev2, self.TRANSCRIPTS))

    def test_genome_protocol_by_contig(self):
        raise unittest.SkipTest('TODO')

    def test_genome_protocol_by_split(self):
        raise unittest.SkipTest('TODO')

    def test_genome_protocol_by_flanking(self):
        raise unittest.SkipTest('TODO')

    def test_genome_protocol_by_mixed_call_methods(self):
        raise unittest.SkipTest('TODO')

    def test_mixed_protocol(self):
        raise unittest.SkipTest('TODO')

    def test_transcriptome_protocol(self):
        raise unittest.SkipTest('TODO')


class TestBreakpointPrediction(unittest.TestCase):
    def setUp(self):
        self.ust = usTranscript([Exon(101, 200), Exon(301, 400), Exon(501, 600)], strand=STRAND.POS)

    def test_exonic_five_prime(self):
        b = Breakpoint('1', 350, orient=ORIENT.LEFT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(2, len(breaks))
        self.assertEqual(200, breaks[0].start)
        self.assertEqual(b, breaks[1])

    def test_exonic_five_prime_first_exon(self):
        b = Breakpoint('1', 150, orient=ORIENT.LEFT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(1, len(breaks))
        self.assertEqual(b, breaks[0])

    def test_exonic_three_prime(self):
        b = Breakpoint('1', 350, orient=ORIENT.RIGHT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(2, len(breaks))
        self.assertEqual(501, breaks[1].start)
        self.assertEqual(b, breaks[0])

    def test_exonic_three_prime_last_exon(self):
        b = Breakpoint('1', 550, orient=ORIENT.RIGHT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(1, len(breaks))
        self.assertEqual(b, breaks[0])

    def test_intronic_five_prime(self):
        b = Breakpoint('1', 450, orient=ORIENT.LEFT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(1, len(breaks))
        self.assertEqual(400, breaks[0].start)

    def test_intronic_three_prime(self):
        b = Breakpoint('1', 250, orient=ORIENT.RIGHT)
        breaks = predict_transcriptome_breakpoint(b, self.ust)
        self.assertEqual(1, len(breaks))
        self.assertEqual(301, breaks[0].start)
