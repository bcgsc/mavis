import unittest

from mavis.summary.summary import alphanumeric_choice, filter_by_annotations
from mavis.breakpoint import BreakpointPair, Breakpoint
from mavis.constants import SVTYPE, COLUMNS, CALL_METHOD, STRAND, PROTOCOL


class TestSummary(unittest.TestCase):
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

    def test_alphanumeric_choice(self):
        self.gev1.data[COLUMNS.transcript1] = 'ABC'
        self.gev1.data[COLUMNS.transcript2] = 'AB1'
        self.gev2.data[COLUMNS.transcript1] = 'ZED'
        self.gev2.data[COLUMNS.transcript2] = 'AB1'
        bpp = alphanumeric_choice(self.gev1, self.gev2)
        self.assertEqual('ABC', bpp.data[COLUMNS.transcript1])

    def test_alphanumeric_choice_numbers(self):
        self.gev1.data[COLUMNS.transcript1] = '123'
        self.gev1.data[COLUMNS.transcript2] = '345'
        self.gev2.data[COLUMNS.transcript1] = '567'
        self.gev2.data[COLUMNS.transcript2] = '890'
        bpp = alphanumeric_choice(self.gev1, self.gev2)
        self.assertEqual('123', bpp.data[COLUMNS.transcript1])

    def test_alphanumeric_choice_gene_names(self):
        self.gev1.data[COLUMNS.transcript1] = 'ENST00000367580'
        self.gev1.data[COLUMNS.transcript2] = 'ENST00000367580'
        self.gev2.data[COLUMNS.transcript1] = 'ENST00000367579'
        self.gev2.data[COLUMNS.transcript2] = 'ENST00000367579'
        bpp = alphanumeric_choice(self.gev1, self.gev2)
        self.assertEqual('ENST00000367579', bpp.data[COLUMNS.transcript1])


class TestFilterByAnnotations(unittest.TestCase):
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
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.fusion_cdna_coding_end: None,
                COLUMNS.fusion_cdna_coding_start: None
            }
        )
        self.gev2 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 100),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.break1_call_method: CALL_METHOD.CONTIG,
                COLUMNS.break2_call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.fusion_cdna_coding_start: None,
                COLUMNS.fusion_cdna_coding_end: None
            }
        )
        self.best_transcripts = {'ABCA': True, 'ABCD': True}

    def test_filter_by_annotations_two_best_transcripts(self):
        self.gev1.data[COLUMNS.gene1] = 'ABC'
        self.gev1.data[COLUMNS.gene2] = 'ABC'
        self.gev1.data[COLUMNS.transcript1] = 'ABCA'
        self.gev1.data[COLUMNS.transcript2] = 'ABCA'
        self.gev2.data[COLUMNS.gene1] = 'ABC'
        self.gev2.data[COLUMNS.gene2] = 'ABC'
        self.gev2.data[COLUMNS.transcript1] = 'ABCD'
        self.gev2.data[COLUMNS.transcript2] = 'ABCD'
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev1, bpp)
        self.assertEqual('ABCA', bpp.data[COLUMNS.transcript1])

    def test_filter_by_annotations_two_transcripts(self):
        self.gev1.data[COLUMNS.gene1] = 'XYZ'
        self.gev1.data[COLUMNS.gene2] = 'XYS'
        self.gev1.data[COLUMNS.transcript1] = 'XYZB'
        self.gev1.data[COLUMNS.transcript2] = 'XYSZ'
        self.gev2.data[COLUMNS.gene1] = 'XYZ'
        self.gev2.data[COLUMNS.gene2] = 'XYS'
        self.gev2.data[COLUMNS.transcript1] = 'XYZA'
        self.gev2.data[COLUMNS.transcript2] = 'XYSB'
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev2, bpp)
        self.assertEqual('XYZA', bpp.data[COLUMNS.transcript1])

    def test_filter_by_annotations_two_fusion_cdna(self):
        self.gev1.data[COLUMNS.gene1] = 'XYZ'
        self.gev1.data[COLUMNS.gene2] = 'XYS'
        self.gev1.data[COLUMNS.transcript1] = 'XYZB'
        self.gev1.data[COLUMNS.transcript2] = 'XYSZ'
        self.gev2.data[COLUMNS.gene1] = 'XYZ'
        self.gev2.data[COLUMNS.gene2] = 'XYS'
        self.gev2.data[COLUMNS.transcript1] = 'XYZB'
        self.gev2.data[COLUMNS.transcript2] = 'XYSZ'
        self.gev1.data[COLUMNS.fusion_cdna_coding_start] = 1
        self.gev1.data[COLUMNS.fusion_cdna_coding_end] = 20
        self.gev2.data[COLUMNS.fusion_cdna_coding_start] = 1
        self.gev2.data[COLUMNS.fusion_cdna_coding_end] = 40
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev2, bpp)

    def test_filter_by_annotations_one_transcript(self):
        self.gev1.data[COLUMNS.gene1] = None
        self.gev1.data[COLUMNS.gene2] = 'XYS'
        self.gev1.data[COLUMNS.transcript1] = None
        self.gev1.data[COLUMNS.transcript2] = 'XYSZ'
        self.gev2.data[COLUMNS.gene1] = 'XYZ'
        self.gev2.data[COLUMNS.gene2] = 'XYS'
        self.gev2.data[COLUMNS.transcript1] = 'XYZA'
        self.gev2.data[COLUMNS.transcript2] = 'XYSB'
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev2, bpp)

    def test_filter_by_annotations_one_best_transcripts(self):
        self.gev1.data[COLUMNS.gene1] = 'XYZ'
        self.gev1.data[COLUMNS.gene2] = 'ABC'
        self.gev1.data[COLUMNS.transcript1] = 'XYZB'
        self.gev1.data[COLUMNS.transcript2] = 'ABCA'
        self.gev2.data[COLUMNS.gene1] = 'XYZ'
        self.gev2.data[COLUMNS.gene2] = 'ABC'
        self.gev2.data[COLUMNS.transcript1] = 'XYZA'
        self.gev2.data[COLUMNS.transcript2] = 'ABCB'
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev1, bpp)
        self.assertEqual('XYZB', bpp.data[COLUMNS.transcript1])

    def test_filter_by_annotations_no_transcripts(self):
        self.gev1.data[COLUMNS.gene1] = None
        self.gev1.data[COLUMNS.gene2] = None
        self.gev1.data[COLUMNS.transcript1] = None
        self.gev1.data[COLUMNS.transcript2] = None
        self.gev2.data[COLUMNS.gene1] = None
        self.gev2.data[COLUMNS.gene2] = None
        self.gev2.data[COLUMNS.transcript1] = None
        self.gev2.data[COLUMNS.transcript2] = None
        self.gev1.break1.strand = STRAND.POS
        bpp = filter_by_annotations(self.gev1, self.gev2, self.best_transcripts)
        self.assertEqual(self.gev1, bpp)
        self.assertEqual(None, bpp.data[COLUMNS.transcript1])

    def test_combine_events(self):
        raise unittest.SkipTest('TODO')

    def test_filtering_events_contigs(self):
        raise unittest.SkipTest('TODO')

    def test_filtering_events_none(self):
        raise unittest.SkipTest('TODO')

    def test_filtering_events_flanking(self):
        raise unittest.SkipTest('TODO')

    def test_filtering_events_spanning(self):
        raise unittest.SkipTest('TODO')

    def test_filtering_events_split(self):
        raise unittest.SkipTest('TODO')
