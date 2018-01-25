import unittest

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CALL_METHOD, COLUMNS, PROTOCOL, STRAND, SVTYPE
from mavis.summary.summary import filter_by_annotations


class TestFilterByAnnotations(unittest.TestCase):
    def setUp(self):
        self.gev1 = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            data={
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
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
                COLUMNS.call_method: CALL_METHOD.CONTIG,
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
        result, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        bpp = result[0]
        print(bpp.data)
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
        bpps, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        print(bpps)
        bpp = bpps[0]
        print(bpp, bpp.data)
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
        result, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        bpp = result[0]
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
        result, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        bpp = result[0]
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
        result, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        bpp = result[0]
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
        result, removed = filter_by_annotations([self.gev1, self.gev2], self.best_transcripts)
        bpp = result[0]
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

    def test_get_pairing_state(self):
        raise unittest.SkipTest('TODO')
