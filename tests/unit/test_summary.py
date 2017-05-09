import unittest

from mavis.summary.summary import alphanumeric_choice
from mavis.breakpoint import BreakpointPair, Breakpoint
from mavis.constants import SVTYPE, COLUMNS, CALL_METHOD, ORIENT, STRAND, PROTOCOL
from mavis.annotate.genomic import usTranscript, Transcript, Exon, Gene


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


class TestCompareBppAnnotations(unittest.TestCase):
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

        self.annotations = {}
        gene = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        gene2 = Gene('1', 10000, 12000, name='ABCA', strand=STRAND.POS)
        self.a_ust = usTranscript(name='ENST00000248690', gene=gene2, exons=[(10101, 10200), (10301, 10400), (10501, 10600)])
        self.c_ust = usTranscript(name='ENST00000248797', gene=gene2, exons=[(10101, 10200), (10301, 10400), (10501, 10600), (11050, 11500)], is_best_transcript=True)

        self.ust = usTranscript(name='ENST00000367580', gene=gene, exons=[(1001, 1100), (1401, 1500), (3001, 3999)])
        self.b_ust = usTranscript(name='ENST00000367579', gene=gene, exons=[(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)], is_best_transcript=True)
        gene.unspliced_transcripts.append(self.ust)
        gene.unspliced_transcripts.append(self.b_ust)
        gene2.unspliced_transcripts.append(self.a_ust)
        gene2.unspliced_transcripts.append(self.c_ust)
        for trans in (self.a_ust, self.ust, self.b_ust, self.c_ust):
            for spl in trans.generate_splicing_patterns():
                t = Transcript(trans, spl)
                trans.transcripts.append(t)
        # for spl in self.a_ust.generate_splicing_patterns():
        #     t = Transcript(self.a_ust, spl)
        #     self.a_ust.transcripts.append(t)
        # for spl in self.c_ust.generate_splicing_patterns():
        #     t = Transcript(self.c_ust, spl)
        #     self.c_ust.transcripts.append(t)
        # for spl in self.b_ust.generate_splicing_patterns():
        #     t = Transcript(self.b_ust, spl)
        #     self.b_ust.transcripts.append(t)
        self.annotations[gene.chr] = [gene, gene2]


    def test_compare_bbp_annotations_two_best_transcripts(self):
        print(self.annotations['1'][0].transcripts[1].exons)
        self.assertTrue(False)

        self.ust = usTranscript([Exon(101, 200), Exon(301, 400), Exon(501, 600)], strand=STRAND.POS)
        raise unittest.SkipTest('TODO')

    # def test_compare_bpp_annotations_two_transcripts(self):
    #     raise unittest.SkipTest('TODO')

    # def test_compare_bbp_annotations_two_fusion_cdna(self):
    #     raise unittest.SkipTest('TODO')

    # def test_compare_bbp_annotations_one_transcripts(self):
    #     raise unittest.SkipTest('TODO')

    # def test_compare_bbp_annotations_one_best_transcripts(self):
    #     raise unittest.SkipTest('TODO')

    # def test_compare_bbp_annotations_no_transcripts(self):
    #     raise unittest.SkipTest('TODO')

    # def test_aggregate_events(self):
    #     raise unittest.SkipTest('TODO')

    # def test_filtering_events_contigs(self):
    #     raise unittest.SkipTest('TODO')

    # def test_filtering_events_none(self):
    #     raise unittest.SkipTest('TODO')

    # def test_filtering_events_flanking(self):
    #     raise unittest.SkipTest('TODO')

    # def test_filtering_events_spanning(self):
    #     raise unittest.SkipTest('TODO')

    # def test_filtering_events_split(self):
    #     raise unittest.SkipTest('TODO')
