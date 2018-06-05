import os
import unittest

from mavis.annotate.base import BioInterval, ReferenceName
from mavis.annotate.file_io import load_reference_genes, load_reference_genome
from mavis.annotate.genomic import Exon, Gene, Template, Transcript, PreTranscript
from mavis.annotate.protein import calculate_orf, Domain, DomainRegion, translate, Translation
from mavis.annotate.variant import _gather_annotations, _gather_breakpoint_annotations, annotate_events, Annotation, flatten_fusion_transcript, overlapping_transcripts
from mavis.annotate.fusion import determine_prime, FusionTranscript
from mavis.annotate.constants import SPLICE_TYPE
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, PRIME, PROTOCOL, reverse_complement, STRAND, SVTYPE
from mavis.error import NotSpecifiedError
from mavis.interval import Interval

from . import MockLongString, MockObject, get_example_genes
from ..util import get_data


REFERENCE_ANNOTATIONS = None
REFERENCE_GENOME = None
REF_CHR = 'fake'
ALT_REF_CHR = 'ref2'


def setUpModule():
    global REFERENCE_ANNOTATIONS, REFERENCE_GENOME, REF_CHR, EXAMPLE_GENES
    EXAMPLE_GENES = get_example_genes()
    REFERENCE_ANNOTATIONS = load_reference_genes(get_data('mock_reference_annotations.tsv'))
    count = sum([len(l) for l in REFERENCE_ANNOTATIONS.values()])
    print('loaded annotations', count)
    assert(count >= 6)  # make sure this is the file we expect
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    assert(REF_CHR in REFERENCE_GENOME)
    print('loaded the reference genome', get_data('mock_reference_genome.fa'))


class TestTemplate(unittest.TestCase):

    def test_template_hashing(self):
        t = Template('1', 1, 10)
        d = {'1': 1, '2': 2, 1: '5'}
        self.assertEqual('1', t.name)
        self.assertEqual(1, d[t.name])
        self.assertEqual(1, d[t])


class TestFusionTranscript(unittest.TestCase):

    def setUp(self):
        self.x = Interval(100, 199)  # C
        self.y = Interval(500, 599)  # G
        self.z = Interval(1200, 1299)  # T
        self.w = Interval(1500, 1599)  # C
        self.s = Interval(1700, 1799)  # G
        # introns: 99, 300, 600, 200, 100, ...
        reference_sequence = 'A' * 99 + 'C' * 100 + 'A' * 300 + 'G' * 100
        reference_sequence += 'A' * 600 + 'T' * 100 + 'A' * 200 + 'C' * 100
        reference_sequence += 'A' * 100 + 'G' * 100 + 'A' * 200 + 'T' * 100

        self.a = Interval(2000, 2099)  # T
        self.b = Interval(2600, 2699)  # C
        self.c = Interval(3000, 3099)  # G
        self.d = Interval(3300, 3399)  # T
        reference_sequence += 'A' * 500 + 'C' * 100 + 'A' * 300 + 'G' * 100
        reference_sequence += 'A' * 200 + 'T' * 100 + 'A' * 200
        self.reference_sequence = reference_sequence

        self.b1 = Interval(600, 699)  # A
        self.b2 = Interval(800, 899)  # G
        self.b3 = Interval(1100, 1199)  # T
        self.b4 = Interval(1400, 1499)  # A
        self.b5 = Interval(1700, 1799)  # G
        self.b6 = Interval(2100, 2199)  # A
        alternate_sequence = 'C' * 599 + 'A' * 100 + 'C' * 100 + 'G' * 100
        alternate_sequence += 'C' * 200 + 'T' * 100 + 'C' * 200 + 'A' * 100
        alternate_sequence += 'C' * 200 + 'G' * 100 + 'C' * 300 + 'A' * 100
        alternate_sequence += 'C' * 200
        self.alternate_sequence = alternate_sequence

    def test__pull_exons_left_pos_intronic(self):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 700, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (700 - 600 + 1)
        self.assertEqual(expt, seq)
        self.assertEqual(2, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)

    def test__pull_exons_left_pos_intronic_splice(self):
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 201, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 100 + 'A' * 2
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(False, e.end_splice_site.intact)

    def test__pull_exons_left_pos_exonic(self):
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        print('transcriptt exons:', t.exons)
        b = Breakpoint(REF_CHR, 199, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(False, e.end_splice_site.intact)

    def test__pull_exons_left_pos_exonic_splice(self):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 101, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 2
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(2, e.end)
        self.assertEqual(False, e.start_splice_site.intact)
        self.assertEqual(False, e.end_splice_site.intact)

    def test__pull_exons_right_pos_intronic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 1600, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))

        b = Breakpoint(REF_CHR, 300, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (499 - 300 + 1) + 'G' * 100 + 'A' * (1199 - 600 + 1) + 'T' * 100
        expt += 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100

        self.assertEqual(expt, seq)
        self.assertEqual(4, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)

    def test__pull_exons_right_pos_intronic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'AA' + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)

    def test__pull_exons_right_pos_exonic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 1201, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'T' * 99 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)

    def test__pull_exons_right_pos_exonic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint(REF_CHR, 1298, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'TT' + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.start_splice_site.intact)
        self.assertEqual(False, e.end_splice_site.intact)

    def test__pull_exons_right_neg_intronic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b = Breakpoint(REF_CHR, 700, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (1199 - 700 + 1) + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100
        expt += 'A' * (1699 - 1600 + 1) + 'G' * 100
        expt = reverse_complement(expt)
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual('C' * 100, seq[e.start - 1:e.end])
        e = new_exons[1][0]
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual('G' * 100, seq[e.start - 1:e.end])

    def test__pull_exons_right_neg_intronic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b = Breakpoint(REF_CHR, 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'AA' + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        expt = reverse_complement(expt)
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual('C' * 100, seq[e.start - 1:e.end])
        e = new_exons[1][0]
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(True, e.end_splice_site.intact)
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual('G' * 100, seq[e.start - 1:e.end])
        e = new_exons[2][0]
        self.assertEqual(True, e.start_splice_site.intact)
        self.assertEqual(False, e.end_splice_site.intact)
        self.assertEqual(501, e.start)
        self.assertEqual(600, e.end)
        self.assertEqual('A' * 100, seq[e.start - 1:e.end])

    def test_build_single_transcript_indel(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 599, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)

        self.assertEqual(expt, ft.seq)
        self.assertEqual(5, len(ft.exons))

        for i, ex in enumerate(t.exons):
            n = ft.exons[i]
            self.assertEqual(ex, ft.exon_mapping[n.position])

        self.assertEqual(1, ft.exons[0].start)
        self.assertEqual(100, ft.exons[0].end)

        splice_pattern = [(True, True), (True, False), (False, True), (True, True), (True, True)]
        char_pattern = [x * 100 for x in ['C', 'G', 'T', 'C', 'G']]

        for i in range(0, len(splice_pattern)):
            s, t = splice_pattern[i]
            ex = ft.exons[i]
            self.assertEqual(s, ex.start_splice_site.intact)
            self.assertEqual(t, ex.end_splice_site.intact)
            self.assertEqual(char_pattern[i], ft.seq[ex.start - 1:ex.end])

    def test_build_single_transcript_inversion(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)
        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'ATCGTC' + 'A' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        exons = [(1, 100), (401, 500), (1407, 1506), (1607, 1706)]
        for i in range(len(exons)):
            self.assertEqual(exons[i][0], ft.exons[i].start)
            self.assertEqual(exons[i][1], ft.exons[i].end)
        self.assertEqual(expt, ft.seq)
        self.assertEqual(4, len(ft.exons))

    def test_build_single_transcript_inversion_transcriptome(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.TRANS)
        ft = FusionTranscript.build(ann, ref)
        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'ATCGTC' + 'A' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        exons = [
            Exon(1, 100, strand=STRAND.POS),
            Exon(401, 500, intact_end_splice=False, strand=STRAND.POS),
            Exon(501, 1406, intact_start_splice=False, intact_end_splice=False, strand=STRAND.POS),
            Exon(1407, 1506, intact_start_splice=False, strand=STRAND.POS),
            Exon(1607, 1706, strand=STRAND.POS)
        ]
        print(ft.exons)
        for i in range(len(exons)):
            self.assertEqual(exons[i].start, ft.exons[i].start)
            self.assertEqual(exons[i].end, ft.exons[i].end)
            self.assertEqual(exons[i].start_splice_site.intact, ft.exons[i].start_splice_site.intact)
            self.assertEqual(exons[i].end_splice_site.intact, ft.exons[i].end_splice_site.intact)
        self.assertEqual(expt, ft.seq)
        self.assertEqual(5, len(ft.exons))

    def test_build_single_transcript_inversion_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b1 = Breakpoint(REF_CHR, 1300, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.s) + 'T' * (1699 - 1600 + 1) + 'G' * len(self.w) + 'T' * (1499 - 1300 + 1)
        expt += 'T' * len(self.z) + 'GACGAT' + 'T' * (1199 - 600 + 1) + 'C' * len(self.y)
        expt += 'T' * (499 - 200 + 1) + 'G' * len(self.x)

        exons = [(1, 100), (201, 300), (1207, 1306), (1607, 1706)]

        for i in range(len(exons)):
            self.assertEqual(exons[i][0], ft.exons[i].start)
            self.assertEqual(exons[i][1], ft.exons[i].end)
        self.assertEqual(expt, ft.seq)
        self.assertEqual(4, len(ft.exons))

    def test_build_single_transcript_duplication_pos(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)
        self.assertEqual(STRAND.POS, ft.get_strand())

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'T' * len(self.z) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        self.assertEqual(expt, ft.seq)
        exons = [(1, 100), (401, 500), (1101, 1200), (1209, 1308), (1509, 1608), (1709, 1808)]
        for i in range(len(exons)):
            self.assertEqual(exons[i][0], ft.exons[i].start)
            self.assertEqual(exons[i][1], ft.exons[i].end)

        self.assertEqual(6, len(ft.exons))
        self.assertTrue(ft.exons[2].start_splice_site.intact)
        self.assertTrue(ft.exons[3].end_splice_site.intact)
        self.assertFalse(ft.exons[2].end_splice_site.intact)
        self.assertFalse(ft.exons[3].start_splice_site.intact)

    def test_build_single_transcript_duplication_pos_transcriptome(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.TRANS)
        ft = FusionTranscript.build(ann, ref)
        self.assertEqual(STRAND.POS, ft.get_strand())

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'T' * len(self.z) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        self.assertEqual(expt, ft.seq)
        exons = [
            Exon(1, 100, strand=STRAND.POS),
            Exon(401, 500, strand=STRAND.POS),
            Exon(1101, 1200, intact_end_splice=False, strand=STRAND.POS),
            Exon(1201, 1208, intact_start_splice=False, intact_end_splice=False, strand=STRAND.POS),
            Exon(1209, 1308, intact_start_splice=False, strand=STRAND.POS),
            Exon(1509, 1608, strand=STRAND.POS),
            Exon(1709, 1808, strand=STRAND.POS)]
        print(ft.exons)
        for i in range(len(exons)):
            self.assertEqual(exons[i].start, ft.exons[i].start)
            self.assertEqual(exons[i].end, ft.exons[i].end)
            self.assertEqual(exons[i].start_splice_site.intact, ft.exons[i].start_splice_site.intact)
            self.assertEqual(exons[i].end_splice_site.intact, ft.exons[i].end_splice_site.intact)

        self.assertEqual(7, len(ft.exons))

    def test_build_single_transcript_duplication_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'T' * len(self.z) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        expt = reverse_complement(expt)
        self.assertEqual(expt, ft.seq)

        exons = [(1, 100), (201, 300), (501, 600), (609, 708), (1309, 1408), (1709, 1808)]

        for i in range(len(exons)):
            self.assertEqual(exons[i][0], ft.exons[i].start)
            self.assertEqual(exons[i][1], ft.exons[i].end)

        self.assertEqual(6, len(ft.exons))
        self.assertTrue(ft.exons[2].start_splice_site.intact)
        self.assertTrue(ft.exons[3].end_splice_site.intact)
        self.assertFalse(ft.exons[2].end_splice_site.intact)
        self.assertFalse(ft.exons[3].start_splice_site.intact)
        self.assertEqual(3, ft.exon_number(ft.exons[2]))
        self.assertEqual(3, ft.exon_number(ft.exons[3]))

    def test_build_two_transcript_inversion_5prime_pos(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.NEG)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGACTC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)
        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'ATCGACTC' + 'G' * len(self.b) + 'T' * (2599 - 2100 + 1) + 'A' * len(self.a)
        self.assertEqual(expt, ft.seq)
        self.assertEqual(4, len(ft.exons))
        self.assertTrue(ft.exons[3].end_splice_site.intact)
        self.assertFalse(ft.exons[2].start_splice_site.intact)
        self.assertTrue(ft.exons[2].end_splice_site.intact)
        self.assertEqual(2, ft.exon_number(ft.exons[1]))
        self.assertEqual(3, ft.exon_number(ft.exons[2]))

    def test_build_two_transcript_inversion_5prime_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGACTC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)
        expt = 'T' * len(self.a) + 'A' * (2599 - 2100 + 1) + 'C' * len(self.b) + 'ATCGACTC'
        expt += 'T' * (1199 - 600 + 1) + 'C' * len(self.y) + 'T' * (499 - 200 + 1) + 'G' * len(self.x)

        self.assertEqual(4, len(ft.exons))
        self.assertEqual(2, ft.exon_number(ft.exons[1]))
        self.assertEqual(4, ft.exon_number(ft.exons[2]))
        self.assertEqual(expt, ft.seq)

    def test_build_two_transcript_duplication_pos(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGAC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)
        expt = 'T' * len(self.a) + 'A' * (2599 - 2100 + 1) + 'C' * len(self.b) + 'ATCGAC'
        expt += 'T' * len(self.z) + 'A' * (1499 - 1300 + 1) + 'C' * len(self.w)
        expt += 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)

        self.assertEqual(5, len(ft.exons))
        self.assertEqual(2, ft.exon_number(ft.exons[1]))
        self.assertEqual(3, ft.exon_number(ft.exons[2]))
        self.assertEqual(expt, ft.seq)

    def test_build_two_transcript_duplication_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.NEG)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGAC')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.s) + 'T' * (1699 - 1600 + 1) + 'G' * len(self.w) + 'T' * (1499 - 1300 + 1)
        expt += 'A' * len(self.z) + 'GTCGAT' + 'G' * len(self.b) + 'T' * (2599 - 2100 + 1)
        expt += 'A' * len(self.a)

        self.assertEqual(5, len(ft.exons))
        self.assertEqual(2, ft.exon_number(ft.exons[1]))
        self.assertEqual(3, ft.exon_number(ft.exons[2]))
        self.assertEqual(expt, ft.seq)

    def test_build_two_transcript_deletion_pos(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2700, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGTGT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1) + 'AACGTGT'
        expt += 'A' * (2999 - 2700 + 1) + 'G' * len(self.c) + 'A' * (3299 - 3100 + 1) + 'T' * len(self.d)

        self.assertEqual(expt, ft.seq)
        self.assertTrue(4, len(ft.exons))

    def test_build_two_transcript_deletion_pos_transcriptome(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2700, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGTGT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.TRANS)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1) + 'AACGTGT'
        expt += 'A' * (2999 - 2700 + 1) + 'G' * len(self.c) + 'A' * (3299 - 3100 + 1) + 'T' * len(self.d)

        self.assertEqual(expt, ft.seq)
        self.assertTrue(5, len(ft.exons))

    def test_build_two_transcript_deletion_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)

        t2 = PreTranscript(exons=[self.a, self.b, self.c, self.d], strand=STRAND.NEG)
        print('t1 exons', t1.exons)
        print('t2 exons', t2.exons)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGAGTGT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME)

        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.s) + 'T' * (1699 - 1600 + 1) + 'G' * len(self.w) + 'T' * (1499 - 1300 + 1)
        expt += 'A' * len(self.z) + 'ACACTCGTT' + 'G' * len(self.b) + 'T' * (2599 - 2100 + 1)
        expt += 'A' * len(self.a)

        self.assertEqual(expt, ft.seq)
        self.assertTrue(5, len(ft.exons))
        self.assertEqual(3, ft.exon_number(ft.exons[2]))
        self.assertEqual(3, ft.exon_number(ft.exons[3]))

    def test_build_two_transcript_translocation(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        t2 = PreTranscript(exons=[self.b1, self.b2, self.b3, self.b4, self.b5, self.b6], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint('ref2', 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='GCAACATAT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence), 'ref2': MockObject(seq=self.alternate_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME)
        self.assertEqual(b1, ann.break1)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'GCAACATAT' + 'C' * (1399 - 1200 + 1) + 'A' * len(self.b4) + 'C' * (1699 - 1500 + 1)
        expt += 'G' * len(self.b5) + 'C' * (2099 - 1800 + 1) + 'A' * len(self.b6)

        self.assertEqual(expt, ft.seq)
        self.assertTrue(5, len(ft.exons))
        self.assertTrue(2, ft.exon_number(ft.exons[1]))
        self.assertTrue(4, ft.exon_number(ft.exons[2]))

    def test_build_two_transcript_translocation_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        t2 = PreTranscript(exons=[self.b1, self.b2, self.b3, self.b4, self.b5, self.b6], strand=STRAND.NEG)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(ALT_REF_CHR, 1199, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='TCTACATAT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence), ALT_REF_CHR: MockObject(seq=self.alternate_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME)
        self.assertEqual(b1, ann.break1)
        self.assertEqual(b2, ann.break2)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.s) + 'T' * (1699 - 1600 + 1) + 'G' * len(self.w) + 'T' * (1499 - 1300 + 1)
        expt += 'A' * len(self.z) + 'ATATGTAGA' + 'A' * len(self.b3) + 'G' * (1099 - 900 + 1)
        expt += 'C' * len(self.b2) + 'G' * (799 - 700 + 1) + 'T' * len(self.b1)

        self.assertEqual(expt, ft.seq)
        self.assertEqual(6, len(ft.exons))
        self.assertTrue(3, ft.exon_number(ft.exons[2]))
        self.assertTrue(3, ft.exon_number(ft.exons[3]))

    def test_build_two_transcript_inverted_translocation(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        t2 = PreTranscript(exons=[self.b1, self.b2, self.b3, self.b4, self.b5, self.b6], strand=STRAND.POS)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(ALT_REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='GATACATAT')
        ref = {REF_CHR: MockObject(seq=self.reference_sequence), ALT_REF_CHR: MockObject(seq=self.alternate_sequence)}
        ann = Annotation(bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME)
        self.assertEqual(b1, ann.break1)
        self.assertEqual(b2, ann.break2)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.s) + 'T' * (1699 - 1600 + 1) + 'G' * len(self.w) + 'T' * (1499 - 1300 + 1)
        expt += 'A' * len(self.z) + 'ATATGTATC' + 'C' * (1399 - 1200 + 1) + 'A' * len(self.b4)
        expt += 'C' * (1699 - 1500 + 1) + 'G' * len(self.b5) + 'C' * (2099 - 1800 + 1) + 'A' * len(self.b6)

        self.assertEqual(expt, ft.seq)
        self.assertEqual(6, len(ft.exons))
        self.assertTrue(3, ft.exon_number(ft.exons[2]))
        self.assertTrue(4, ft.exon_number(ft.exons[3]))


class TestSequenceFetching(unittest.TestCase):

    def setUp(self):
        self.gene = Gene(REF_CHR, 1, 900, strand=STRAND.POS)

        self.pre_transcript = PreTranscript(exons=[(101, 200), (301, 400), (501, 600), (701, 800)], gene=self.gene)
        self.gene.transcripts.append(self.pre_transcript)

        self.transcript = Transcript(self.pre_transcript, self.pre_transcript.generate_splicing_patterns()[0])
        self.pre_transcript.transcripts.append(self.transcript)

        self.translation = Translation(51, 350, self.transcript)
        self.transcript.translations.append(self.translation)

        self.spliced_seq = 'GGTGAATTTCTAGTTTGCCTTTTCAGCTAGGGATTAGCTTTTTAGGGGTCCCAATG' \
            'CCTAGGGAGATTTCTAGGTCCTCTGTTCCTTGCTGACCTCCAATAATCAGAAAATGCTGTGAAGGAAAAAC' \
            'AAAATGAAATTGCATTGTTTCTACCGGCCCTTTATCAAGCCCTGGCCACCATGATAGTCATGAATTCCAAT' \
            'TGTGTTGAAATCACTTCAATGTGTTTCTCTTCTTTCTGGGAGCTTACACACTCAAGTTCTGGATGCTTTGA' \
            'TTGCTATCAGAAGCCGTTAAATAGCTACTTATAAATAGCATTGAGTTATCAGTACTTTCATGTCTTGATAC' \
            'ATTTCTTCTTGAAAATGTTCATGCTTGCTGATTTGTCTGTTTGTTGAGAGGAGAATGTTC'

        self.domain = Domain(name=REF_CHR, regions=[(11, 20), (51, 60)], translation=self.translation)
        self.translation.domains.append(self.domain)

    def test_fetch_gene_seq_from_ref(self):
        expt = str(REFERENCE_GENOME[REF_CHR][0:900].seq).upper()
        self.assertEqual(expt, self.gene.get_seq(REFERENCE_GENOME))
        # gene seq should be the same if gene in on reverse strand b/c gene seq always given on pos
        self.gene.strand = STRAND.NEG
        self.assertEqual(expt, self.gene.get_seq(REFERENCE_GENOME))

    def test_fetch_gene_seq_from_stored(self):
        expt = 'AAA'
        self.gene.seq = expt
        self.assertEqual(expt, self.gene.get_seq(REFERENCE_GENOME))

    def test_fetch_gene_seq_force_uncached(self):
        expt = str(REFERENCE_GENOME[REF_CHR][0:900].seq).upper()
        self.gene.seq = 'AAA'
        self.assertEqual(expt, self.gene.get_seq(REFERENCE_GENOME, ignore_cache=True))

    def test_fetch_us_transcript_seq_from_ref(self):
        expt = str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper()
        self.assertEqual(expt, self.pre_transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_us_transcript_seq_from_ref_revcomp(self):
        self.gene.strand = STRAND.NEG
        expt = reverse_complement(str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper())
        self.assertEqual(expt, self.pre_transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_us_transcript_seq_from_stored(self):
        expt = 'AAA'
        self.pre_transcript.seq = expt
        self.assertEqual(expt, self.pre_transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_us_transcript_seq_from_parent_gene(self):
        self.gene.seq = 'A' * len(self.gene)
        self.assertEqual('A' * len(self.pre_transcript), self.pre_transcript.get_seq())

    def test_fetch_us_transcript_seq_from_parent_gene_revcomp(self):
        self.gene.seq = 'A' * len(self.gene)
        self.gene.strand = STRAND.NEG
        self.assertEqual('T' * len(self.pre_transcript), self.pre_transcript.get_seq())

    def test_fetch_us_transcript_seq_force_uncached(self):
        expt = str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper()
        self.pre_transcript.seq = 'AAA'
        self.assertEqual(expt, self.pre_transcript.get_seq(REFERENCE_GENOME, ignore_cache=True))

    def test_fetch_transcript_seq_from_ref(self):
        self.assertEqual(self.spliced_seq, self.transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_transcript_seq_from_ref_revcomp(self):
        self.gene.strand = STRAND.NEG
        self.assertEqual(reverse_complement(self.spliced_seq), self.transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_transcript_seq_from_stored(self):
        expt = 'AAA'
        self.transcript.seq = expt
        self.assertEqual(expt, self.transcript.get_seq(REFERENCE_GENOME))

    def test_fetch_transcript_seq_from_parent_ust(self):
        self.pre_transcript.seq = 'A' * len(self.pre_transcript)
        self.assertEqual('A' * len(self.transcript), self.transcript.get_seq())

    def test_fetch_transcript_seq_from_parent_gene(self):
        self.gene.seq = 'A' * len(self.gene)
        self.assertEqual('A' * len(self.transcript), self.transcript.get_seq())

    def test_fetch_transcript_seq_force_uncached(self):
        self.transcript.seq = 'AAA'
        self.assertEqual(self.spliced_seq, self.transcript.get_seq(REFERENCE_GENOME, ignore_cache=True))

    def test_fetch_translation_aa_seq_from_ref(self):
        cds = self.spliced_seq[self.translation.start - 1:self.translation.end]
        self.assertEqual(translate(cds), self.translation.get_aa_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_ref(self):
        cds = self.spliced_seq[self.translation.start - 1:self.translation.end]
        self.assertEqual(cds, self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_ref_revcomp(self):
        self.gene.strand = STRAND.NEG
        cdna = reverse_complement(self.spliced_seq)
        cds = cdna[self.translation.start - 1:self.translation.end]
        self.assertEqual(cds, self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_stored(self):
        expt = 'AAA'
        self.translation.seq = expt
        self.assertEqual(expt, self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_parent_transcript(self):
        self.transcript.seq = 'A' * len(self.transcript)
        self.assertEqual('A' * len(self.translation), self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_parent_ust(self):
        self.pre_transcript.seq = 'A' * len(self.pre_transcript)
        self.assertEqual('A' * len(self.translation), self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_from_parent_gene(self):
        self.gene.seq = 'A' * len(self.gene)
        self.assertEqual('A' * len(self.translation), self.translation.get_seq(REFERENCE_GENOME))

    def test_fetch_translation_cds_seq_force_uncached(self):
        self.translation.seq = 'AAA'
        cds = self.spliced_seq[self.translation.start - 1:self.translation.end]
        self.assertEqual(cds, self.translation.get_seq(REFERENCE_GENOME, ignore_cache=True))

    def test_fetch_domain_seq_from_ref(self):
        seqs = ['VPC*PPIIRK', 'C*NHFNVFLF']
        self.assertEqual(seqs, self.domain.get_seqs(REFERENCE_GENOME))


class TestStrandInheritance(unittest.TestCase):

    def setUp(self):
        self.gene = Gene('1', 1, 500, strand=STRAND.POS)
        pre_transcript = PreTranscript(gene=self.gene, exons=[(1, 100), (200, 300), (400, 500)])
        self.gene.unspliced_transcripts.append(pre_transcript)
        for spl in pre_transcript.generate_splicing_patterns():
            t = Transcript(pre_transcript, spl)
            pre_transcript.spliced_transcripts.append(t)
            tl = Translation(51, 250, t)
            t.translations.append(tl)

    def test_strand_gene(self):
        self.assertEqual(STRAND.POS, self.gene.get_strand())

    def test_strand_us_transcript(self):
        self.assertEqual(STRAND.POS, self.gene.unspliced_transcripts[0].get_strand())

    def test_strand_spl_transcript(self):
        self.assertEqual(STRAND.POS, self.gene.spliced_transcripts[0].get_strand())

    def test_strand_translation(self):
        self.assertEqual(STRAND.POS, self.gene.spliced_transcripts[0].translations[0].get_strand())


class TestCoordinateCoversion(unittest.TestCase):

    def setUp(self):
        self.gene = Gene('1', 15, 700, strand=STRAND.POS)

        self.pre_transcript = PreTranscript(gene=self.gene, exons=[(101, 200), (301, 400), (501, 600)])
        self.gene.unspliced_transcripts.append(self.pre_transcript)
        assert(1 == len(self.pre_transcript.generate_splicing_patterns()))

        spl = self.pre_transcript.generate_splicing_patterns()[0]
        self.transcript = Transcript(self.pre_transcript, spl)
        self.pre_transcript.spliced_transcripts.append(self.transcript)

        self.translation = Translation(51, 251, self.transcript)
        self.transcript.translations.append(self.translation)

        self.rev_gene = Gene('1', 15, 700, strand=STRAND.NEG)
        self.rev_ust = PreTranscript(gene=self.rev_gene, exons=[(101, 200), (301, 400), (501, 600)])
        self.gene.unspliced_transcripts.append(self.rev_ust)
        assert(1 == len(self.rev_ust.generate_splicing_patterns()))

        spl = self.rev_ust.generate_splicing_patterns()[0]
        self.rev_transcript = Transcript(self.rev_ust, spl)
        self.rev_ust.spliced_transcripts.append(self.rev_transcript)

        self.rev_translation = Translation(51, 251, self.rev_transcript)
        self.rev_transcript.translations.append(self.rev_translation)

    def test_cdna_to_genomic(self):
        self.assertEqual(150, self.transcript.convert_cdna_to_genomic(50))
        self.assertEqual(550, self.transcript.convert_cdna_to_genomic(250))

    def test_cdna_to_genomic_before(self):
        self.assertEqual(100, self.transcript.convert_cdna_to_genomic(-1))
        self.assertEqual(51, self.transcript.convert_cdna_to_genomic(-50))

    def test_cdna_to_genomic_after(self):
        self.assertEqual(650, self.transcript.convert_cdna_to_genomic(350))

    def test_cdna_to_genomic_revcomp(self):
        self.assertEqual(551, self.rev_transcript.convert_cdna_to_genomic(50))
        self.assertEqual(151, self.rev_transcript.convert_cdna_to_genomic(250))

    def test_genomic_to_cdna(self):
        self.assertEqual(50, self.transcript.convert_genomic_to_cdna(150))
        self.assertEqual(249, self.transcript.convert_genomic_to_cdna(549))

    def test_genomic_to_cdna_before(self):
        self.assertEqual((1, -1), self.transcript.convert_genomic_to_nearest_cdna(100))

    def test_genomic_to_cdna_after(self):
        self.assertEqual((300, 1), self.transcript.convert_genomic_to_nearest_cdna(601))

    def test_genomic_to_cdna_revcomp(self):
        self.assertEqual(50, self.rev_transcript.convert_genomic_to_cdna(551))
        self.assertEqual(250, self.rev_transcript.convert_genomic_to_cdna(151))

    def test_aa_to_cdna(self):
        self.assertEqual(Interval(51, 53), self.translation.convert_aa_to_cdna(1))
        self.assertEqual(Interval(249, 251), self.translation.convert_aa_to_cdna(67))

    def test_cdna_to_aa(self):
        self.assertEqual(1, self.translation.convert_cdna_to_aa(51))
        self.assertEqual(67, self.translation.convert_cdna_to_aa(251))
        with self.assertRaises(IndexError):
            self.translation.convert_cdna_to_aa(50)
        with self.assertRaises(IndexError):
            self.translation.convert_cdna_to_aa(252)

    def test_genomic_to_cds(self):
        self.assertEqual(1, self.translation.convert_genomic_to_cds(151))
        self.assertEqual(201, self.translation.convert_genomic_to_cds(551))

    def test_genomic_to_cds_3prime_utr(self):
        self.assertEqual(-1, self.translation.convert_genomic_to_cds(150))

    def test_genomic_to_cds_5prime_utr(self):
        self.assertEqual(202, self.translation.convert_genomic_to_cds(552))

    def test_genomic_to_cds_notation(self):
        self.assertEqual('1', self.translation.convert_genomic_to_cds_notation(151))
        self.assertEqual('201', self.translation.convert_genomic_to_cds_notation(551))

    def test_genomic_to_cds_notation_3prime_utr(self):
        self.assertEqual('-1', self.translation.convert_genomic_to_cds_notation(150))

    def test_genomic_to_cds_notation_5prime_utr(self):
        self.assertEqual('*1', self.translation.convert_genomic_to_cds_notation(552))

    def test_genomic_to_cds_notation_intronic_pos(self):
        self.assertEqual('50+2', self.translation.convert_genomic_to_cds_notation(202))

    def test_genomic_to_cds_notation_intronic_neg(self):
        self.assertEqual('51-2', self.translation.convert_genomic_to_cds_notation(299))

    def test_genomic_to_nearest_cdna_exonic(self):
        self.assertEqual((1, 0), self.transcript.convert_genomic_to_nearest_cdna(101))
        self.assertEqual((300, 0), self.transcript.convert_genomic_to_nearest_cdna(600))
        self.assertEqual((101, 0), self.transcript.convert_genomic_to_nearest_cdna(301))

    def test_genomic_to_nearest_cdna_intronic_pos(self):
        self.assertEqual((100, 10), self.transcript.convert_genomic_to_nearest_cdna(210))

    def test_genomic_to_nearest_cdna_intronic_neg(self):
        self.assertEqual((101, -2), self.transcript.convert_genomic_to_nearest_cdna(299))

    def test_genomic_to_nearest_cdna_rev_exonic(self):
        self.assertEqual((300, 0), self.rev_transcript.convert_genomic_to_nearest_cdna(101))
        self.assertEqual((1, 0), self.rev_transcript.convert_genomic_to_nearest_cdna(600))
        self.assertEqual((101, 0), self.rev_transcript.convert_genomic_to_nearest_cdna(400))

    def test_genomic_to_nearest_cdna_rev_intronic_pos(self):
        self.assertEqual((201, -10), self.rev_transcript.convert_genomic_to_nearest_cdna(210))

    def test_genomic_to_nearest_cdna_rev_intronic_neg(self):
        self.assertEqual((200, 2), self.rev_transcript.convert_genomic_to_nearest_cdna(299))


class TestUSTranscript(unittest.TestCase):

    def test___init__implicit_start(self):
        t = PreTranscript(gene=None, exons=[(1, 100), (200, 300), (400, 500)], strand=STRAND.POS)
        self.assertEqual(1, t.start)
        self.assertEqual(t.start, t.start)
        self.assertEqual(500, t.end)
        self.assertEqual(t.end, t.end)
        self.assertEqual(1, t[0])
        self.assertEqual(500, t[1])
        self.assertFalse(Interval.overlaps((0, 0), t))
        self.assertTrue(Interval.overlaps((1, 1), t))
        self.assertTrue(Interval.overlaps((1, 50), t))

    def test___init__strand_mismatch(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)

        with self.assertRaises(AssertionError):
            PreTranscript([(1, 100)], gene=g, strand=STRAND.NEG)

    def test___init__overlapping_exon_error(self):
        with self.assertRaises(AttributeError):
            PreTranscript(exons=[Exon(1, 15), Exon(10, 20)])

    def test_exon_number(self):
        t = PreTranscript(gene=None, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        for i, e in enumerate(t.exons):
            self.assertEqual(i + 1, t.exon_number(e))

        t = PreTranscript(gene=None, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        for i, e in enumerate(sorted(t.exons, key=lambda x: x.start, reverse=True)):
            self.assertEqual(i + 1, t.exon_number(e))


class TestDomain(unittest.TestCase):

    def test___init__region_error(self):
        with self.assertRaises(AttributeError):
            Domain('name', [(1, 3), (4, 3)])

    def test_get_seq_from_ref(self):
        ref = {'1': MockObject(seq='CCCTAATCCCCTTT')}
        g = Gene('1', 1, 16, strand=STRAND.NEG)
        t = PreTranscript(exons=[(2, 5), (7, 15)], gene=g)
        tl = Translation(4, 11, t, [])
        d = Domain('name', [(1, 2)], translation=tl)
        self.assertEqual([translate('GGGGAT')], d.get_seqs(ref))

    def test_get_seq_from_translation_seq(self):
        t = PreTranscript(exons=[(2, 5), (7, 15)], seq='CCCTAATCCCCTTT', strand=STRAND.NEG)
        tl = Translation(4, 11, t, [])
        d = Domain('name', [(1, 2)], translation=tl)
        self.assertEqual([translate('TAATCC')], d.get_seqs())

    def test_align_seq(self):
        regions = [
            DomainRegion(216, 261, 'DVNECITGSHSCRLGESCINTVGSFRCQRDSSCGTGYELTEDNSCK'),
            DomainRegion(262, 307, 'DIDECESGIHNCLPDFICQNTLGSFRCRPKLQCKSGFIQDALGNCI'),
            DomainRegion(308, 355, 'DINECLSISAPCPIGHTCINTEGSYTCQKNVPNCGRGYHLNEEGTRCV'),
            DomainRegion(356, 398, 'DVDECAPPAEPCGKGHRCVNSPGSFRCECKTGYYFDGISRMCV'),
            DomainRegion(399, 440, 'DVNECQRYPGRLCGHKCENTLGSYLCSCSVGFRLSVDGRSCE'),
            DomainRegion(441, 480, 'DINECSSSPCSQECANVYGSYQCYCRRGYQLSDVDGVTCE'),
            DomainRegion(481, 524, 'DIDECALPTGGHICSYRCINIPGSFQCSCPSSGYRLAPNGRNCQ'),
            DomainRegion(525, 578, 'DIDECVTGIHNCSINETCFNIQGGFRCLAFECPENYRRSAATLQQEKTDTVRCI')
        ]
        print(regions)
        refseq = 'MERAAPSRRVPLPLLLLGGLALLAAGVDADVLLEACCADGHRMATHQKDCSLPYATESKE' \
            'CRMVQEQCCHSQLEELHCATGISLANEQDRCATPHGDNASLEATFVKRCCHCCLLGRAAQ' \
            'AQGQSCEYSLMVGYQCGQVFQACCVKSQETGDLDVGGLQETDKIIEVEEEQEDPYLNDRC' \
            'RGGGPCKQQCRDTGDEVVCSCFVGYQLLSDGVSCEDVNECITGSHSCRLGESCINTVGSF' \
            'RCQRDSSCGTGYELTEDNSCKDIDECESGIHNCLPDFICQNTLGSFRCRPKLQCKSGFIQ' \
            'DALGNCIDINECLSISAPCPIGHTCINTEGSYTCQKNVPNCGRGYHLNEEGTRCVDVDEC' \
            'APPAEPCGKGHRCVNSPGSFRCECKTGYYFDGISRMCVDVNECQRYPGRLCGHKCENTLG' \
            'SYLCSCSVGFRLSVDGRSCEDINECSSSPCSQECANVYGSYQCYCRRGYQLSDVDGVTCE' \
            'DIDECALPTGGHICSYRCINIPGSFQCSCPSSGYRLAPNGRNCQDIDECVTGIHNCSINE' \
            'TCFNIQGGFRCLAFECPENYRRSAATLQQEKTDTVRCIKSCRPNDVTCVFDPVHTISHTV' \
            'ISLPTFREFTRPEEIIFLRAITPPHPASQANIIFDITEGNLRDSFDIIKRYMDGMTVGVV' \
            'RQVRPIVGPFHAVLKLEMNYVVGGVVSHRNVVNVHIFVSEYWF'

        d = Domain('name', regions)
        self.assertTrue(len(refseq) >= 578)
        match, total, temp = d.align_seq(refseq)
        self.assertEqual(sum([len(d.seq) for d in regions]), total)
        self.assertEqual(total, match)
        self.assertEqual(len(regions), len(temp))
        for dr1, dr2 in zip(temp, regions):
            self.assertEqual(dr1.start, dr2.start)
            self.assertEqual(dr1.end, dr2.end)
            self.assertEqual(dr1.seq, dr2.seq)

        refseq = 'MHRPPRHMGNKAMEPMDSPLMSAIPRLRPLQPMGRPPMQLLMDSLPLVILLQLPPRHTASLSRGMALVLMIPPL' \
            'LQSPPPRPPMQLSLHMALSLLIQPMGSSQQPLHLQAIPLHSRLVMIRAVTLSRTPMGNRAAMDSRVAMVNKAAMGSSLP' \
            'LVTHPKLDPTAKLQVNIANRAAATGSRTLLMTQSEEELGAIT*'

        d = 'QAAAQQGYSAYTAQPTQGYAQTTQAYGQQSYGTYGQPTDVSYTQAQTTATYGQTAYATSYGQPPTGYTTPTAPQAYSQP' \
            'VQGYGTGAYDTTTATVTTTQASYAAQSAYGTQPAYPAYGQQPAATAPTSYSSTQPTSYDQSSYSQQNTYGQPSSYGQQS' \
            'SYGQQSSYGQQPPTSYPPQTGSYSQAPSQYSQQSSSYGQQSSFRQ'

        dom = Domain('name', [DomainRegion(1, len(d), d)])
        with self.assertRaises(UserWarning):
            dom.align_seq(refseq)

        seq = 'MASTDYSTYSQAAAQQGYSAYTAQPTQGYAQTTQAYGQQSYGTYGQPTDVSYTQAQTTATYGQTAYATSYGQPPTGY' \
            'TTPTAPQAYSQPVQGYGTGAYDTTTATVTTTQASYAAQSAYGTQPAYPAYGQQPAATAPTRPQDGNKPTETSQPQSSTG' \
            'GYNQPSLGYGQSNYSYPQVPGSYPMQPVTAPPSYPPTSYSSTQPTSYDQSSYSQQNTYGQPSSYGQQSSYGQQSSYGQQ' \
            'PPTSYPPQTGSYSQAPSQYSQQSSSYGQQNPSYDSVRRGAWGNNMNSGLNKSPPLGGAQTISKNTEQRPQPDPYQILGP' \
            'TSSRLANPGSGQIQLWQFLLELLSDSANASCITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYYDKNI' \
            'MTKVHGKRYAYKFDFHGIAQALQPHPTESSMYKYPSDISYMPSYHAHQQKVNFVPPHPSSMPVTSSSFFGAASQYWTSP' \
            'TGGIYPNPNVPRHPNTHVPSHLGSYY'
        d = 'IYVQGLNDSVTLDDLADFFKQCGVVKMNKRTGQPMIHIYLDKETGKPKGDATVSYEDPPTAKAAVEWFDGKDFQGSKLK'
        dom = Domain('name', [DomainRegion(1, len(d), d)])

        with self.assertRaises(UserWarning):
            m, t, regions = dom.align_seq(seq)


class TestBioInterval(unittest.TestCase):

    def test___eq__(self):
        a = BioInterval(REF_CHR, 1, 2)
        b = BioInterval(REF_CHR, 1, 2)
        c = BioInterval('test2', 1, 2)
        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, None)
        self.assertNotEqual(a, c)


class TestGene(unittest.TestCase):

    def test___hash__(self):
        g1 = Gene(REF_CHR, 1, 2, 'name1', STRAND.POS)
        g2 = Gene(REF_CHR, 1, 2, 'name2', STRAND.POS)
        h = set([g1, g2])
        self.assertEqual(2, len(h))

    def test___eq__(self):
        g1 = Gene(REF_CHR, 1, 2, 'name1', STRAND.POS)
        g2 = Gene(REF_CHR, 1, 2, 'name2', STRAND.POS)
        self.assertNotEqual(g1, g2)
        g3 = Gene('test2', 1, 2, 'name1', STRAND.POS)
        self.assertNotEqual(g1, g3)
        self.assertNotEqual(g3, g1)
        self.assertNotEqual(g1, None)
        self.assertNotEqual(None, g1)

    def test_get_seq(self):
        ref = {'1': MockObject(seq='AACCCTTTGGG')}
        g = Gene('1', 3, 8, strand=STRAND.POS)
        self.assertEqual('CCCTTT', g.get_seq(ref))
        g = Gene(REF_CHR, 2836, 4144, strand=STRAND.POS)
        seq = (
            'GCAACTATATAATCTGTGGGAATATCTCCTTTTACACCTAGCCCTACTTCTGTCTGGCTACAGTCATTTATCTGGCTTTGGGAAATGTGACCACAGAATCAGATAT'
            'ATACATGAGATTAAATAATACATGTGTATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAATCATACA'
            'GTAAATAGATTTTTGAATTTAATCTAGTACCTAAATAATCAGAGTAGGGAGGTTAGATATTAAAATCAGGCTAAAGATATAGGCAACATGGATCTAGAAAACATGG'
            'ATTGCATGGCCATTTCACTTAGAGTTCATGGGCTTGGAATCTCTATTAACATAACTTTTACAATGTTAGAATTTGTTCCCATATTAATGAGGGAAAAACAAACAAT'
            'TACCCTGAGTATCTGAAGCTCCAGATCTCATTTTCCAGTCAAAATCTCTGATAGGTAAACAACCTGAAAAAGTAGCCACAACTCACTGAGGTGATAACCTCATTTG'
            'CTTAAGAGAATGTAATTGTTTTTATGATTTTTTTTATCCCAGGAAAACATTGAAAAAAAGTTTAGAGATGATGAAGTATATGAAAACTATAATATTTATACTTTAG'
            'AGATGTGatatttatttataattgtattagtatttaaatataGATTAGCATTTTACATTCCAATTTTCAATGTGTAACAGAATATTTTAGATATTGGGGTTGTTTT'
            'TTAGTTGAAATAATAAGCGGTTTTACCGAGTTGCCAGTAGTGGTTTAACATTGAAGATAATTTAACATTCATGATTTTGTGAGTTTAATTTATTAGCTCTATAAGG'
            'GTTGTTTAAGTACTCTGAAGGCTTTATTTGTTAGTCCGATAATTAAAATGTTCATAAAGATAATTCAACATATTAAATTTGTAAATGTAGTTTAAAATCTTTAAGG'
            'GAGTTTAATTAACTAAGTTGTAAATGGACAAAACATTAATCAAAGTCCCCCTTAAAAATAATTTTTAATGTACTAGATTTATAAATAGAACAACAAGATTTCTAAT'
            'TTAAACTCAAAAATTTTTTAAATTGGTTAACAATTTAACATAATATGCTGCACATTAATTCAGAATATGAAATCTTATATGTAGTCCTTTTTACATTCAAGAATCA'
            'CATCGATAAACATCACAAAATGACTACTGGTAACCACTATGAAACTCTTTAAGCGGTAGGTCCTGTATGAATTTTACTCCTCATGATTTGAAGATTATGCATAAAT'
            'TCCTTCTTCCTGTTATTTTGTTTCCAATTTAGTCTTT').upper()
        self.assertEqual(seq, g.get_seq(REFERENCE_GENOME))


class TestAnnotationGathering(unittest.TestCase):

    def test_overlapping_transcripts(self):
        b = Breakpoint('C', 1000, strand=STRAND.POS)
        g = Gene('C', 1, 9999, 'gene1', STRAND.POS)
        t = PreTranscript(exons=[(100, 199), (500, 699), (1200, 1300)], gene=g)
        g.transcripts.append(t)
        self.assertTrue(Interval.overlaps(b, t))
        t = PreTranscript(exons=[(100, 199), (500, 699), (800, 900)], gene=g)
        g.transcripts.append(t)
        h = Gene('C', 1, 9999, 'gene1', STRAND.NEG)
        t = PreTranscript(exons=[(100, 199), (500, 699), (1200, 1300)], gene=h, strand=STRAND.NEG)
        h.transcripts.append(t)
        d = {'C': [g, h]}
        tlist = overlapping_transcripts(d, b)
        self.assertEqual(1, len(tlist))

    def test_breakpoint_within_gene(self):
        b = Breakpoint(REF_CHR, 150, 150)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(1, len(pos))
        self.assertEqual(1, len(neg))
        self.assertEqual(STRAND.POS, pos[0].get_strand())
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)
        self.assertEqual(STRAND.NEG, neg[0].get_strand())

    def test_breakpoint_overlapping_gene(self):
        b = Breakpoint(REF_CHR, 150, 230)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

        b = Breakpoint(REF_CHR, 150, 225, strand=STRAND.POS)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(100, pos[0].start)
        self.assertEqual(200, pos[0].end)
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

        b = Breakpoint(REF_CHR, 375, 425, strand=STRAND.POS)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(300, pos[0].start)
        self.assertEqual(400, pos[0].end)
        self.assertEqual(401, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

    def test_breakpoint_overlapping_mutliple_genes_and_intergenic(self):
        b = Breakpoint(REF_CHR, 150, 275)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(2, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(249, neg[0].end)

    def test_breakpoint_overlapping_mutliple_pos_genes(self):
        b = Breakpoint(REF_CHR, 575, 625)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

    def test_breakpoint_overlapping_mutliple_genes(self):
        b = Breakpoint(REF_CHR, 300, 350)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(1, len(pos))
        self.assertEqual(1, len(neg))

    def test_intrachromosomal(self):
        b1 = Breakpoint(REF_CHR, 150, 225, strand=STRAND.POS)
        b2 = Breakpoint(REF_CHR, 375, 425, strand=STRAND.POS)
        bpp = BreakpointPair(b1, b2, protocol=PROTOCOL.GENOME, event_type=SVTYPE.DEL)
        ann_list = sorted(_gather_annotations(REFERENCE_ANNOTATIONS, bpp),
                          key=lambda x: (x.break1, x.break2))
        self.assertEqual(5, len(ann_list))
        first = ann_list[0]
        self.assertEqual(1, len(first.encompassed_genes))
        self.assertEqual(0, len(first.genes_proximal_to_break1))
        self.assertEqual(1, len(first.genes_proximal_to_break2))
        self.assertEqual(0, len(first.genes_overlapping_break1))
        self.assertEqual(0, len(first.genes_overlapping_break2))
        near, dist = list(first.genes_proximal_to_break2)[0]
        self.assertEqual(50, dist)
        self.assertEqual(2, len(ann_list[1].encompassed_genes))

    def test_interchromosomal(self):
        raise unittest.SkipTest('TODO')

    def test_intrachromosomal_within_gene_inversion(self):
        raise unittest.SkipTest('TODO')
        g = Gene(REF_CHR, 1000, 3000, strand=STRAND.POS)
        t = PreTranscript(gene=g, exons=[(1001, 1100), (1501, 1600), (2001, 2100), (2501, 2600)])
        g.transcripts.append(t)
        ref = {
            REF_CHR: [g]
        }
        b1 = Breakpoint(REF_CHR, 1250, strand=STRAND.POS)
        b2 = Breakpoint(REF_CHR, 2250, strand=STRAND.NEG)
        bpp = BreakpointPair(b1, b2)
        ann_list = sorted(_gather_annotations(ref, bpp),
                          key=lambda x: (x.break1, x.break2))
        self.assertEqual(1, len(ann_list))
        self.assertEqual(ann_list[0].transcript1, ann_list[0].transcript2)


class TestAnnotate(unittest.TestCase):

    def test_reference_name_eq(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        self.assertEqual(first, second)

    def test_reference_name_set(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        d = {first, second}
        self.assertEqual(1, len(d))

    def test_reference_name_dict(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        d = {first: 1}
        d[second] = 2
        print(d)
        self.assertEqual(1, len(d))
        d = {first: 1, second: 2}
        self.assertEqual(1, len(d))

    def test_loading_json_annotations(self):
        annotations = load_reference_genes(get_data('mock_reference_annotations.json'))
        self.assertEqual(1, len(annotations.keys()))
        self.assertEqual(1, len(list(annotations.values())[0]))

    def test_loading_annotations_not_found(self):
        with self.assertRaises(FileNotFoundError):
            load_reference_genes('file.other')

    def test_determine_prime(self):
        tneg = PreTranscript(exons=[(3, 4)], strand=STRAND.NEG)
        tpos = PreTranscript(exons=[(3, 4)], strand=STRAND.POS)
        bleft = Breakpoint(REF_CHR, 1, 2, orient=ORIENT.LEFT)
        bright = Breakpoint(REF_CHR, 1, 2, orient=ORIENT.RIGHT)
        # positive left should be five prime
        self.assertEqual(PRIME.FIVE, determine_prime(tpos, bleft))
        # positive right should be three prime
        self.assertEqual(PRIME.THREE, determine_prime(tpos, bright))
        # negative left should be three prime
        self.assertEqual(PRIME.THREE, determine_prime(tneg, bleft))
        # negative right should be five prime
        self.assertEqual(PRIME.FIVE, determine_prime(tneg, bright))

        with self.assertRaises(NotSpecifiedError):
            bleft.orient = ORIENT.NS
            determine_prime(tpos, bleft)

        with self.assertRaises(NotSpecifiedError):
            determine_prime(tneg, bleft)

        with self.assertRaises(NotSpecifiedError):
            tpos.strand = STRAND.NS
            determine_prime(tpos, bright)

    def test_calculate_orf_nested(self):
        seq = 'ATGAACACGCAGGAACACCCACGCTCCCGCTCCTCCTACT' \
            'CCAGGGGGAGACCGAACCCGAGAGCGACATCCGGAGCTGG' \
            'AAGCCGCTGCAACGGGGGCGCCGGCCTCCCTGCCCCGCAG' \
            'TCTCCCGCCCTTTCCCAGACACCGGAGGATCCCAGACCAG' \
            'CCCAGTATTCCAGACGGAGAGGGTGCCCGAGTCCTCCGCA' \
            'GGGCCCCGCCTCCCTAGGCGCTCGGCGCTGGAGACCTAGG' \
            'GGAGCGCGCGAGAGAAACGCTGAGCTCCACCCGGCTGGGG' \
            'CGGAGGCTGAGCACTATACCCGGGACTCGAACCCCGCGCC' \
            'CTCGTCCCCCAGCTTTAGGCATCCCACGGCGAGCCCTGGC' \
            'CCAACGCCCCCGCGCCGCGTCTTCCCCCGCACCCCCGCCC' \
            'GGAGCCGAGTCCCGGGTTGCCAGCGCGCCCATCCCGCCCC' \
            'CAACAGGTCTCGGGCAGCGGGGAGCGCGCTGCCTGTCCCG' \
            'CAGGCGCTGGGAGGCGAGCCGGAGCGAAGCCGAGGGCCGG' \
            'CCGCGCGCGCGGGAGGGGTCGCCGAGGGCTGCGAGCCCGG' \
            'CCACTTACGCGCGGGCTCTGGGGAACCCCATGGAGAGGAG' \
            'CACGTCCAGCGCCGATCCATGCTTGATGGTGCCGGGGCGC' \
            'TGTTGGCGGTTCCTCCGGGGGGTGACTTTGCTGTACAGCT' \
            'CCTCTCTCGCAGCCATGCCGAGCGGACTGGGGTGGCCGTA' \
            'CTGAGCCATGGCTCAGCGGCCAGCCATCGCCAGGGAAGGA' \
            'GCCCGAGCCGGATGGTTCGGGAAGGAGGAGGGAGTCGGGC' \
            'TGGGCCAGGGGAGGGGGCTCGGGGACCCAGAGCCAGGCAG' \
            'GAGGCGGCGGCAGCGGAGGCGGCGCTCCGAGCTTCCCCAG' \
            'TGTCGGGACTCTGA'
        orfs = calculate_orf(seq)
        for orf in orfs:
            self.assertEqual('ATG', seq[orf.start - 1:orf.start + 2])
        orfs = sorted(orfs)
        self.assertEqual(2, len(orfs))
        self.assertEqual(Interval(1, 894), orfs[0])
        self.assertEqual(Interval(590, 724), orfs[1])

        seq = 'AAGGAGAGAAAATGGCGTCCACGGATTACAGTACCTATAGCCAAGCTGCAGCGCAGCAGGGCTACAGTGCTTACACCGCCCAGCCCACTCAAGGATATGC' \
              'ACAGACCACCCAGGCATATGGGCAACAAAGCTATGGAACCTATGGACAGCCCACTGATGTCAGCTATACCCAGGCTCAGACCACTGCAACCTATGGGCAG' \
              'ACCGCCTATGCAACTTCTTATGGACAGCCTCCCACTGGTTATACTACTCCAACTGCCCCCCAGGCATACAGCCAGCCTGTCCAGGGGTATGGCACTGGTG' \
              'CTTATGATACCACCACTGCTACAGTCACCACCACCCAGGCCTCCTATGCAGCTCAGTCTGCATATGGCACTCAGCCTGCTTATCCAGCCTATGGGCAGCA' \
              'GCCAGCAGCCACTGCACCTACAAGCTATTCCTCTACACAGCCGACTAGTTATGATCAGAGCAGTTACTCTCAGCAGAACACCTATGGGCAACCGAGCAGC' \
              'TATGGACAGCAGAGTAGCTATGGTCAACAAAGCAGCTATGGGCAGCAGCCTCCCACTAGTTACCCACCCCAAACTGGATCCTACAGCCAAGCTCCAAGTC' \
              'AATATAGCCAACAGAGCAGCAGCTACGGGCAGCAGAGTTCATTCCGACAGGACCACCCCAGTAGCATGGGTGTTTATGGGCAGGAGTCTGGAGGATTTTC' \
              'CGGACCAGGAGAGAACCGGAGCATGAGTGGCCCTGATAACCGGGGCAGGGGAAGAGGGGGATTTGATCGTGGAGGCATGAGCAGAGGTGGGCGGGGAGGA' \
              'GGACGCGGTGGAATGGGCAGCGCTGGAGAGCGAGGTGGCTTCAATAAGCCTGGTGGACCCATGGATGAAGGACCAGATCTTGATCTAGGCCCACCTGTAG' \
              'ATCCAGATGAAGACTCTGACAACAGTGCAATTTATGTACAAGGATTAAATGACAGTGTGACTCTAGATGATCTGGCAGACTTCTTTAAGCAGTGTGGGGT' \
              'TGTTAAGATGAACAAGAGAACTGGGCAACCCATGATCCACATCTACCTGGACAAGGAAACAGGAAAGCCCAAAGGCGATGCCACAGTGTCCTATGAAGAC' \
              'CCACCCACTGCCAAGGCTGCCGTGGAATGGTTTGATGGGAAAGATTTTCAAGGGAGCAAACTTAAAGTCTCCCTTGCTCGGAAGAAGCCTCCAATGAACA' \
              'GTATGCGGGGTGGTCTGCCACCCCGTGAGGGCAGAGGCATGCCACCACCACTCCGTGGAGGTCCAGGAGGCCCAGGAGGTCCTGGGGGACCCATGGGTCG' \
              'CATGGGAGGCCGTGGAGGAGATAGAGGAGGCTTCCCTCCAAGAGGACCCCGGGGTTCCCGAGGGAACCCCTCTGGAGGAGGAAACGTCCAGCACCGAGCT' \
              'GGAGACTGGCAGTGTCCCAATCCGGGTTGTGGAAACCAGAACTTCGCCTGGAGAACAGAGTGCAACCAGTGTAAGGCCCCAAAGCCTGAAGGCTTCCTCC' \
              'CGCCACCCTTTCCGCCCCCGGGTGGTGATCGTGGCAGAGGTGGCCCTGGTGGCATGCGGGGAGGAAGAGGTGGCCTCATGGATCGTGGTGGTCCCGGTGG' \
              'AATGTTCAGAGGTGGCCGTGGTGGAGACAGAGGTGGCTTCCGTGGTGGCCGGGGCATGGACCGAGGTGGCTTTGGTGGAGGAAGACGAGGTGGCCCTGGG' \
              'GGGCCCCCTGGACCTTTGATGGAACAGATGGGAGGAAGAAGAGGAGGACGTGGAGGACCTGGAAAAATGGATAAAGGCGAGCACCGTCAGGAGCGCAGAG' \
              'ATCGGCCCTACTAGATGCAGAGACCCCGCAGAGCTGCATTGACTACCAGATTTATTTTTTAAACCAGAAAATGTTTTAAATTTATAATTCCATATTTATA' \
              'ATGTTGGCCACAACATTATGATTATTCCTTGTCTGTACTTTAGTATTTTTCACCATTTGTGAAGAAACATTAAAACAAGTTAAATGGTA'

        orfs = calculate_orf(seq)
        for orf in orfs:
            self.assertEqual('ATG', seq[orf.start - 1:orf.start + 2])


class TestAnnotateEvents(unittest.TestCase):

    def test_annotate_events(self):
        reference_annotations = load_reference_genes(get_data('mock_reference_annotations.full.tsv'))
        b1 = Breakpoint('fakereference9', 658, orient=ORIENT.RIGHT, strand=STRAND.POS)
        b2 = Breakpoint('fakereference9', 10237, orient=ORIENT.RIGHT, strand=STRAND.NEG)
        bpp = BreakpointPair(
            b1, b2, stranded=True, opposing_strands=True, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME)
        annotations = annotate_events(
            [bpp], reference_genome=REFERENCE_GENOME, annotations=reference_annotations, filters=[])
        self.assertEqual(4, len(annotations))
        self.assertEqual(STRAND.POS, annotations[0].transcript1.get_strand())
        self.assertEqual(STRAND.NEG, annotations[0].transcript2.get_strand())
        self.assertEqual('ENST00000375851', annotations[0].transcript1.name)
        self.assertEqual(None, annotations[0].transcript2.name)
        for ann in annotations:
            print(ann.transcript1, ann.transcript2)
        annotations = annotate_events(
            [bpp], reference_genome=REFERENCE_GENOME, annotations=reference_annotations)
        self.assertEqual(2, len(annotations))
