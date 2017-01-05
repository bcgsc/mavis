import unittest
from structural_variant.annotate import *
from structural_variant.constants import reverse_complement
from structural_variant.constants import STRAND
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from tests import REFERENCE_ANNOTATIONS_FILE, MockSeq


REFERENCE_ANNOTATIONS = None


def setUpModule():
    global REFERENCE_ANNOTATIONS
    REFERENCE_ANNOTATIONS = load_reference_genes(REFERENCE_ANNOTATIONS_FILE)
    print('loaded {} annotations', sum([len(l) for l in REFERENCE_ANNOTATIONS.values()]))


class TestFusionTranscript(unittest.TestCase):

    def setUp(self):
        self.x = Exon(100, 199)  # C
        self.y = Exon(500, 599)  # G
        self.z = Exon(1200, 1299)  # T
        self.w = Exon(1500, 1599)  # C
        self.s = Exon(1700, 1799)  # G
        # introns: 99, 300, 600, 200, 100, ...
        reference_sequence = 'A' * 99 + 'C' * 100 + 'A' * 300 + 'G' * 100
        reference_sequence += 'A' * 600 + 'T' * 100 + 'A' * 200 + 'C' * 100
        reference_sequence += 'A' * 100 + 'G' * 100 + 'A' * 200
        self.reference_sequence = reference_sequence

    def test_determine_prime(self):
        tneg = Transcript(1, 2, genomic_start=3, genomic_end=4, strand=STRAND.NEG)
        tpos = Transcript(1, 2, genomic_start=3, genomic_end=4, strand=STRAND.POS)
        bleft = Breakpoint('test', 1, 2, orient=ORIENT.LEFT)
        bright = Breakpoint('test', 1, 2, orient=ORIENT.RIGHT)
        # positive left should be five prime
        self.assertEqual(PRIME.FIVE, FusionTranscript.determine_prime(tpos, bleft))
        # positive right should be three prime
        self.assertEqual(PRIME.THREE, FusionTranscript.determine_prime(tpos, bright))
        # negative left should be three prime
        self.assertEqual(PRIME.THREE, FusionTranscript.determine_prime(tneg, bleft))
        # negative right should be five prime
        self.assertEqual(PRIME.FIVE, FusionTranscript.determine_prime(tneg, bright))

        with self.assertRaises(AttributeError):
            bleft.orient = ORIENT.NS
            FusionTranscript.determine_prime(tpos, bleft)

        with self.assertRaises(AttributeError):
            FusionTranscript.determine_prime(tneg, bleft)

        with self.assertRaises(AttributeError):
            tpos._strand = STRAND.NS
            FusionTranscript.determine_prime(tpos, bright)

    def test__pull_exons_left_pos_intronic(self):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 700, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (700 - 600 + 1)
        self.assertEqual(expt, seq)
        self.assertEqual(2, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)

    def test__pull_exons_left_pos_intronic_splice(self):
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 201, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 100 + 'A' * 2
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(False, e.intact_end_splice)

    def test__pull_exons_left_pos_exonic(self):
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 199, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(False, e.intact_end_splice)

    def test__pull_exons_left_pos_exonic_splice(self):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 101, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'C' * 2
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(1, e.start)
        self.assertEqual(2, e.end)
        self.assertEqual(False, e.intact_start_splice)
        self.assertEqual(False, e.intact_end_splice)

    def test__pull_exons_right_pos_intronic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 1600, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        self.assertEqual(expt, seq)
        self.assertEqual(1, len(new_exons))

        b = Breakpoint('test', 300, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (499 - 300 + 1) + 'G' * 100 + 'A' * (1199 - 600 + 1) + 'T' * 100
        expt += 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100

        self.assertEqual(expt, seq)
        self.assertEqual(4, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)

    def test__pull_exons_right_pos_intronic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'AA' + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)

    def test__pull_exons_right_pos_exonic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 1201, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'T' * 99 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)

    def test__pull_exons_right_pos_exonic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b = Breakpoint('test', 1298, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'TT' + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(False, e.intact_start_splice)
        self.assertEqual(False, e.intact_end_splice)

    def test__pull_exons_right_neg_intronic(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b = Breakpoint('test', 700, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'A' * (1199 - 700 + 1) + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100
        expt += 'A' * (1699 - 1600 + 1) + 'G' * 100
        expt = reverse_complement(expt)
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual('C' * 100, seq[e.start - 1:e.end])
        e = new_exons[1][0]
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual('G' * 100, seq[e.start - 1:e.end])

    def test__pull_exons_right_neg_intronic_splice(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b = Breakpoint('test', 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, self.reference_sequence)
        expt = 'AA' + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        expt = reverse_complement(expt)
        self.assertEqual(expt, seq)
        self.assertEqual(3, len(new_exons))
        e = new_exons[0][0]
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)
        self.assertEqual(1, e.start)
        self.assertEqual(100, e.end)
        self.assertEqual('C' * 100, seq[e.start - 1:e.end])
        e = new_exons[1][0]
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(True, e.intact_end_splice)
        self.assertEqual(201, e.start)
        self.assertEqual(300, e.end)
        self.assertEqual('G' * 100, seq[e.start - 1:e.end])
        e = new_exons[2][0]
        self.assertEqual(True, e.intact_start_splice)
        self.assertEqual(False, e.intact_end_splice)
        self.assertEqual(501, e.start)
        self.assertEqual(600, e.end)
        self.assertEqual('A' * 100, seq[e.start - 1:e.end])

    def test_build_single_transcript_indel(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint('test', 599, orient=ORIENT.LEFT)
        b2 = Breakpoint('test', 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_sequence='ATCGATCG')
        ref = {'test': MockSeq(self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DEL)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)

        self.assertEqual(expt, ft.sequence)
        self.assertEqual(5, len(ft.exons))

        for i, ex in enumerate(t.exons):
            n = ft.exons[i]
            print(i, ex, (ex.start, ex.end), n, (n.start, n.end), ft.exon_mapping[n], (ft.exon_mapping[n].start, ft.exon_mapping[n].end))
            self.assertEqual(ex, ft.exon_mapping[n])

        self.assertEqual(1, ft.exons[0].start)
        self.assertEqual(100, ft.exons[0].end)

        splice_pattern = [(True, True), (True, False), (False, True), (True, True), (True, True)]
        char_pattern = [x * 100 for x in ['C', 'G', 'T', 'C', 'G']]

        for i in range(0, len(splice_pattern)):
            s, t = splice_pattern[i]
            ex = ft.exons[i]
            self.assertEqual(s, ex.intact_start_splice)
            self.assertEqual(t, ex.intact_end_splice)
            temp = ft.sequence[ex.start - 1:ex.end]
            self.assertEqual(char_pattern[i], ft.sequence[ex.start - 1:ex.end])

    def test_build_single_transcript_inversion(self):
        pass  # TODO

    def test_build_single_transcript_duplication_pos(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.POS)
        b1 = Breakpoint('test', 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint('test', 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_sequence='ATCGATCG')
        ref = {'test': MockSeq(self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'T' * len(self.z) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        print(expt)
        print(ft.sequence)
        self.assertEqual(expt, ft.sequence)
        self.assertEqual(6, len(ft.exons))
        self.assertTrue(ft.exons[2].intact_start_splice)
        self.assertTrue(ft.exons[3].intact_end_splice)
        self.assertFalse(ft.exons[2].intact_end_splice)
        self.assertFalse(ft.exons[3].intact_start_splice)

    def test_build_single_transcript_duplication_neg(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = Transcript(50, 249, exons=[self.x, self.y, self.z, self.w, self.s], strand=STRAND.NEG)
        b1 = Breakpoint('test', 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint('test', 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_sequence='ATCGATCG')
        ref = {'test': MockSeq(self.reference_sequence)}
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP)
        ft = FusionTranscript.build(ann, ref)

        expt = 'C' * len(self.x) + 'A' * (499 - 200 + 1) + 'G' * len(self.y) + 'A' * (1199 - 600 + 1)
        expt += 'T' * len(self.z) + 'ATCGATCG' + 'T' * len(self.z)
        expt += 'A' * (1499 - 1300 + 1) + 'C' * len(self.w) + 'A' * (1699 - 1600 + 1) + 'G' * len(self.s)
        expt = reverse_complement(expt)
        print(expt)
        print(ft.sequence)
        self.assertEqual(expt, ft.sequence)
        self.assertEqual(6, len(ft.exons))
        self.assertTrue(ft.exons[2].intact_start_splice)
        self.assertTrue(ft.exons[3].intact_end_splice)
        self.assertFalse(ft.exons[2].intact_end_splice)
        self.assertFalse(ft.exons[3].intact_start_splice)
        self.assertEqual(3, ft.exon_number(ft.exons[2]))
        self.assertEqual(3, ft.exon_number(ft.exons[3]))
    
    def test_build_two_transcript_inversion_5prime_pos(self):
        pass  # TODO

    def test_build_two_transcript_inversion_5prime_neg(self):
        pass  # TODO

    def test_build_two_transcript_translocation(self):
        pass  # TODO

    def test_build_two_transcript_inverted_translocation(self):
        pass  # TODO

    def test_build_two_transcript_duplication_pos(self):
        pass  # TODO

    def test_build_two_transcript_duplication_neg(self):
        pass  # TODO

    def test_build_two_transcript_deletion_pos(self):
        pass  # TODO

    def test_build_two_transcript_deletion_neg(self):
        pass  # TODO


class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.x = Exon(100, 199)  # C
        self.y = Exon(500, 599)  # G
        self.z = Exon(1200, 1299)  # T
        self.w = Exon(1500, 1599)  # C
        self.s = Exon(1700, 1799)  # G
        # introns: 99, 300, 600, 200, 100, ...
        reference_sequence = 'A' * 99 + 'C' * 100 + 'A' * 300 + 'G' * 100
        reference_sequence += 'A' * 600 + 'T' * 100 + 'A' * 200 + 'C' * 100
        reference_sequence += 'A' * 100 + 'G' * 100 + 'A' * 200
        self.reference_sequence = reference_sequence

    def test___init__(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        t = Transcript(gene=g, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(1, len(g.transcripts))
        self.assertEqual(g, t.gene)

        t = Transcript(gene=None, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(None, t.gene)
        t = Transcript(1, 10, genomic_start=1, genomic_end=10, domains=[Domain('name', [])])
        self.assertEqual(1, len(t.translations))

    def test___init__implicit_genomic_start(self):
        t = Transcript(gene=None, cds_start=1, cds_end=10, exons=[(1, 100), (200, 300), (400, 500)])
        self.assertEqual(1, t.start)
        self.assertEqual(t.start, t.genomic_start)
        self.assertEqual(500, t.end)
        self.assertEqual(t.end, t.genomic_end)
        self.assertEqual(1, t[0])
        self.assertEqual(500, t[1])
        self.assertFalse(Interval.overlaps((0, 0), t))
        self.assertTrue(Interval.overlaps((1, 1), t))
        self.assertTrue(Interval.overlaps((1, 50), t))

    def test___init__strand_mismatch(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)

        with self.assertRaises(AttributeError):
            t = Transcript(genomic_start=1, genomic_end=100, gene=g, strand=STRAND.NEG)


    def test___init__overlapping_exon_error(self):
        with self.assertRaises(AttributeError):
            Transcript(exons=[Exon(1, 15), Exon(10, 20)])

    def test_strand(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        t = Transcript(genomic_start=1, genomic_end=100, gene=g)
        self.assertEqual(STRAND.POS, t.strand)
        t = Transcript(genomic_start=1, genomic_end=100, strand=STRAND.POS)
        self.assertEqual(STRAND.POS, t.strand)

    def test_genomic_length(self):
        t = Transcript(genomic_start=1, genomic_end=2)
        self.assertEqual(2, t.genomic_length())

    def test_convert_cdna_to_genomic(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        self.assertEqual(1, len(t.translations))
        tl = t.translations[0]
        self.assertEqual(50, tl.start)
        self.assertEqual(249, tl.end)
        self.assertEqual(50, t.convert_cdna_to_genomic(tl.start))
        self.assertEqual(449, t.convert_cdna_to_genomic(tl.end))

        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        self.assertEqual(1, len(t.translations))
        tl = t.translations[0]
        self.assertEqual(50, tl.start)
        self.assertEqual(249, tl.end)
        self.assertEqual(450, t.convert_cdna_to_genomic(tl.start))
        self.assertEqual(51, t.convert_cdna_to_genomic(tl.end))

    def test_convert_aa_to_cdna(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        self.assertEqual(1, len(t.translations))
        tl = t.translations[0]
        self.assertEqual(50, tl.start)
        self.assertEqual(249, tl.end)
        self.assertEqual(Interval(56, 58), tl.convert_aa_to_cdna(3))

    def test_genomic_utr_regions(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        tl = t.translations[0]
        self.assertEqual([Interval(1, 50), Interval(449, 499)], tl.genomic_utr_regions())
    
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        tl = t.translations[0]
        self.assertEqual([Interval(450, 499), Interval(1, 51)], tl.genomic_utr_regions())
    
    def test_exon_number(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        for i, e in enumerate(t.exons):
            self.assertEqual(i + 1, t.exon_number(e))

        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        for i, e in enumerate(sorted(t.exons, key=lambda x: x.start, reverse=True)):
            self.assertEqual(i + 1, t.exon_number(e))
    
    def test_splicing_patterns(self):
        t = Transcript(1, 2, genomic_start=3, genomic_end=4)
        self.assertEqual(1, len(t.splicing_patterns()))

    def test_splicing_patterns_35(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        self.y.intact_end_splice = False
        self.z.intact_start_splice = False
        patterns = ft.splicing_patterns()
        self.assertEqual(1, len(patterns))
        self.assertEqual([self.x.end, self.y.start, self.z.end, self.w.start], patterns[0])

    def test_splicing_patterns_5_last_exon(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        self.w.intact_start_splice = False
        patterns = ft.splicing_patterns()
        self.assertEqual(1, len(patterns))
        self.assertEqual([self.x.end, self.y.start, self.y.end, self.z.start], patterns[0])

    def test_splicing_patterns_5(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        self.y.intact_start_splice = False
        patterns = ft.splicing_patterns()
        self.assertEqual(2, len(patterns))
        patterns = sorted(patterns)
        self.assertEqual([self.x.end, self.z.start, self.z.end, self.w.start], patterns[0])
        self.assertEqual([self.y.end, self.z.start, self.z.end, self.w.start], patterns[1])

    def test_splicing_patterns_3(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        self.y.intact_end_splice = False
        patterns = ft.splicing_patterns()
        self.assertEqual(2, len(patterns))
        patterns = sorted(patterns)
        self.assertEqual([self.x.end, self.z.start, self.z.end, self.w.start], patterns[1])
        self.assertEqual([self.x.end, self.y.start, self.z.end, self.w.start], patterns[0])
    
    def test_splicing_patterns_53(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        self.y.intact_end_splice = False
        self.y.intact_start_splice = False
        patterns = ft.splicing_patterns()
        self.assertEqual(1, len(patterns))
        self.assertEqual([self.x.end, self.z.start, self.z.end, self.w.start], patterns[0])

    def test_splicing_patterns_normal(self):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        ft = Transcript(exons=[self.x, self.y, self.z, self.w])
        patterns = ft.splicing_patterns()
        self.assertEqual(1, len(patterns))
        self.assertEqual([self.x.end, self.y.start, self.y.end, self.z.start, self.z.end, self.w.start], patterns[0])


class TestDomain(unittest.TestCase):

    def test_key(self):
        d = Domain('name', [])
        self.assertEqual(('name', None), d.key)

    def test___init__region_error(self):
        with self.assertRaises(AttributeError):
            Domain('name', [(1, 3), (4, 3)])

    def test___init__with_translation(self):
        t = Transcript(1, 2, genomic_start=3, genomic_end=4)
        for tl in t.translations:
            print(tl, tl.start, tl.end, tl.splicing_pattern)
        self.assertEqual(1, len(t.translations))
        t = t.translations[0]
        Domain('name', [], translation=t)
        self.assertEqual(1, len(t.domains))


class TestIntergenicRegion(unittest.TestCase):

    def test_key(self):
        b = IntergenicRegion('1', 1, 4, STRAND.NS)
        self.assertEqual(('1', 1, 4, STRAND.NS), b.key)


class TestBioInterval(unittest.TestCase):

    def test___eq__(self):
        a = BioInterval('test', 1, 2)
        b = BioInterval('test', 1, 2)
        c = BioInterval('test2', 1, 2)
        d = BioInterval('test', 3, 6)
        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertNotEqual(a, None)
        self.assertNotEqual(a, c)


class TestGene(unittest.TestCase):

    def test___hash__(self):
        g1 = Gene('test', 1, 2, 'name1', STRAND.POS)
        g2 = Gene('test', 1, 2, 'name2', STRAND.POS)
        h = set([g1, g2])
        self.assertEqual(2, len(h))

    def test___eq__(self):
        g1 = Gene('test', 1, 2, 'name1', STRAND.POS)
        g2 = Gene('test', 1, 2, 'name2', STRAND.POS)
        self.assertNotEqual(g1, g2)
        g3 = Gene('test2', 1, 2, 'name1', STRAND.POS)
        self.assertNotEqual(g1, g3)
        self.assertNotEqual(g3, g1)
        self.assertNotEqual(g1, None)
        self.assertNotEqual(None, g1)


class TestExon(unittest.TestCase):

    def test_end_splice_site(self):
        e = Exon(100, 199)
        self.assertEqual(2, SPLICE_SITE_RADIUS)
        self.assertEqual(Interval(198, 201), e.end_splice_site)

    def test_start_splice_site(self):
        e = Exon(100, 199)
        self.assertEqual(2, SPLICE_SITE_RADIUS)
        self.assertEqual(Interval(98, 101), e.start_splice_site)


class TestAnnotate(unittest.TestCase):
    def test_overlapping_transcripts(self):
        b = Breakpoint('C', 1000, strand=STRAND.POS)
        g = Gene('C', 1, 9999, 'gene1', STRAND.POS)
        t = Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=g)
        self.assertTrue(Interval.overlaps(b, t))
        Transcript(1, 10, exons=[(100, 199), (500, 699), (800, 900)], gene=g)
        h = Gene('C', 1, 9999, 'gene1', STRAND.NEG)
        Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=h, strand=STRAND.NEG)
        d = {'C': [g, h]}
        tlist = overlapping_transcripts(d, b)
        self.assertEqual(1, len(tlist))

    def test_gather_breakpoint_annotations_within_gene(self):
        b = Breakpoint('test', 150, 150)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(1, len(pos))
        self.assertEqual(1, len(neg))
        self.assertEqual(STRAND.POS, pos[0].strand)
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)
        self.assertEqual(STRAND.NEG, neg[0].strand)

    def test_gather_breakpoint_annotations_overlapping_gene(self):
        b = Breakpoint('test', 150, 230)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

        b = Breakpoint('test', 150, 225, strand=STRAND.POS)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(100, pos[0].start)
        self.assertEqual(200, pos[0].end)
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

        b = Breakpoint('test', 375, 425, strand=STRAND.POS)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(300, pos[0].start)
        self.assertEqual(400, pos[0].end)
        self.assertEqual(401, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

    def test_gather_breakpoint_annotations_overlapping_mutliple_genes_and_intergenic(self):
        b = Breakpoint('test', 150, 275)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(201, pos[1].start)
        self.assertEqual(b.end, pos[1].end)
        self.assertEqual(2, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(249, neg[0].end)

    def test_gather_breakpoint_annotations_overlapping_mutliple_pos_genes(self):
        b = Breakpoint('test', 575, 625)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(2, len(pos))
        self.assertEqual(1, len(neg))
        self.assertEqual(b.start, neg[0].start)
        self.assertEqual(b.end, neg[0].end)

    def test_gather_breakpoint_annotations_overlapping_mutliple_genes(self):
        b = Breakpoint('test', 300, 350)
        pos, neg = gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        self.assertEqual(1, len(pos))
        self.assertEqual(1, len(neg))

    def test_gather_annotations_intrachromosomal(self):
        b1 = Breakpoint('test', 150, 225, strand=STRAND.POS)
        b2 = Breakpoint('test', 375, 425, strand=STRAND.POS)
        bpp = BreakpointPair(b1, b2)
        ann_list = sorted(gather_annotations(REFERENCE_ANNOTATIONS, bpp),
                          key=lambda x: (x.break1, x.break2))
        self.assertEqual(5, len(ann_list))
        first = ann_list[0]
        print('original', bpp)
        for ann in ann_list:
            print(ann, [(g.start, g.end) for g in ann.encompassed_genes])
        self.assertEqual(1, len(first.encompassed_genes))
        self.assertEqual(0, len(first.genes_proximal_to_break1))
        self.assertEqual(1, len(first.genes_proximal_to_break2))
        self.assertEqual(0, len(first.genes_overlapping_break1))
        self.assertEqual(0, len(first.genes_overlapping_break2))
        near, dist = list(first.genes_proximal_to_break2)[0]
        self.assertEqual(50, dist)
        self.assertEqual(2, len(ann_list[1].encompassed_genes))

    def test_gather_annotations_interchromosomal(self):
        pass  # TODO

    def test_gather_annotations_intrachromosomal_within_gene_inversion(self):
        g = Gene('test', 1000, 3000, strand=STRAND.POS)
        t = Transcript(50, 450, gene=g, exons=[(1001, 1100), (1501, 1600), (2001, 2100), (2501, 2600)])
        ref = {
            'test': [g]
        }
        b1 = Breakpoint('test', 1250, strand=STRAND.POS)
        b2 = Breakpoint('test', 2250, strand=STRAND.NEG)
        bpp = BreakpointPair(b1, b2)
        ann_list = sorted(gather_annotations(ref, bpp),
                          key=lambda x: (x.break1, x.break2))
        print(ann_list)
        self.assertEqual(1, len(ann_list))
        self.assertEqual(ann_list[0].transcript1, ann_list[0].transcript2)
