import unittest
from structural_variant.annotate import *
from structural_variant.constants import STRAND
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from tests import REFERENCE_ANNOTATIONS_FILE


REFERENCE_ANNOTATIONS = None


def setUpModule():
    global REFERENCE_ANNOTATIONS
    REFERENCE_ANNOTATIONS = load_reference_genes(REFERENCE_ANNOTATIONS_FILE)
    print('loaded {} annotations', sum([len(l) for l in REFERENCE_ANNOTATIONS.values()]))


class TestTranscript(unittest.TestCase):

    def test_transcript_constructor_with_reference_object(self):
        g = Gene('1', 1, 9999, name='KRAS', strand='+')
        t = Transcript(gene=g, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(1, len(g.transcripts))
        self.assertEqual(g, t.gene)

    def test_transcript_constructor_without_reference_object(self):
        t = Transcript(gene=None, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(None, t.gene)

    def test_transcript_constructor_implicit_genomic_start(self):
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


class TestAnnotate(unittest.TestCase):
    def test_overlapping_transcripts(self):
        b = Breakpoint('X', 1000, strand=STRAND.POS)
        g = Gene('X', 1, 9999, 'gene1', STRAND.POS)
        t = Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=g)
        self.assertTrue(Interval.overlaps(b, t))
        Transcript(1, 10, exons=[(100, 199), (500, 699), (800, 900)], gene=g)
        h = Gene('X', 1, 9999, 'gene1', STRAND.NEG)
        Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=h, strand=STRAND.NEG)
        d = {'X': [g, h]}
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
                          key=lambda x: (x.breakpoint_pair.break1, x.breakpoint_pair.break2))
        self.assertEqual(5, len(ann_list))
        first = ann_list[0]
        self.assertEqual(1, len(first.encompassed_genes))
        self.assertEqual(0, len(first.nearest_gene_break1))
        self.assertEqual(1, len(first.nearest_gene_break2))
        self.assertEqual(0, len(first.genes_at_break1))
        self.assertEqual(0, len(first.genes_at_break2))
        near, dist = list(first.nearest_gene_break2)[0]
        self.assertEqual(50, dist)
        self.assertEqual(2, len(ann_list[1].encompassed_genes))

    def test_gather_annotations_interchromosomal(self):
        pass
