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

    def test___init__(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        t = Transcript(gene=g, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(1, len(g.transcripts))
        self.assertEqual(g, t.gene)

        t = Transcript(gene=None, cds_start=1, cds_end=10, genomic_start=1, genomic_end=100)
        self.assertEqual(None, t.gene)
        t = Transcript(1, 10, 1, 10, domains=[Domain('name', [])])
        self.assertEqual(1, len(t.domains))

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

    def test___init__cds_error(self):
        with self.assertRaises(AttributeError):
            Transcript(9, 2, 3, 4)

        with self.assertRaises(AttributeError):
            Transcript('&', 2, 3, 4)

    def test___init__overlapping_exon_error(self):
        with self.assertRaises(AttributeError):
            Transcript(1, 10, exons=[Exon(1, 15), Exon(10, 20)])

    def test_strand(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        t = Transcript(genomic_start=1, genomic_end=100, gene=g)
        self.assertEqual(STRAND.POS, t.strand)
        t = Transcript(genomic_start=1, genomic_end=100, strand=STRAND.POS)
        self.assertEqual(STRAND.POS, t.strand)

    def test_genomic_length(self):
        t = Transcript(1, 2, 3, 4)
        self.assertEqual(2, t.genomic_length())

    def test_convert_cdna_to_genomic(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        self.assertEqual(50, t.cds_start)
        self.assertEqual(249, t.cds_end)
        self.assertEqual(50, t.convert_cdna_to_genomic(t.cds_start))
        self.assertEqual(449, t.convert_cdna_to_genomic(t.cds_end))

        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        self.assertEqual(50, t.cds_start)
        self.assertEqual(249, t.cds_end)
        self.assertEqual(450, t.convert_cdna_to_genomic(t.cds_start))
        self.assertEqual(51, t.convert_cdna_to_genomic(t.cds_end))

    def test_convert_aa_to_cdna(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        self.assertEqual(50, t.cds_start)
        self.assertEqual(249, t.cds_end)
        self.assertEqual(Interval(56, 58), t.convert_aa_to_cdna(3))

    def test_genomic_utr_regions(self):
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        self.assertEqual([Interval(1, 50), Interval(449, 499)], t.genomic_utr_regions())

        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        self.assertEqual([Interval(450, 499), Interval(1, 51)], t.genomic_utr_regions())


class TestDomain(unittest.TestCase):

    def test_key(self):
        d = Domain('name', [])
        self.assertEqual(('name', None), d.key)

    def test___init__region_error(self):
        with self.assertRaises(AttributeError):
            Domain('name', [(1, 3), (4, 3)])

    def test___init__with_transcript(self):
        t = Transcript(1, 2, 3, 4)
        Domain('name', [], transcript=t)
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
                          key=lambda x: (x.break1, x.break2))
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
