import unittest
from structural_variant.annotate import *
from structural_variant.constants import STRAND
from structural_variant.breakpoint import Breakpoint
from tests import REFERENCE_ANNOTATIONS_FILE


REFERENCE_ANNOTATIONS = None


def setUpModule():
    global REFERENCE_ANNOTATIONS
    REFERENCE_ANNOTATIONS = load_reference_genes(REFERENCE_ANNOTATIONS_FILE)


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
