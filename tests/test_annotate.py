import unittest
from structural_variant.annotate import *
from structural_variant.constants import STRAND
from structural_variant.breakpoint import Breakpoint

REFERENCE_ANNOTATIONS = None


def setUpModule():
    global REFERENCE_ANNOTATIONS
    # REFERENCE_ANNOTATIONS = load_reference_genes(
    #     '/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv'
    # )


class TestTranscript(unittest.TestCase):

    def test_transcript_constructor_with_reference_object(self):
        g = Gene(chr='1', name='KRAS', strand='+')
        t = Transcript(gene=g, cds_start=1, cds_end=10)
        self.assertEqual(1, len(g.transcripts))
        self.assertEqual(g, t.gene)

    def test_transcript_constructor_without_reference_object(self):
        t = Transcript(gene=None, cds_start=1, cds_end=10)
        self.assertEqual(None, t.gene)


class TestAnnotate(unittest.TestCase):

    def test_overlapping_transcripts(self):
        b = Breakpoint('X', 1000, strand=STRAND.POS)
        g = Gene('X', 'gene1', STRAND.POS)
        Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=g)
        Transcript(1, 10, exons=[(100, 199), (500, 699), (800, 900)], gene=g)
        h = Gene('X', 'gene1', STRAND.NEG)
        Transcript(1, 10, exons=[(100, 199), (500, 699), (1200, 1300)], gene=h, strand=STRAND.NEG)
        d = {'X': [g, h]}
        tlist = overlapping_transcripts(d, b)
        self.assertEqual(1, len(tlist))
