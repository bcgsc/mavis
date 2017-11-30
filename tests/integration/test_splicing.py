import os
import unittest

from mavis.annotate.constants import SPLICE_SITE_RADIUS, SPLICE_TYPE
from mavis.annotate.file_io import load_annotations, load_reference_genome
from mavis.annotate.genomic import Exon, PreTranscript
from mavis.annotate.splicing import predict_splice_sites
from mavis.annotate.variant import annotate_events
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import PROTOCOL, reverse_complement, STRAND, SVTYPE
from mavis.interval import Interval

from . import DATA_DIR, MockLongString, MockObject, get_example_genes

EXAMPLE_GENES = None


def setUpModule():
    global EXAMPLE_GENES
    EXAMPLE_GENES = get_example_genes()


class TestSplicingPatterns(unittest.TestCase):

    def setUp(self):
        self.setup_by_strand(STRAND.POS)

    def setup_by_strand(self, strand):
        self.ex1 = Exon(100, 199, strand=strand)  # C
        self.ex2 = Exon(500, 599, strand=strand)  # G
        self.ex3 = Exon(1200, 1299, strand=strand)  # T
        self.ex4 = Exon(1500, 1599, strand=strand)  # C
        self.ex5 = Exon(1700, 1799, strand=strand)  # G
        self.ex6 = Exon(2000, 2099, strand=strand)  # C
        # introns: 99, 300, 600, 200, 100, ...
        reference_sequence = 'a' * 99 + 'C' * 100 + 'a' * 300 + 'G' * 100
        reference_sequence += 'a' * 600 + 'T' * 100 + 'a' * 200 + 'C' * 100
        reference_sequence += 'a' * 100 + 'G' * 100 + 'a' * 200 + 'C' * 100
        self.reference_sequence = reference_sequence
        self.pre_transcript = PreTranscript(exons=[self.ex1, self.ex2, self.ex3, self.ex4, self.ex5, self.ex6], strand=strand)

    def test_single_exon(self):
        t = PreTranscript([(3, 4)], strand=STRAND.POS)
        patt = t.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(0, len(patt[0]))
        self.assertEqual(SPLICE_TYPE.NORMAL, patt[0].splice_type)

    def test_normal_pattern_pos(self):
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.NORMAL, patt[0].splice_type)

    def test_normal_pattern_neg(self):
        self.setup_by_strand(STRAND.NEG)
        self.assertTrue(self.pre_transcript.is_reverse)
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            sorted([s.pos for s in patt[0]])
        )
        self.assertEqual(SPLICE_TYPE.NORMAL, patt[0].splice_type)

    def test_abrogate_a_pos(self):
        self.ex2.start_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(2, len(patt))

        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[1]]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_a_neg(self):
        self.setup_by_strand(STRAND.NEG)
        self.ex2.start_splice_site.intact = False
        patt = sorted(self.pre_transcript.generate_splicing_patterns())
        self.assertEqual(2, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            sorted([s.pos for s in patt[0]])
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)
        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            sorted([s.pos for s in patt[1]])
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_a_last_exon(self):
        self.ex6.start_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_abrogate_d_first_exon(self):
        self.ex1.end_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_abrogate_ad(self):
        self.ex2.start_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(2, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[1]]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_da(self):
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_multiple_exons_or_multiple_introns_abrogate_ada(self):
        self.ex2.start_splice_site.intact = False
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(2, len(patt))

        self.assertEqual(
            [
                self.ex1.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[1]]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_RETAIN, patt[1].splice_type)

    def test_multiple_exons_or_multiple_introns_abrogate_dad(self):
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        self.ex3.end_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(2, len(patt))

        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[0]]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_RETAIN, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex1.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            [s.pos for s in patt[1]]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_SKIP, patt[1].splice_type)

    def test_complex(self):
        self.ex2.end_splice_site.intact = False
        self.ex4.end_splice_site.intact = False
        patt = self.pre_transcript.generate_splicing_patterns()
        self.assertEqual(4, len(patt))
        self.assertTrue(SPLICE_TYPE.COMPLEX in [p.splice_type for p in patt])


class TestExonSpliceSites(unittest.TestCase):

    def test_end_splice_site(self):
        e = Exon(100, 199, strand=STRAND.POS)
        self.assertEqual(2, SPLICE_SITE_RADIUS)
        self.assertEqual(Interval(198, 201), e.end_splice_site)

    def test_start_splice_site(self):
        e = Exon(100, 199, strand=STRAND.POS)
        self.assertEqual(2, SPLICE_SITE_RADIUS)
        self.assertEqual(Interval(98, 101), e.start_splice_site)


class TestPredictSpliceSites(unittest.TestCase):

    def test_gimap4(self):
        gimap4 = EXAMPLE_GENES['GIMAP4']
        donors = predict_splice_sites(gimap4.seq)
        for d in donors:
            print(d)
        self.assertEqual(5, len(donors))

    def test_gimap4_reverse(self):
        gimap4 = EXAMPLE_GENES['GIMAP4']
        gimap4_seq = reverse_complement(gimap4.seq)
        donors = predict_splice_sites(gimap4_seq, True)
        for d in donors:
            self.assertEqual(d.seq, gimap4_seq[d.start - 1:d.end])
        self.assertEqual(5, len(donors))

    def test_fusion_with_novel_splice_site(self):
        raise unittest.SkipTest('TODO: dependent functionality not yet implemented')
        bpp = BreakpointPair(
            Breakpoint('7', 150268089, 150268089, 'L', '+'),
            Breakpoint('8', 79715940, 79715940, 'L', '-'),
            event_type=SVTYPE.ITRANS,
            protocol=PROTOCOL.GENOME,
            untemplated_seq=''
        )
        gimap4 = EXAMPLE_GENES['GIMAP4']
        il7 = EXAMPLE_GENES['IL7']
        ref_genome = {
            gimap4.chr: MockObject(seq=MockLongString(gimap4.seq, offset=gimap4.start - 1)),
            il7.chr: MockObject(seq=MockLongString(il7.seq, offset=il7.start - 1))
        }
        annotations = annotate_events([bpp], {gimap4.chr: [gimap4], il7.chr: [il7]}, ref_genome)
        self.assertEqual(1, len(annotations))
        ann = annotations[0]
        print(ann, ann.transcript1, ann.transcript2)
        print(ann.fusion)
        print(ann.fusion.transcripts[0].splicing_pattern, ann.fusion.transcripts[0].splicing_pattern.splice_type)
        for ex in ann.fusion.transcripts[0].exons:
            print(ex, len(ex))
        self.assertTrue(False)
