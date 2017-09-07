from mavis.annotate.constants import SPLICE_SITE_RADIUS, SPLICE_TYPE
from mavis.annotate.splicing import predict_splice_sites
from mavis.annotate.genomic import usTranscript, Exon
from mavis.annotate.variant import annotate_events
from mavis.annotate.file_io import load_reference_genome, load_annotations
from mavis.constants import STRAND, reverse_complement, SVTYPE, PROTOCOL
from mavis.interval import Interval
from mavis.breakpoint import Breakpoint, BreakpointPair
import unittest
import os
from . import DATA_DIR, MockLongString, MockSeq

REF_GENOME = {}
HUGO_GENES = {}
ANNOTATIONS = None


def setUpModule():
    global ANNOTATIONS
    temp = load_reference_genome(os.path.join(DATA_DIR, 'novel_exon_test_reference.fa'))
    ANNOTATIONS = load_annotations(os.path.join(DATA_DIR, 'novel_exon_test_annotations.tab'))
    for chr, gene_list in ANNOTATIONS.items():
        for gene in gene_list:
            assert(len(gene.aliases) == 1)
            hugo = gene.aliases[0].lower()
            HUGO_GENES[hugo] = gene
            offset = gene.start - 1
            if gene.chr in REF_GENOME:
                raise AssertionError('conflicting sequences')
            REF_GENOME[gene.chr] = MockSeq(MockLongString(str(temp[hugo].seq), offset=offset))


class TestSplicingPatterns(unittest.TestCase):

    def setUp(self):
        self.ex1 = Exon(100, 199, strand=STRAND.POS)  # C
        self.ex2 = Exon(500, 599, strand=STRAND.POS)  # G
        self.ex3 = Exon(1200, 1299, strand=STRAND.POS)  # T
        self.ex4 = Exon(1500, 1599, strand=STRAND.POS)  # C
        self.ex5 = Exon(1700, 1799, strand=STRAND.POS)  # G
        self.ex6 = Exon(2000, 2099, strand=STRAND.POS)  # C
        # introns: 99, 300, 600, 200, 100, ...
        reference_sequence = 'A' * 99 + 'C' * 100 + 'A' * 300 + 'G' * 100
        reference_sequence += 'A' * 600 + 'T' * 100 + 'A' * 200 + 'C' * 100
        reference_sequence += 'A' * 100 + 'G' * 100 + 'A' * 200 + 'C' * 100
        self.reference_sequence = reference_sequence
        self.ust = usTranscript(exons=[self.ex1, self.ex2, self.ex3, self.ex4, self.ex5, self.ex6], strand=STRAND.POS)

    def test_single_exon(self):
        t = usTranscript([(3, 4)], strand=STRAND.POS)
        patt = t.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(0, len(patt[0]))
        self.assertEqual(SPLICE_TYPE.NORMAL, patt[0].splice_type)

    def test_normal_pattern(self):
        for strand in [STRAND.POS, STRAND.NEG]:
            self.ust.strand = strand
            patt = self.ust.generate_splicing_patterns()
            self.assertEqual(1, len(patt))
            self.assertEqual(
                [
                    self.ex1.end, self.ex2.start,
                    self.ex2.end, self.ex3.start,
                    self.ex3.end, self.ex4.start,
                    self.ex4.end, self.ex5.start,
                    self.ex5.end, self.ex6.start
                ],
                patt[0]
            )
            self.assertEqual(SPLICE_TYPE.NORMAL, patt[0].splice_type)

    def test_abrogate_A_pos(self):
        self.ex2.start_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(2, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[1]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_A_neg(self):
        self.ex1 = Exon(100, 199, strand=STRAND.NEG)  # C
        self.ex2 = Exon(500, 599, strand=STRAND.NEG)  # G
        self.ex3 = Exon(1200, 1299, strand=STRAND.NEG)  # T
        self.ex4 = Exon(1500, 1599, strand=STRAND.NEG)  # C
        self.ex5 = Exon(1700, 1799, strand=STRAND.NEG)  # G
        self.ex6 = Exon(2000, 2099, strand=STRAND.NEG)  # C
        self.ex2.start_splice_site.intact = False
        self.ust = usTranscript(exons=[self.ex1, self.ex2, self.ex3, self.ex4, self.ex5, self.ex6], strand=STRAND.NEG)
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(2, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)
        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[1]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_A_last_exon(self):
        self.ex6.start_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_abrogate_D_first_exon(self):
        self.ex1.end_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_abrogate_AD(self):
        self.ex2.start_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(2, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex2.end, self.ex3.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[1]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[1].splice_type)

    def test_abrogate_DA(self):
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(1, len(patt))
        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.RETAIN, patt[0].splice_type)

    def test_multiple_exons_or_multiple_introns_abrogate_ADA(self):
        self.ex2.start_splice_site.intact = False
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(2, len(patt))

        self.assertEqual(
            [
                self.ex1.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_SKIP, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex3.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[1]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_RETAIN, patt[1].splice_type)

    def test_multiple_exons_or_multiple_introns_abrogate_DAD(self):
        self.ex2.end_splice_site.intact = False
        self.ex3.start_splice_site.intact = False
        self.ex3.end_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
        self.assertEqual(2, len(patt))

        self.assertEqual(
            [
                self.ex1.end, self.ex2.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[0]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_RETAIN, patt[0].splice_type)

        self.assertEqual(
            [
                self.ex1.end, self.ex4.start,
                self.ex4.end, self.ex5.start,
                self.ex5.end, self.ex6.start
            ],
            patt[1]
        )
        self.assertEqual(SPLICE_TYPE.MULTI_SKIP, patt[1].splice_type)

    def test_complex(self):
        self.ex2.end_splice_site.intact = False
        self.ex4.end_splice_site.intact = False
        patt = self.ust.generate_splicing_patterns()
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
        gimap4 = HUGO_GENES['gimap4']
        gimap4_seq = gimap4.get_seq(REF_GENOME)
        print(gimap4_seq[:100])
        donors = predict_splice_sites(gimap4_seq)
        for d in donors:
            print(d)
        self.assertEqual(5, len(donors))

    def test_gimap4_reverse(self):
        gimap4 = HUGO_GENES['gimap4']
        gimap4_seq = gimap4.get_seq(REF_GENOME)
        gimap4_seq = reverse_complement(gimap4_seq)
        print(gimap4_seq[:100])
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
        annotations = annotate_events([bpp], ANNOTATIONS, REF_GENOME)
        self.assertEqual(1, len(annotations))
        ann = annotations[0]
        print(ann, ann.transcript1, ann.transcript2)
        print(ann.fusion)
        print(ann.fusion.transcripts[0].splicing_pattern, ann.fusion.transcripts[0].splicing_pattern.splice_type)
        for ex in ann.fusion.transcripts[0].exons:
            print(ex, len(ex))
        self.assertTrue(False)
