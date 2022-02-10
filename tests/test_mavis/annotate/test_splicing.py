import argparse

import pytest
from mavis.annotate.constants import SPLICE_SITE_RADIUS
from mavis.annotate.genomic import Exon, PreTranscript
from mavis.annotate.splicing import predict_splice_sites
from mavis.annotate.variant import annotate_events
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import PROTOCOL, SPLICE_TYPE, STRAND, SVTYPE, reverse_complement
from mavis.interval import Interval

from ..mock import MockLongString, MockObject, get_example_genes

EXAMPLE_GENES = None


def setUpModule():
    global EXAMPLE_GENES
    EXAMPLE_GENES = get_example_genes()


@pytest.fixture
def neg_splicing_pattern():
    n = argparse.Namespace()
    n.ex1 = Exon(100, 199, strand=STRAND.NEG)  # C
    n.ex2 = Exon(500, 599, strand=STRAND.NEG)  # G
    n.ex3 = Exon(1200, 1299, strand=STRAND.NEG)  # T
    n.ex4 = Exon(1500, 1599, strand=STRAND.NEG)  # C
    n.ex5 = Exon(1700, 1799, strand=STRAND.NEG)  # G
    n.ex6 = Exon(2000, 2099, strand=STRAND.NEG)  # C
    # introns: 99, 300, 600, 200, 100, ...
    reference_sequence = 'a' * 99 + 'C' * 100 + 'a' * 300 + 'G' * 100
    reference_sequence += 'a' * 600 + 'T' * 100 + 'a' * 200 + 'C' * 100
    reference_sequence += 'a' * 100 + 'G' * 100 + 'a' * 200 + 'C' * 100
    n.reference_sequence = reference_sequence
    n.pre_transcript = PreTranscript(
        exons=[n.ex1, n.ex2, n.ex3, n.ex4, n.ex5, n.ex6], strand=STRAND.NEG
    )
    return n


@pytest.fixture
def pos_splicing_pattern():
    n = argparse.Namespace()
    n.ex1 = Exon(100, 199, strand=STRAND.POS)  # C
    n.ex2 = Exon(500, 599, strand=STRAND.POS)  # G
    n.ex3 = Exon(1200, 1299, strand=STRAND.POS)  # T
    n.ex4 = Exon(1500, 1599, strand=STRAND.POS)  # C
    n.ex5 = Exon(1700, 1799, strand=STRAND.POS)  # G
    n.ex6 = Exon(2000, 2099, strand=STRAND.POS)  # C
    # introns: 99, 300, 600, 200, 100, ...
    reference_sequence = 'a' * 99 + 'C' * 100 + 'a' * 300 + 'G' * 100
    reference_sequence += 'a' * 600 + 'T' * 100 + 'a' * 200 + 'C' * 100
    reference_sequence += 'a' * 100 + 'G' * 100 + 'a' * 200 + 'C' * 100
    n.reference_sequence = reference_sequence
    n.pre_transcript = PreTranscript(
        exons=[n.ex1, n.ex2, n.ex3, n.ex4, n.ex5, n.ex6], strand=STRAND.POS
    )
    return n


class TestSplicingPatterns:
    def test_single_exon(self):
        t = PreTranscript([(3, 4)], strand=STRAND.POS)
        patt = t.generate_splicing_patterns()
        assert len(patt) == 1
        assert len(patt[0]) == 0
        assert patt[0].splice_type == SPLICE_TYPE.NORMAL

    def test_normal_pattern_pos(self, pos_splicing_pattern):
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 1
        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex2.start,
            pos_splicing_pattern.ex2.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.NORMAL

    def test_normal_pattern_neg(self, neg_splicing_pattern):
        assert neg_splicing_pattern.pre_transcript.is_reverse
        patt = neg_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 1
        assert sorted([s.pos for s in patt[0]]) == [
            neg_splicing_pattern.ex1.end,
            neg_splicing_pattern.ex2.start,
            neg_splicing_pattern.ex2.end,
            neg_splicing_pattern.ex3.start,
            neg_splicing_pattern.ex3.end,
            neg_splicing_pattern.ex4.start,
            neg_splicing_pattern.ex4.end,
            neg_splicing_pattern.ex5.start,
            neg_splicing_pattern.ex5.end,
            neg_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.NORMAL

    def test_abrogate_a_pos(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.start_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 2

        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.SKIP

        assert [s.pos for s in patt[1]] == [
            pos_splicing_pattern.ex2.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[1].splice_type == SPLICE_TYPE.RETAIN

    def test_abrogate_a_neg(self, neg_splicing_pattern):
        neg_splicing_pattern.ex2.start_splice_site.intact = False
        patt = sorted(neg_splicing_pattern.pre_transcript.generate_splicing_patterns())
        assert len(patt) == 2
        assert sorted([s.pos for s in patt[0]]) == [
            neg_splicing_pattern.ex1.end,
            neg_splicing_pattern.ex3.start,
            neg_splicing_pattern.ex3.end,
            neg_splicing_pattern.ex4.start,
            neg_splicing_pattern.ex4.end,
            neg_splicing_pattern.ex5.start,
            neg_splicing_pattern.ex5.end,
            neg_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.SKIP
        assert sorted([s.pos for s in patt[1]]) == [
            neg_splicing_pattern.ex2.end,
            neg_splicing_pattern.ex3.start,
            neg_splicing_pattern.ex3.end,
            neg_splicing_pattern.ex4.start,
            neg_splicing_pattern.ex4.end,
            neg_splicing_pattern.ex5.start,
            neg_splicing_pattern.ex5.end,
            neg_splicing_pattern.ex6.start,
        ]
        assert patt[1].splice_type == SPLICE_TYPE.RETAIN

    def test_abrogate_a_last_exon(self, pos_splicing_pattern):
        pos_splicing_pattern.ex6.start_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 1
        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex2.start,
            pos_splicing_pattern.ex2.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.RETAIN

    def test_abrogate_d_first_exon(self, pos_splicing_pattern):
        pos_splicing_pattern.ex1.end_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 1
        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex2.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.RETAIN

    def test_abrogate_ad(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.start_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 2
        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.SKIP

        assert [s.pos for s in patt[1]] == [
            pos_splicing_pattern.ex2.end,
            pos_splicing_pattern.ex3.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[1].splice_type == SPLICE_TYPE.RETAIN

    def test_abrogate_da(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.end_splice_site.intact = False
        pos_splicing_pattern.ex3.start_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 1
        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex2.start,
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.RETAIN

    def test_multiple_exons_or_multiple_introns_abrogate_ada(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.start_splice_site.intact = False
        pos_splicing_pattern.ex2.end_splice_site.intact = False
        pos_splicing_pattern.ex3.start_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 2

        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.MULTI_SKIP

        assert [s.pos for s in patt[1]] == [
            pos_splicing_pattern.ex3.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[1].splice_type == SPLICE_TYPE.MULTI_RETAIN

    def test_multiple_exons_or_multiple_introns_abrogate_dad(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.end_splice_site.intact = False
        pos_splicing_pattern.ex3.start_splice_site.intact = False
        pos_splicing_pattern.ex3.end_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 2

        assert [s.pos for s in patt[0]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex2.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[0].splice_type == SPLICE_TYPE.MULTI_RETAIN

        assert [s.pos for s in patt[1]] == [
            pos_splicing_pattern.ex1.end,
            pos_splicing_pattern.ex4.start,
            pos_splicing_pattern.ex4.end,
            pos_splicing_pattern.ex5.start,
            pos_splicing_pattern.ex5.end,
            pos_splicing_pattern.ex6.start,
        ]
        assert patt[1].splice_type == SPLICE_TYPE.MULTI_SKIP

    def test_complex(self, pos_splicing_pattern):
        pos_splicing_pattern.ex2.end_splice_site.intact = False
        pos_splicing_pattern.ex4.end_splice_site.intact = False
        patt = pos_splicing_pattern.pre_transcript.generate_splicing_patterns()
        assert len(patt) == 4
        assert SPLICE_TYPE.COMPLEX in [p.splice_type for p in patt]


class TestExonSpliceSites:
    def test_end_splice_site(self):
        e = Exon(100, 199, strand=STRAND.POS)
        assert SPLICE_SITE_RADIUS == 2
        print(e.end_splice_site)
        assert Interval(198, 201) == e.end_splice_site

    def test_start_splice_site(self):
        e = Exon(100, 199, strand=STRAND.POS)
        assert SPLICE_SITE_RADIUS == 2
        print(e.start_splice_site)
        assert Interval(98, 101) == e.start_splice_site


class TestPredictSpliceSites:
    def test_gimap4(self):
        gimap4 = EXAMPLE_GENES['GIMAP4']
        donors = predict_splice_sites(gimap4.seq)
        for d in donors:
            print(d)
        assert len(donors) == 5

    def test_gimap4_reverse(self):
        gimap4 = EXAMPLE_GENES['GIMAP4']
        gimap4_seq = reverse_complement(gimap4.seq)
        donors = predict_splice_sites(gimap4_seq, True)
        for d in donors:
            assert gimap4_seq[d.start - 1 : d.end] == d.seq
        assert len(donors) == 5

    @pytest.mark.skip(reason='TODO: dependent functionality not yet implemented')
    def test_fusion_with_novel_splice_site(self):
        bpp = BreakpointPair(
            Breakpoint('7', 150268089, 150268089, 'L', '+'),
            Breakpoint('8', 79715940, 79715940, 'L', '-'),
            event_type=SVTYPE.ITRANS,
            protocol=PROTOCOL.GENOME,
            untemplated_seq='',
        )
        gimap4 = EXAMPLE_GENES['GIMAP4']
        il7 = EXAMPLE_GENES['IL7']
        ref_genome = {
            gimap4.chr: MockObject(seq=MockLongString(gimap4.seq, offset=gimap4.start - 1)),
            il7.chr: MockObject(seq=MockLongString(il7.seq, offset=il7.start - 1)),
        }
        annotations = annotate_events([bpp], {gimap4.chr: [gimap4], il7.chr: [il7]}, ref_genome)
        assert len(annotations) == 1
        ann = annotations[0]
        print(ann, ann.transcript1, ann.transcript2)
        print(ann.fusion)
        print(
            ann.fusion.transcripts[0].splicing_pattern,
            ann.fusion.transcripts[0].splicing_pattern.splice_type,
        )
        for ex in ann.fusion.transcripts[0].exons:
            print(ex, len(ex))
        assert False
