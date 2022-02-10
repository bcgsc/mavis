import argparse
import unittest

import pytest
from mavis.annotate.base import BioInterval, ReferenceName
from mavis.annotate.file_io import load_annotations, load_reference_genome
from mavis.annotate.fusion import FusionTranscript, determine_prime
from mavis.annotate.genomic import Exon, Gene, PreTranscript, Template, Transcript
from mavis.annotate.protein import Domain, DomainRegion, Translation, calculate_orf, translate
from mavis.annotate.variant import (
    Annotation,
    _gather_annotations,
    _gather_breakpoint_annotations,
    annotate_events,
    overlapping_transcripts,
)
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, PRIME, PROTOCOL, STRAND, SVTYPE, reverse_complement
from mavis.error import NotSpecifiedError
from mavis.interval import Interval

from ...util import get_data
from ..mock import MockObject, get_example_genes

REFERENCE_ANNOTATIONS = None
REFERENCE_GENOME = None
REF_CHR = 'fake'
ALT_REF_CHR = 'ref2'


def setUpModule():
    global REFERENCE_ANNOTATIONS, REFERENCE_GENOME, REF_CHR, EXAMPLE_GENES
    EXAMPLE_GENES = get_example_genes()
    REFERENCE_ANNOTATIONS = load_annotations(get_data('mock_reference_annotations2.json'))
    count = sum([len(genes) for genes in REFERENCE_ANNOTATIONS.values()])
    print('loaded annotations', count)
    assert count >= 6  # make sure this is the file we expect
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    assert REF_CHR in REFERENCE_GENOME
    print('loaded the reference genome', get_data('mock_reference_genome.fa'))


class TestTemplate:
    def test_template_hashing(self):
        t = Template('1', 1, 10)
        d = {'1': 1, '2': 2, 1: '5'}
        assert t.name == '1'
        assert d[t.name] == 1
        assert d[t] == 1


@pytest.fixture
def intervals():
    n = argparse.Namespace()
    n.x = Interval(100, 199)  # C
    n.y = Interval(500, 599)  # G
    n.z = Interval(1200, 1299)  # T
    n.w = Interval(1500, 1599)  # C
    n.s = Interval(1700, 1799)  # G
    # introns: 99, 300, 600, 200, 100, ...
    reference_sequence = 'A' * 99 + 'C' * 100 + 'A' * 300 + 'G' * 100
    reference_sequence += 'A' * 600 + 'T' * 100 + 'A' * 200 + 'C' * 100
    reference_sequence += 'A' * 100 + 'G' * 100 + 'A' * 200 + 'T' * 100

    n.a = Interval(2000, 2099)  # T
    n.b = Interval(2600, 2699)  # C
    n.c = Interval(3000, 3099)  # G
    n.d = Interval(3300, 3399)  # T
    reference_sequence += 'A' * 500 + 'C' * 100 + 'A' * 300 + 'G' * 100
    reference_sequence += 'A' * 200 + 'T' * 100 + 'A' * 200
    n.reference_sequence = reference_sequence

    n.b1 = Interval(600, 699)  # A
    n.b2 = Interval(800, 899)  # G
    n.b3 = Interval(1100, 1199)  # T
    n.b4 = Interval(1400, 1499)  # A
    n.b5 = Interval(1700, 1799)  # G
    n.b6 = Interval(2100, 2199)  # A
    alternate_sequence = 'C' * 599 + 'A' * 100 + 'C' * 100 + 'G' * 100
    alternate_sequence += 'C' * 200 + 'T' * 100 + 'C' * 200 + 'A' * 100
    alternate_sequence += 'C' * 200 + 'G' * 100 + 'C' * 300 + 'A' * 100
    alternate_sequence += 'C' * 200
    n.alternate_sequence = alternate_sequence
    return n


class TestFusionTranscript:
    def test__pull_exons_left_pos_intronic(self, intervals):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 700, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (700 - 600 + 1)
        )
        assert seq == expt
        assert len(new_exons) == 2
        e = new_exons[0][0]
        assert e.start == 1
        assert e.end == 100
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True

    def test__pull_exons_left_pos_intronic_splice(self, intervals):
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 201, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'C' * 100 + 'A' * 2
        assert seq == expt
        assert len(new_exons) == 1
        e = new_exons[0][0]
        assert e.start == 1
        assert e.end == 100
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is False

    def test__pull_exons_left_pos_exonic(self, intervals):
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        print('transcriptt exons:', t.exons)
        b = Breakpoint(REF_CHR, 199, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'C' * 100
        assert seq == expt
        assert len(new_exons) == 1
        e = new_exons[0][0]
        assert e.start == 1
        assert e.end == 100
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is False

    def test__pull_exons_left_pos_exonic_splice(self, intervals):
        # 100-199, 500-599, 1200-1299, 1500-1599, 1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 101, orient=ORIENT.LEFT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'C' * 2
        assert seq == expt
        assert len(new_exons) == 1
        e = new_exons[0][0]
        assert e.start == 1
        assert e.end == 2
        assert e.start_splice_site.intact is False
        assert e.end_splice_site.intact is False

    def test__pull_exons_right_pos_intronic(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 1600, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'A' * (1699 - 1600 + 1) + 'G' * len(intervals.s)
        assert seq == expt
        assert len(new_exons) == 1

        b = Breakpoint(REF_CHR, 300, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'A' * (499 - 300 + 1) + 'G' * 100 + 'A' * (1199 - 600 + 1) + 'T' * 100
        expt += 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100

        assert seq == expt
        assert len(new_exons) == 4
        e = new_exons[0][0]
        assert e.start == 201
        assert e.end == 300
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True

    def test__pull_exons_right_pos_intronic_splice(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = (
            'AA'
            + 'T' * 100
            + 'A' * (1499 - 1300 + 1)
            + 'C' * 100
            + 'A' * (1699 - 1600 + 1)
            + 'G' * 100
        )
        assert seq == expt
        assert len(new_exons) == 3
        e = new_exons[0][0]
        assert e.start_splice_site.intact is False
        assert e.end_splice_site.intact is True

    def test__pull_exons_right_pos_exonic(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 1201, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'T' * 99 + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        assert seq == expt
        assert len(new_exons) == 3
        e = new_exons[0][0]
        assert e.start_splice_site.intact is False
        assert e.end_splice_site.intact is True

    def test__pull_exons_right_pos_exonic_splice(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b = Breakpoint(REF_CHR, 1298, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'TT' + 'A' * (1499 - 1300 + 1) + 'C' * 100 + 'A' * (1699 - 1600 + 1) + 'G' * 100
        assert seq == expt
        assert len(new_exons) == 3
        e = new_exons[0][0]
        assert e.start_splice_site.intact is False
        assert e.end_splice_site.intact is False

    def test__pull_exons_right_neg_intronic(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        b = Breakpoint(REF_CHR, 700, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = 'A' * (1199 - 700 + 1) + 'T' * 100 + 'A' * (1499 - 1300 + 1) + 'C' * 100
        expt += 'A' * (1699 - 1600 + 1) + 'G' * 100
        expt = reverse_complement(expt)
        assert seq == expt
        assert len(new_exons) == 3
        e = new_exons[0][0]
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True
        assert e.start == 1
        assert e.end == 100
        assert seq[e.start - 1 : e.end] == 'C' * 100
        e = new_exons[1][0]
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True
        assert e.start == 201
        assert e.end == 300
        assert seq[e.start - 1 : e.end] == 'G' * 100

    def test__pull_exons_right_neg_intronic_splice(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        b = Breakpoint(REF_CHR, 1198, orient=ORIENT.RIGHT)
        seq, new_exons = FusionTranscript._pull_exons(t, b, intervals.reference_sequence)
        expt = (
            'AA'
            + 'T' * 100
            + 'A' * (1499 - 1300 + 1)
            + 'C' * 100
            + 'A' * (1699 - 1600 + 1)
            + 'G' * 100
        )
        expt = reverse_complement(expt)
        assert seq == expt
        assert len(new_exons) == 3
        e = new_exons[0][0]
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True
        assert e.start == 1
        assert e.end == 100
        assert seq[e.start - 1 : e.end] == 'C' * 100
        e = new_exons[1][0]
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is True
        assert e.start == 201
        assert e.end == 300
        assert seq[e.start - 1 : e.end] == 'G' * 100
        e = new_exons[2][0]
        assert e.start_splice_site.intact is True
        assert e.end_splice_site.intact is False
        assert e.start == 501
        assert e.end == 600
        assert seq[e.start - 1 : e.end] == 'A' * 100

    def test_build_single_transcript_indel(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 599, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'ATCGATCG'
            + 'T' * len(intervals.z)
        )
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )

        assert ft.seq == expt
        assert len(ft.exons) == 5

        for i, ex in enumerate(t.exons):
            n = ft.exons[i]
            assert ft.exon_mapping[n.position] == ex

        assert ft.exons[0].start == 1
        assert ft.exons[0].end == 100

        splice_pattern = [(True, True), (True, False), (False, True), (True, True), (True, True)]
        char_pattern = [x * 100 for x in ['C', 'G', 'T', 'C', 'G']]

        for i in range(0, len(splice_pattern)):
            s, t = splice_pattern[i]
            ex = ft.exons[i]
            assert ex.start_splice_site.intact == s
            assert ex.end_splice_site.intact == t
            assert ft.seq[ex.start - 1 : ex.end] == char_pattern[i]

    def test_build_single_transcript_inversion(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)
        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += 'ATCGTC' + 'A' * len(intervals.z)
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )
        exons = [(1, 100), (401, 500), (1407, 1506), (1607, 1706)]
        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i][0]
            assert ft.exons[i].end == exons[i][1]
        assert ft.seq == expt
        assert len(ft.exons) == 4

    def test_build_single_transcript_inversion_transcriptome(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.TRANS
        )
        ft = FusionTranscript.build(ann, ref)
        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += 'ATCGTC' + 'A' * len(intervals.z)
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )
        exons = [
            Exon(1, 100, strand=STRAND.POS),
            Exon(401, 500, intact_end_splice=False, strand=STRAND.POS),
            Exon(501, 1406, intact_start_splice=False, intact_end_splice=False, strand=STRAND.POS),
            Exon(1407, 1506, intact_start_splice=False, strand=STRAND.POS),
            Exon(1607, 1706, strand=STRAND.POS),
        ]
        print(ft.exons)
        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i].start
            assert ft.exons[i].end == exons[i].end
            assert ft.exons[i].start_splice_site.intact == exons[i].start_splice_site.intact
            assert ft.exons[i].end_splice_site.intact == exons[i].end_splice_site.intact
        assert ft.seq == expt
        assert len(ft.exons) == 5

    def test_build_single_transcript_inversion_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        b1 = Breakpoint(REF_CHR, 1300, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGTC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.s)
            + 'T' * (1699 - 1600 + 1)
            + 'G' * len(intervals.w)
            + 'T' * (1499 - 1300 + 1)
        )
        expt += 'T' * len(intervals.z) + 'GACGAT' + 'T' * (1199 - 600 + 1) + 'C' * len(intervals.y)
        expt += 'T' * (499 - 200 + 1) + 'G' * len(intervals.x)

        exons = [(1, 100), (201, 300), (1207, 1306), (1607, 1706)]

        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i][0]
            assert ft.exons[i].end == exons[i][1]
        assert ft.seq == expt
        assert len(ft.exons) == 4

    def test_build_single_transcript_duplication_pos(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)
        assert ft.get_strand() == STRAND.POS

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += 'T' * len(intervals.z) + 'ATCGATCG' + 'T' * len(intervals.z)
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )
        assert ft.seq == expt
        exons = [(1, 100), (401, 500), (1101, 1200), (1209, 1308), (1509, 1608), (1709, 1808)]
        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i][0]
            assert ft.exons[i].end == exons[i][1]

        assert len(ft.exons) == 6
        assert ft.exons[2].start_splice_site.intact
        assert ft.exons[3].end_splice_site.intact
        assert not ft.exons[2].end_splice_site.intact
        assert not ft.exons[3].start_splice_site.intact

    def test_build_single_transcript_duplication_pos_transcriptome(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.TRANS
        )
        ft = FusionTranscript.build(ann, ref)
        assert ft.get_strand() == STRAND.POS

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += 'T' * len(intervals.z) + 'ATCGATCG' + 'T' * len(intervals.z)
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )
        assert ft.seq == expt
        exons = [
            Exon(1, 100, strand=STRAND.POS),
            Exon(401, 500, strand=STRAND.POS),
            Exon(1101, 1200, intact_end_splice=False, strand=STRAND.POS),
            Exon(1201, 1208, intact_start_splice=False, intact_end_splice=False, strand=STRAND.POS),
            Exon(1209, 1308, intact_start_splice=False, strand=STRAND.POS),
            Exon(1509, 1608, strand=STRAND.POS),
            Exon(1709, 1808, strand=STRAND.POS),
        ]
        print(ft.exons)
        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i].start
            assert ft.exons[i].end == exons[i].end
            assert ft.exons[i].start_splice_site.intact == exons[i].start_splice_site.intact
            assert ft.exons[i].end_splice_site.intact == exons[i].end_splice_site.intact

        assert len(ft.exons) == 7

    def test_build_single_transcript_duplication_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        t = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 1299, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGATCG')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += 'T' * len(intervals.z) + 'ATCGATCG' + 'T' * len(intervals.z)
        expt += (
            'A' * (1499 - 1300 + 1)
            + 'C' * len(intervals.w)
            + 'A' * (1699 - 1600 + 1)
            + 'G' * len(intervals.s)
        )
        expt = reverse_complement(expt)
        assert ft.seq == expt

        exons = [(1, 100), (201, 300), (501, 600), (609, 708), (1309, 1408), (1709, 1808)]

        for i in range(len(exons)):
            assert ft.exons[i].start == exons[i][0]
            assert ft.exons[i].end == exons[i][1]

        assert len(ft.exons) == 6
        assert ft.exons[2].start_splice_site.intact
        assert ft.exons[3].end_splice_site.intact
        assert not ft.exons[2].end_splice_site.intact
        assert not ft.exons[3].start_splice_site.intact
        assert ft.exon_number(ft.exons[2]) == 3
        assert ft.exon_number(ft.exons[3]) == 3

    def test_build_two_transcript_inversion_5prime_pos(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.NEG
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGACTC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)
        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += (
            'ATCGACTC' + 'G' * len(intervals.b) + 'T' * (2599 - 2100 + 1) + 'A' * len(intervals.a)
        )
        assert ft.seq == expt
        assert len(ft.exons) == 4
        assert ft.exons[3].end_splice_site.intact
        assert not ft.exons[2].start_splice_site.intact
        assert ft.exons[2].end_splice_site.intact
        assert ft.exon_number(ft.exons[1]) == 2
        assert ft.exon_number(ft.exons[2]) == 3

    def test_build_two_transcript_inversion_5prime_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.POS
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='ATCGACTC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.INV, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)
        expt = (
            'T' * len(intervals.a) + 'A' * (2599 - 2100 + 1) + 'C' * len(intervals.b) + 'ATCGACTC'
        )
        expt += (
            'T' * (1199 - 600 + 1)
            + 'C' * len(intervals.y)
            + 'T' * (499 - 200 + 1)
            + 'G' * len(intervals.x)
        )

        assert len(ft.exons) == 4
        assert ft.exon_number(ft.exons[1]) == 2
        assert ft.exon_number(ft.exons[2]) == 4
        assert ft.seq == expt

    def test_build_two_transcript_duplication_pos(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.POS
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGAC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)
        expt = 'T' * len(intervals.a) + 'A' * (2599 - 2100 + 1) + 'C' * len(intervals.b) + 'ATCGAC'
        expt += 'T' * len(intervals.z) + 'A' * (1499 - 1300 + 1) + 'C' * len(intervals.w)
        expt += 'A' * (1699 - 1600 + 1) + 'G' * len(intervals.s)

        assert len(ft.exons) == 5
        assert ft.exon_number(ft.exons[1]) == 2
        assert ft.exon_number(ft.exons[2]) == 3
        assert ft.seq == expt

    def test_build_two_transcript_duplication_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.NEG
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='ATCGAC')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.s)
            + 'T' * (1699 - 1600 + 1)
            + 'G' * len(intervals.w)
            + 'T' * (1499 - 1300 + 1)
        )
        expt += 'A' * len(intervals.z) + 'GTCGAT' + 'G' * len(intervals.b) + 'T' * (2599 - 2100 + 1)
        expt += 'A' * len(intervals.a)

        assert len(ft.exons) == 5
        assert ft.exon_number(ft.exons[1]) == 2
        assert ft.exon_number(ft.exons[2]) == 3
        assert ft.seq == expt

    def test_build_two_transcript_deletion_pos(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.POS
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2700, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGTGT')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
            + 'AACGTGT'
        )
        expt += (
            'A' * (2999 - 2700 + 1)
            + 'G' * len(intervals.c)
            + 'A' * (3299 - 3100 + 1)
            + 'T' * len(intervals.d)
        )

        assert ft.seq == expt
        assert 4, len(ft.exons)

    def test_build_two_transcript_deletion_pos_transcriptome(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.POS
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2700, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGTGT')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.TRANS
        )
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
            + 'AACGTGT'
        )
        expt += (
            'A' * (2999 - 2700 + 1)
            + 'G' * len(intervals.c)
            + 'A' * (3299 - 3100 + 1)
            + 'T' * len(intervals.d)
        )

        assert ft.seq == expt
        assert 5, len(ft.exons)

    def test_build_two_transcript_deletion_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # a:2000-2099, b:2600-2699, c:3000-3099, d:3300-3399
        #   TTTTTTTTT    CCCCCCCCC    GGGGGGGGG    TTTTTTTTT
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )

        t2 = PreTranscript(
            exons=[intervals.a, intervals.b, intervals.c, intervals.d], strand=STRAND.NEG
        )
        print('t1 exons', t1.exons)
        print('t2 exons', t2.exons)
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2699, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='AACGAGTGT')
        ref = {REF_CHR: MockObject(seq=intervals.reference_sequence)}
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME
        )

        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.s)
            + 'T' * (1699 - 1600 + 1)
            + 'G' * len(intervals.w)
            + 'T' * (1499 - 1300 + 1)
        )
        expt += (
            'A' * len(intervals.z) + 'ACACTCGTT' + 'G' * len(intervals.b) + 'T' * (2599 - 2100 + 1)
        )
        expt += 'A' * len(intervals.a)

        assert ft.seq == expt
        assert 5, len(ft.exons)
        assert ft.exon_number(ft.exons[2]) == 3
        assert ft.exon_number(ft.exons[3]) == 3

    def test_build_two_transcript_translocation(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.POS,
        )
        t2 = PreTranscript(
            exons=[
                intervals.b1,
                intervals.b2,
                intervals.b3,
                intervals.b4,
                intervals.b5,
                intervals.b6,
            ],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1199, orient=ORIENT.LEFT)
        b2 = Breakpoint('ref2', 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='GCAACATAT')
        ref = {
            REF_CHR: MockObject(seq=intervals.reference_sequence),
            'ref2': MockObject(seq=intervals.alternate_sequence),
        }
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME
        )
        assert ann.break1 == b1
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.x)
            + 'A' * (499 - 200 + 1)
            + 'G' * len(intervals.y)
            + 'A' * (1199 - 600 + 1)
        )
        expt += (
            'GCAACATAT'
            + 'C' * (1399 - 1200 + 1)
            + 'A' * len(intervals.b4)
            + 'C' * (1699 - 1500 + 1)
        )
        expt += 'G' * len(intervals.b5) + 'C' * (2099 - 1800 + 1) + 'A' * len(intervals.b6)

        assert ft.seq == expt
        assert 5, len(ft.exons)
        assert 2, ft.exon_number(ft.exons[1])
        assert 4, ft.exon_number(ft.exons[2])

    def test_build_two_transcript_translocation_neg(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        t2 = PreTranscript(
            exons=[
                intervals.b1,
                intervals.b2,
                intervals.b3,
                intervals.b4,
                intervals.b5,
                intervals.b6,
            ],
            strand=STRAND.NEG,
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(ALT_REF_CHR, 1199, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='TCTACATAT')
        ref = {
            REF_CHR: MockObject(seq=intervals.reference_sequence),
            ALT_REF_CHR: MockObject(seq=intervals.alternate_sequence),
        }
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME
        )
        assert ann.break1 == b1
        assert ann.break2 == b2
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.s)
            + 'T' * (1699 - 1600 + 1)
            + 'G' * len(intervals.w)
            + 'T' * (1499 - 1300 + 1)
        )
        expt += (
            'A' * len(intervals.z) + 'ATATGTAGA' + 'A' * len(intervals.b3) + 'G' * (1099 - 900 + 1)
        )
        expt += 'C' * len(intervals.b2) + 'G' * (799 - 700 + 1) + 'T' * len(intervals.b1)

        assert ft.seq == expt
        assert len(ft.exons) == 6
        assert 3, ft.exon_number(ft.exons[2])
        assert 3, ft.exon_number(ft.exons[3])

    def test_build_two_transcript_inverted_translocation(self, intervals):
        # x:100-199, y:500-599, z:1200-1299, w:1500-1599, s:1700-1799
        #   CCCCCCC    GGGGGGG    TTTTTTTTT    CCCCCCCCC    GGGGGGGGG
        # 1:600-699, 2:800-899, 3:1100-1199, 4:1400-1499, 5:1700-1799 6:2100-2199
        #   AAAAAAA    GGGGGGG,   TTTTTTTTT,   AAAAAAAAA,   GGGGGGGGG   AAAAAAAAA
        t1 = PreTranscript(
            exons=[intervals.x, intervals.y, intervals.z, intervals.w, intervals.s],
            strand=STRAND.NEG,
        )
        t2 = PreTranscript(
            exons=[
                intervals.b1,
                intervals.b2,
                intervals.b3,
                intervals.b4,
                intervals.b5,
                intervals.b6,
            ],
            strand=STRAND.POS,
        )
        b1 = Breakpoint(REF_CHR, 1200, orient=ORIENT.RIGHT)
        b2 = Breakpoint(ALT_REF_CHR, 1200, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='GATACATAT')
        ref = {
            REF_CHR: MockObject(seq=intervals.reference_sequence),
            ALT_REF_CHR: MockObject(seq=intervals.alternate_sequence),
        }
        ann = Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.TRANS, protocol=PROTOCOL.GENOME
        )
        assert ann.break1 == b1
        assert ann.break2 == b2
        ft = FusionTranscript.build(ann, ref)

        expt = (
            'C' * len(intervals.s)
            + 'T' * (1699 - 1600 + 1)
            + 'G' * len(intervals.w)
            + 'T' * (1499 - 1300 + 1)
        )
        expt += (
            'A' * len(intervals.z) + 'ATATGTATC' + 'C' * (1399 - 1200 + 1) + 'A' * len(intervals.b4)
        )
        expt += (
            'C' * (1699 - 1500 + 1)
            + 'G' * len(intervals.b5)
            + 'C' * (2099 - 1800 + 1)
            + 'A' * len(intervals.b6)
        )

        assert ft.seq == expt
        assert len(ft.exons) == 6
        assert 3, ft.exon_number(ft.exons[2])
        assert 4, ft.exon_number(ft.exons[3])


@pytest.fixture
def mock_ann_obj():
    n = argparse.Namespace()
    n.gene = Gene(REF_CHR, 1, 900, strand=STRAND.POS)

    n.pre_transcript = PreTranscript(
        exons=[(101, 200), (301, 400), (501, 600), (701, 800)], gene=n.gene
    )
    n.gene.transcripts.append(n.pre_transcript)

    n.transcript = Transcript(n.pre_transcript, n.pre_transcript.generate_splicing_patterns()[0])
    n.pre_transcript.transcripts.append(n.transcript)

    n.translation = Translation(51, 350, n.transcript)
    n.transcript.translations.append(n.translation)

    n.spliced_seq = (
        'GGTGAATTTCTAGTTTGCCTTTTCAGCTAGGGATTAGCTTTTTAGGGGTCCCAATG'
        'CCTAGGGAGATTTCTAGGTCCTCTGTTCCTTGCTGACCTCCAATAATCAGAAAATGCTGTGAAGGAAAAAC'
        'AAAATGAAATTGCATTGTTTCTACCGGCCCTTTATCAAGCCCTGGCCACCATGATAGTCATGAATTCCAAT'
        'TGTGTTGAAATCACTTCAATGTGTTTCTCTTCTTTCTGGGAGCTTACACACTCAAGTTCTGGATGCTTTGA'
        'TTGCTATCAGAAGCCGTTAAATAGCTACTTATAAATAGCATTGAGTTATCAGTACTTTCATGTCTTGATAC'
        'ATTTCTTCTTGAAAATGTTCATGCTTGCTGATTTGTCTGTTTGTTGAGAGGAGAATGTTC'
    )

    n.domain = Domain(name=REF_CHR, regions=[(11, 20), (51, 60)], translation=n.translation)
    n.translation.domains.append(n.domain)
    return n


class TestSequenceFetching:
    def test_fetch_gene_seq_from_ref(self, mock_ann_obj):
        expt = str(REFERENCE_GENOME[REF_CHR][0:900].seq).upper()
        assert mock_ann_obj.gene.get_seq(REFERENCE_GENOME) == expt
        # gene seq should be the same if gene in on reverse strand b/c gene seq always given on pos
        mock_ann_obj.gene.strand = STRAND.NEG
        assert mock_ann_obj.gene.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_gene_seq_from_stored(self, mock_ann_obj):
        expt = 'AAA'
        mock_ann_obj.gene.seq = expt
        assert mock_ann_obj.gene.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_gene_seq_force_uncached(self, mock_ann_obj):
        expt = str(REFERENCE_GENOME[REF_CHR][0:900].seq).upper()
        mock_ann_obj.gene.seq = 'AAA'
        assert mock_ann_obj.gene.get_seq(REFERENCE_GENOME, ignore_cache=True) == expt

    def test_fetch_us_transcript_seq_from_ref(self, mock_ann_obj):
        expt = str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper()
        assert mock_ann_obj.pre_transcript.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_us_transcript_seq_from_ref_revcomp(self, mock_ann_obj):
        mock_ann_obj.gene.strand = STRAND.NEG
        expt = reverse_complement(str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper())
        assert mock_ann_obj.pre_transcript.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_us_transcript_seq_from_stored(self, mock_ann_obj):
        expt = 'AAA'
        mock_ann_obj.pre_transcript.seq = expt
        assert mock_ann_obj.pre_transcript.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_us_transcript_seq_from_parent_gene(self, mock_ann_obj):
        mock_ann_obj.gene.seq = 'A' * len(mock_ann_obj.gene)
        assert mock_ann_obj.pre_transcript.get_seq() == 'A' * len(mock_ann_obj.pre_transcript)

    def test_fetch_us_transcript_seq_from_parent_gene_revcomp(self, mock_ann_obj):
        mock_ann_obj.gene.seq = 'A' * len(mock_ann_obj.gene)
        mock_ann_obj.gene.strand = STRAND.NEG
        assert mock_ann_obj.pre_transcript.get_seq() == 'T' * len(mock_ann_obj.pre_transcript)

    def test_fetch_us_transcript_seq_force_uncached(self, mock_ann_obj):
        expt = str(REFERENCE_GENOME[REF_CHR][100:800].seq).upper()
        mock_ann_obj.pre_transcript.seq = 'AAA'
        assert mock_ann_obj.pre_transcript.get_seq(REFERENCE_GENOME, ignore_cache=True) == expt

    def test_fetch_transcript_seq_from_ref(self, mock_ann_obj):
        assert mock_ann_obj.transcript.get_seq(REFERENCE_GENOME) == mock_ann_obj.spliced_seq

    def test_fetch_transcript_seq_from_ref_revcomp(self, mock_ann_obj):
        mock_ann_obj.gene.strand = STRAND.NEG
        assert mock_ann_obj.transcript.get_seq(REFERENCE_GENOME) == reverse_complement(
            mock_ann_obj.spliced_seq
        )

    def test_fetch_transcript_seq_from_stored(self, mock_ann_obj):
        expt = 'AAA'
        mock_ann_obj.transcript.seq = expt
        assert mock_ann_obj.transcript.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_transcript_seq_from_parent_ust(self, mock_ann_obj):
        mock_ann_obj.pre_transcript.seq = 'A' * len(mock_ann_obj.pre_transcript)
        assert mock_ann_obj.transcript.get_seq() == 'A' * len(mock_ann_obj.transcript)

    def test_fetch_transcript_seq_from_parent_gene(self, mock_ann_obj):
        mock_ann_obj.gene.seq = 'A' * len(mock_ann_obj.gene)
        assert mock_ann_obj.transcript.get_seq() == 'A' * len(mock_ann_obj.transcript)

    def test_fetch_transcript_seq_force_uncached(self, mock_ann_obj):
        mock_ann_obj.transcript.seq = 'AAA'
        assert (
            mock_ann_obj.transcript.get_seq(REFERENCE_GENOME, ignore_cache=True)
            == mock_ann_obj.spliced_seq
        )

    def test_fetch_translation_aa_seq_from_ref(self, mock_ann_obj):
        cds = mock_ann_obj.spliced_seq[
            mock_ann_obj.translation.start - 1 : mock_ann_obj.translation.end
        ]
        assert mock_ann_obj.translation.get_aa_seq(REFERENCE_GENOME) == translate(cds)

    def test_fetch_translation_cds_seq_from_ref(self, mock_ann_obj):
        cds = mock_ann_obj.spliced_seq[
            mock_ann_obj.translation.start - 1 : mock_ann_obj.translation.end
        ]
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == cds

    def test_fetch_translation_cds_seq_from_ref_revcomp(self, mock_ann_obj):
        mock_ann_obj.gene.strand = STRAND.NEG
        cdna = reverse_complement(mock_ann_obj.spliced_seq)
        cds = cdna[mock_ann_obj.translation.start - 1 : mock_ann_obj.translation.end]
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == cds

    def test_fetch_translation_cds_seq_from_stored(self, mock_ann_obj):
        expt = 'AAA'
        mock_ann_obj.translation.seq = expt
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == expt

    def test_fetch_translation_cds_seq_from_parent_transcript(self, mock_ann_obj):
        mock_ann_obj.transcript.seq = 'A' * len(mock_ann_obj.transcript)
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == 'A' * len(
            mock_ann_obj.translation
        )

    def test_fetch_translation_cds_seq_from_parent_ust(self, mock_ann_obj):
        mock_ann_obj.pre_transcript.seq = 'A' * len(mock_ann_obj.pre_transcript)
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == 'A' * len(
            mock_ann_obj.translation
        )

    def test_fetch_translation_cds_seq_from_parent_gene(self, mock_ann_obj):
        mock_ann_obj.gene.seq = 'A' * len(mock_ann_obj.gene)
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME) == 'A' * len(
            mock_ann_obj.translation
        )

    def test_fetch_translation_cds_seq_force_uncached(self, mock_ann_obj):
        mock_ann_obj.translation.seq = 'AAA'
        cds = mock_ann_obj.spliced_seq[
            mock_ann_obj.translation.start - 1 : mock_ann_obj.translation.end
        ]
        assert mock_ann_obj.translation.get_seq(REFERENCE_GENOME, ignore_cache=True) == cds

    def test_fetch_domain_seq_from_ref(self, mock_ann_obj):
        seqs = ['VPC*PPIIRK', 'C*NHFNVFLF']
        assert mock_ann_obj.domain.get_seqs(REFERENCE_GENOME) == seqs


@pytest.fixture
def unstranded_gene():
    gene = Gene('1', 1, 500, strand=STRAND.POS)
    pre_transcript = PreTranscript(gene=gene, exons=[(1, 100), (200, 300), (400, 500)])
    gene.unspliced_transcripts.append(pre_transcript)
    for spl in pre_transcript.generate_splicing_patterns():
        t = Transcript(pre_transcript, spl)
        pre_transcript.spliced_transcripts.append(t)
        tl = Translation(51, 250, t)
        t.translations.append(tl)
    return gene


class TestStrandInheritance:
    def test_strand_gene(self, unstranded_gene):
        assert unstranded_gene.get_strand() == STRAND.POS

    def test_strand_us_transcript(self, unstranded_gene):
        assert unstranded_gene.unspliced_transcripts[0].get_strand() == STRAND.POS

    def test_strand_spl_transcript(self, unstranded_gene):
        assert unstranded_gene.spliced_transcripts[0].get_strand() == STRAND.POS

    def test_strand_translation(self, unstranded_gene):
        assert unstranded_gene.spliced_transcripts[0].translations[0].get_strand() == STRAND.POS


@pytest.fixture
def coord_conv_setup():
    n = argparse.Namespace()
    n.gene = Gene('1', 15, 700, strand=STRAND.POS)

    n.pre_transcript = PreTranscript(gene=n.gene, exons=[(101, 200), (301, 400), (501, 600)])
    n.gene.unspliced_transcripts.append(n.pre_transcript)
    assert 1 == len(n.pre_transcript.generate_splicing_patterns())

    spl = n.pre_transcript.generate_splicing_patterns()[0]
    n.transcript = Transcript(n.pre_transcript, spl)
    n.pre_transcript.spliced_transcripts.append(n.transcript)

    n.translation = Translation(51, 251, n.transcript)
    n.transcript.translations.append(n.translation)

    n.rev_gene = Gene('1', 15, 700, strand=STRAND.NEG)
    n.rev_ust = PreTranscript(gene=n.rev_gene, exons=[(101, 200), (301, 400), (501, 600)])
    n.gene.unspliced_transcripts.append(n.rev_ust)
    assert 1 == len(n.rev_ust.generate_splicing_patterns())

    spl = n.rev_ust.generate_splicing_patterns()[0]
    n.rev_transcript = Transcript(n.rev_ust, spl)
    n.rev_ust.spliced_transcripts.append(n.rev_transcript)

    n.rev_translation = Translation(51, 251, n.rev_transcript)
    n.rev_transcript.translations.append(n.rev_translation)
    return n


class TestCoordinateCoversion:
    def test_cdna_to_genomic(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_cdna_to_genomic(50) == 150
        assert coord_conv_setup.transcript.convert_cdna_to_genomic(250) == 550

    def test_cdna_to_genomic_before(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_cdna_to_genomic(-1) == 100
        assert coord_conv_setup.transcript.convert_cdna_to_genomic(-50) == 51

    def test_cdna_to_genomic_after(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_cdna_to_genomic(350) == 650

    def test_cdna_to_genomic_revcomp(self, coord_conv_setup):
        assert coord_conv_setup.rev_transcript.convert_cdna_to_genomic(50) == 551
        assert coord_conv_setup.rev_transcript.convert_cdna_to_genomic(250) == 151

    def test_genomic_to_cdna(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_cdna(150) == 50
        assert coord_conv_setup.transcript.convert_genomic_to_cdna(549) == 249

    def test_genomic_to_cdna_before(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(100) == (1, -1)

    def test_genomic_to_cdna_after(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(601) == (300, 1)

    def test_genomic_to_cdna_revcomp(self, coord_conv_setup):
        assert coord_conv_setup.rev_transcript.convert_genomic_to_cdna(551) == 50
        assert coord_conv_setup.rev_transcript.convert_genomic_to_cdna(151) == 250

    def test_aa_to_cdna(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_aa_to_cdna(1) == Interval(51, 53)
        assert coord_conv_setup.translation.convert_aa_to_cdna(67) == Interval(249, 251)

    def test_cdna_to_aa(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_cdna_to_aa(51) == 1
        assert coord_conv_setup.translation.convert_cdna_to_aa(251) == 67
        with pytest.raises(IndexError):
            coord_conv_setup.translation.convert_cdna_to_aa(50)
        with pytest.raises(IndexError):
            coord_conv_setup.translation.convert_cdna_to_aa(252)

    def test_genomic_to_cds(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds(151) == 1
        assert coord_conv_setup.translation.convert_genomic_to_cds(551) == 201

    def test_genomic_to_cds_3prime_utr(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds(150) == -1

    def test_genomic_to_cds_5prime_utr(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds(552) == 202

    def test_genomic_to_cds_notation(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(151) == '1'
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(551) == '201'

    def test_genomic_to_cds_notation_3prime_utr(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(150) == '-1'

    def test_genomic_to_cds_notation_5prime_utr(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(552) == '*1'

    def test_genomic_to_cds_notation_intronic_pos(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(202) == '50+2'

    def test_genomic_to_cds_notation_intronic_neg(self, coord_conv_setup):
        assert coord_conv_setup.translation.convert_genomic_to_cds_notation(299) == '51-2'

    def test_genomic_to_nearest_cdna_exonic(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(101) == (1, 0)
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(600) == (300, 0)
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(301) == (101, 0)

    def test_genomic_to_nearest_cdna_intronic_pos(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(210) == (100, 10)

    def test_genomic_to_nearest_cdna_intronic_neg(self, coord_conv_setup):
        assert coord_conv_setup.transcript.convert_genomic_to_nearest_cdna(299) == (101, -2)

    def test_genomic_to_nearest_cdna_rev_exonic(self, coord_conv_setup):
        assert coord_conv_setup.rev_transcript.convert_genomic_to_nearest_cdna(101) == (300, 0)
        assert coord_conv_setup.rev_transcript.convert_genomic_to_nearest_cdna(600) == (1, 0)
        assert coord_conv_setup.rev_transcript.convert_genomic_to_nearest_cdna(400) == (101, 0)

    def test_genomic_to_nearest_cdna_rev_intronic_pos(self, coord_conv_setup):
        assert coord_conv_setup.rev_transcript.convert_genomic_to_nearest_cdna(210) == (201, -10)

    def test_genomic_to_nearest_cdna_rev_intronic_neg(self, coord_conv_setup):
        assert coord_conv_setup.rev_transcript.convert_genomic_to_nearest_cdna(299) == (200, 2)


class TestUSTranscript:
    def test___init__implicit_start(self):
        t = PreTranscript(gene=None, exons=[(1, 100), (200, 300), (400, 500)], strand=STRAND.POS)
        assert t.start == 1
        assert t.start == t.start
        assert t.end == 500
        assert t.end == t.end
        assert t[0] == 1
        assert t[1] == 500
        assert not Interval.overlaps((0, 0), t)
        assert Interval.overlaps((1, 1), t)
        assert Interval.overlaps((1, 50), t)

    def test___init__strand_mismatch(self):
        g = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)

        with pytest.raises(AssertionError):
            PreTranscript([(1, 100)], gene=g, strand=STRAND.NEG)

    def test___init__overlapping_exon_error(self):
        with pytest.raises(AttributeError):
            PreTranscript(exons=[Exon(1, 15), Exon(10, 20)])

    def test_exon_number(self):
        t = PreTranscript(gene=None, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)
        for i, e in enumerate(t.exons):
            assert t.exon_number(e) == i + 1

        t = PreTranscript(gene=None, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.NEG)
        for i, e in enumerate(sorted(t.exons, key=lambda x: x.start, reverse=True)):
            assert t.exon_number(e) == i + 1


class TestDomain:
    def test___init__region_error(self):
        with pytest.raises(AttributeError):
            Domain('name', [(1, 3), (4, 3)])

    def test_get_seq_from_ref(self):
        ref = {'1': MockObject(seq='CCCTAATCCCCTTT')}
        g = Gene('1', 1, 16, strand=STRAND.NEG)
        t = PreTranscript(exons=[(2, 5), (7, 15)], gene=g)
        tl = Translation(4, 11, t, [])
        d = Domain('name', [(1, 2)], translation=tl)
        assert d.get_seqs(ref) == [translate('GGGGAT')]

    def test_get_seq_from_translation_seq(self):
        t = PreTranscript(exons=[(2, 5), (7, 15)], seq='CCCTAATCCCCTTT', strand=STRAND.NEG)
        tl = Translation(4, 11, t, [])
        d = Domain('name', [(1, 2)], translation=tl)
        assert d.get_seqs() == [translate('TAATCC')]

    def test_align_seq(self):
        regions = [
            DomainRegion(216, 261, 'DVNECITGSHSCRLGESCINTVGSFRCQRDSSCGTGYELTEDNSCK'),
            DomainRegion(262, 307, 'DIDECESGIHNCLPDFICQNTLGSFRCRPKLQCKSGFIQDALGNCI'),
            DomainRegion(308, 355, 'DINECLSISAPCPIGHTCINTEGSYTCQKNVPNCGRGYHLNEEGTRCV'),
            DomainRegion(356, 398, 'DVDECAPPAEPCGKGHRCVNSPGSFRCECKTGYYFDGISRMCV'),
            DomainRegion(399, 440, 'DVNECQRYPGRLCGHKCENTLGSYLCSCSVGFRLSVDGRSCE'),
            DomainRegion(441, 480, 'DINECSSSPCSQECANVYGSYQCYCRRGYQLSDVDGVTCE'),
            DomainRegion(481, 524, 'DIDECALPTGGHICSYRCINIPGSFQCSCPSSGYRLAPNGRNCQ'),
            DomainRegion(525, 578, 'DIDECVTGIHNCSINETCFNIQGGFRCLAFECPENYRRSAATLQQEKTDTVRCI'),
        ]
        print(regions)
        refseq = (
            'MERAAPSRRVPLPLLLLGGLALLAAGVDADVLLEACCADGHRMATHQKDCSLPYATESKE'
            'CRMVQEQCCHSQLEELHCATGISLANEQDRCATPHGDNASLEATFVKRCCHCCLLGRAAQ'
            'AQGQSCEYSLMVGYQCGQVFQACCVKSQETGDLDVGGLQETDKIIEVEEEQEDPYLNDRC'
            'RGGGPCKQQCRDTGDEVVCSCFVGYQLLSDGVSCEDVNECITGSHSCRLGESCINTVGSF'
            'RCQRDSSCGTGYELTEDNSCKDIDECESGIHNCLPDFICQNTLGSFRCRPKLQCKSGFIQ'
            'DALGNCIDINECLSISAPCPIGHTCINTEGSYTCQKNVPNCGRGYHLNEEGTRCVDVDEC'
            'APPAEPCGKGHRCVNSPGSFRCECKTGYYFDGISRMCVDVNECQRYPGRLCGHKCENTLG'
            'SYLCSCSVGFRLSVDGRSCEDINECSSSPCSQECANVYGSYQCYCRRGYQLSDVDGVTCE'
            'DIDECALPTGGHICSYRCINIPGSFQCSCPSSGYRLAPNGRNCQDIDECVTGIHNCSINE'
            'TCFNIQGGFRCLAFECPENYRRSAATLQQEKTDTVRCIKSCRPNDVTCVFDPVHTISHTV'
            'ISLPTFREFTRPEEIIFLRAITPPHPASQANIIFDITEGNLRDSFDIIKRYMDGMTVGVV'
            'RQVRPIVGPFHAVLKLEMNYVVGGVVSHRNVVNVHIFVSEYWF'
        )

        d = Domain('name', regions)
        assert len(refseq) >= 578
        match, total, temp = d.align_seq(refseq)
        assert total == sum([len(d.seq) for d in regions])
        assert match == total
        assert len(temp) == len(regions)
        for dr1, dr2 in zip(temp, regions):
            assert dr2.start == dr1.start
            assert dr2.end == dr1.end
            assert dr2.seq == dr1.seq

        refseq = (
            'MHRPPRHMGNKAMEPMDSPLMSAIPRLRPLQPMGRPPMQLLMDSLPLVILLQLPPRHTASLSRGMALVLMIPPL'
            'LQSPPPRPPMQLSLHMALSLLIQPMGSSQQPLHLQAIPLHSRLVMIRAVTLSRTPMGNRAAMDSRVAMVNKAAMGSSLP'
            'LVTHPKLDPTAKLQVNIANRAAATGSRTLLMTQSEEELGAIT*'
        )

        d = (
            'QAAAQQGYSAYTAQPTQGYAQTTQAYGQQSYGTYGQPTDVSYTQAQTTATYGQTAYATSYGQPPTGYTTPTAPQAYSQP'
            'VQGYGTGAYDTTTATVTTTQASYAAQSAYGTQPAYPAYGQQPAATAPTSYSSTQPTSYDQSSYSQQNTYGQPSSYGQQS'
            'SYGQQSSYGQQPPTSYPPQTGSYSQAPSQYSQQSSSYGQQSSFRQ'
        )

        dom = Domain('name', [DomainRegion(1, len(d), d)])
        with pytest.raises(UserWarning):
            dom.align_seq(refseq)

        seq = (
            'MASTDYSTYSQAAAQQGYSAYTAQPTQGYAQTTQAYGQQSYGTYGQPTDVSYTQAQTTATYGQTAYATSYGQPPTGY'
            'TTPTAPQAYSQPVQGYGTGAYDTTTATVTTTQASYAAQSAYGTQPAYPAYGQQPAATAPTRPQDGNKPTETSQPQSSTG'
            'GYNQPSLGYGQSNYSYPQVPGSYPMQPVTAPPSYPPTSYSSTQPTSYDQSSYSQQNTYGQPSSYGQQSSYGQQSSYGQQ'
            'PPTSYPPQTGSYSQAPSQYSQQSSSYGQQNPSYDSVRRGAWGNNMNSGLNKSPPLGGAQTISKNTEQRPQPDPYQILGP'
            'TSSRLANPGSGQIQLWQFLLELLSDSANASCITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYYDKNI'
            'MTKVHGKRYAYKFDFHGIAQALQPHPTESSMYKYPSDISYMPSYHAHQQKVNFVPPHPSSMPVTSSSFFGAASQYWTSP'
            'TGGIYPNPNVPRHPNTHVPSHLGSYY'
        )
        d = 'IYVQGLNDSVTLDDLADFFKQCGVVKMNKRTGQPMIHIYLDKETGKPKGDATVSYEDPPTAKAAVEWFDGKDFQGSKLK'
        dom = Domain('name', [DomainRegion(1, len(d), d)])

        with pytest.raises(UserWarning):
            m, t, regions = dom.align_seq(seq)


class TestBioInterval:
    def test___eq__(self):
        a = BioInterval(REF_CHR, 1, 2)
        b = BioInterval(REF_CHR, 1, 2)
        c = BioInterval('test2', 1, 2)
        assert a == a
        assert b == a
        assert a is not None
        assert c != a


class TestGene:
    def test___hash__(self):
        g1 = Gene(REF_CHR, 1, 2, 'name1', STRAND.POS)
        g2 = Gene(REF_CHR, 1, 2, 'name2', STRAND.POS)
        h = set([g1, g2])
        assert len(h) == 2

    def test___eq__(self):
        g1 = Gene(REF_CHR, 1, 2, 'name1', STRAND.POS)
        g2 = Gene(REF_CHR, 1, 2, 'name2', STRAND.POS)
        assert g2 != g1
        g3 = Gene('test2', 1, 2, 'name1', STRAND.POS)
        assert g3 != g1
        assert g1 != g3
        assert None != g1  # noqa: E711
        assert g1 != None  # noqa: E711

    def test_get_seq(self):
        ref = {'1': MockObject(seq='AACCCTTTGGG')}
        g = Gene('1', 3, 8, strand=STRAND.POS)
        assert g.get_seq(ref) == 'CCCTTT'
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
            'TCCTTCTTCCTGTTATTTTGTTTCCAATTTAGTCTTT'
        ).upper()
        assert g.get_seq(REFERENCE_GENOME) == seq


class TestAnnotationGathering:
    def test_overlapping_transcripts(self):
        b = Breakpoint('C', 1000, strand=STRAND.POS)
        g = Gene('C', 1, 9999, 'gene1', STRAND.POS)
        t = PreTranscript(exons=[(100, 199), (500, 699), (1200, 1300)], gene=g)
        g.transcripts.append(t)
        assert Interval.overlaps(b, t)
        t = PreTranscript(exons=[(100, 199), (500, 699), (800, 900)], gene=g)
        g.transcripts.append(t)
        h = Gene('C', 1, 9999, 'gene1', STRAND.NEG)
        t = PreTranscript(exons=[(100, 199), (500, 699), (1200, 1300)], gene=h, strand=STRAND.NEG)
        h.transcripts.append(t)
        d = {'C': [g, h]}
        tlist = overlapping_transcripts(d, b)
        assert len(tlist) == 1

    def test_breakpoint_within_gene(self):
        b = Breakpoint(REF_CHR, 150, 150)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 1
        assert len(neg) == 1
        assert pos[0].get_strand() == STRAND.POS
        assert neg[0].start == b.start
        assert neg[0].end == b.end
        assert neg[0].get_strand() == STRAND.NEG

    def test_breakpoint_overlapping_gene(self):
        b = Breakpoint(REF_CHR, 150, 230)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 2
        assert pos[1].start == 201
        assert pos[1].end == b.end
        assert len(neg) == 1
        assert neg[0].start == b.start
        assert neg[0].end == b.end

        b = Breakpoint(REF_CHR, 150, 225, strand=STRAND.POS)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 2
        assert pos[0].start == 100
        assert pos[0].end == 200
        assert pos[1].start == 201
        assert pos[1].end == b.end
        assert len(neg) == 1
        assert neg[0].start == b.start
        assert neg[0].end == b.end

        b = Breakpoint(REF_CHR, 375, 425, strand=STRAND.POS)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 2
        assert pos[0].start == 300
        assert pos[0].end == 400
        assert pos[1].start == 401
        assert pos[1].end == b.end
        assert len(neg) == 1
        assert neg[0].start == b.start
        assert neg[0].end == b.end

    def test_breakpoint_overlapping_mutliple_genes_and_intergenic(self):
        b = Breakpoint(REF_CHR, 150, 275)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 2
        assert pos[1].start == 201
        assert pos[1].end == b.end
        assert len(neg) == 2
        assert neg[0].start == b.start
        assert neg[0].end == 249

    def test_breakpoint_overlapping_mutliple_pos_genes(self):
        b = Breakpoint(REF_CHR, 575, 625)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 2
        assert len(neg) == 1
        assert neg[0].start == b.start
        assert neg[0].end == b.end

    def test_breakpoint_overlapping_mutliple_genes(self):
        b = Breakpoint(REF_CHR, 300, 350)
        pos, neg = _gather_breakpoint_annotations(REFERENCE_ANNOTATIONS, b)
        assert len(pos) == 1
        assert len(neg) == 1

    def test_intrachromosomal(self):
        b1 = Breakpoint(REF_CHR, 150, 225, strand=STRAND.POS)
        b2 = Breakpoint(REF_CHR, 375, 425, strand=STRAND.POS)
        bpp = BreakpointPair(b1, b2, protocol=PROTOCOL.GENOME, event_type=SVTYPE.DEL)
        ann_list = sorted(
            _gather_annotations(REFERENCE_ANNOTATIONS, bpp), key=lambda x: (x.break1, x.break2)
        )
        assert len(ann_list) == 5
        first = ann_list[0]
        assert len(first.encompassed_genes) == 1
        assert len(first.genes_proximal_to_break1) == 0
        assert len(first.genes_proximal_to_break2) == 1
        assert len(first.genes_overlapping_break1) == 0
        assert len(first.genes_overlapping_break2) == 0
        near, dist = list(first.genes_proximal_to_break2)[0]
        assert dist == 50
        assert len(ann_list[1].encompassed_genes) == 2

    def test_interchromosomal(self):
        raise unittest.SkipTest('TODO')

    def test_intrachromosomal_within_gene_inversion(self):
        raise unittest.SkipTest('TODO')
        g = Gene(REF_CHR, 1000, 3000, strand=STRAND.POS)
        t = PreTranscript(gene=g, exons=[(1001, 1100), (1501, 1600), (2001, 2100), (2501, 2600)])
        g.transcripts.append(t)
        ref = {REF_CHR: [g]}
        b1 = Breakpoint(REF_CHR, 1250, strand=STRAND.POS)
        b2 = Breakpoint(REF_CHR, 2250, strand=STRAND.NEG)
        bpp = BreakpointPair(b1, b2)
        ann_list = sorted(_gather_annotations(ref, bpp), key=lambda x: (x.break1, x.break2))
        assert len(ann_list) == 1
        assert ann_list[0].transcript2 == ann_list[0].transcript1

    def test_breakpoint_single_gene(self):
        g = Gene(REF_CHR, 1000, 3000, strand=STRAND.POS)
        t = PreTranscript(gene=g, exons=[(1001, 1100), (1501, 1600), (2001, 2100), (2501, 2600)])
        g.transcripts.append(t)
        ref = {REF_CHR: [g]}
        b1 = Breakpoint(REF_CHR, 789, 1050, strand=STRAND.POS)
        b2 = Breakpoint(REF_CHR, 800, strand=STRAND.POS)
        bpp = BreakpointPair(b1, b2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME)
        ann_list = sorted(_gather_annotations(ref, bpp), key=lambda x: (x.break1, x.break2))
        assert len(ann_list) == 3
        for ann in ann_list:
            assert ann.break1.start in ann.transcript1.position
            assert ann.break1.end in ann.transcript1.position
            assert ann.break2.start in ann.transcript2.position
            assert ann.break2.end in ann.transcript2.position


class TestAnnotate:
    def test_reference_name_eq(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        assert second == first

    def test_reference_name_set(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        d = {first, second}
        assert len(d) == 1

    def test_reference_name_dict(self):
        first, second = ReferenceName('chr1'), ReferenceName('1')
        d = {first: 1}
        d[second] = 2
        print(d)
        assert len(d) == 1
        d = {first: 1, second: 2}
        assert len(d) == 1

    def test_loading_json_annotations(self):
        annotations = load_annotations(get_data('mock_reference_annotations.json'))
        assert len(annotations.keys()) == 1
        assert len(list(annotations.values())[0]) == 1

    def test_loading_annotations_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_annotations('file.other')

    def test_determine_prime(self):
        tneg = PreTranscript(exons=[(3, 4)], strand=STRAND.NEG)
        tpos = PreTranscript(exons=[(3, 4)], strand=STRAND.POS)
        bleft = Breakpoint(REF_CHR, 1, 2, orient=ORIENT.LEFT)
        bright = Breakpoint(REF_CHR, 1, 2, orient=ORIENT.RIGHT)
        # positive left should be five prime
        assert determine_prime(tpos, bleft) == PRIME.FIVE
        # positive right should be three prime
        assert determine_prime(tpos, bright) == PRIME.THREE
        # negative left should be three prime
        assert determine_prime(tneg, bleft) == PRIME.THREE
        # negative right should be five prime
        assert determine_prime(tneg, bright) == PRIME.FIVE

        with pytest.raises(NotSpecifiedError):
            bleft.orient = ORIENT.NS
            determine_prime(tpos, bleft)

        with pytest.raises(NotSpecifiedError):
            determine_prime(tneg, bleft)

        with pytest.raises(NotSpecifiedError):
            tpos.strand = STRAND.NS
            determine_prime(tpos, bright)

    def test_calculate_orf_nested(self):
        seq = (
            'ATGAACACGCAGGAACACCCACGCTCCCGCTCCTCCTACT'
            'CCAGGGGGAGACCGAACCCGAGAGCGACATCCGGAGCTGG'
            'AAGCCGCTGCAACGGGGGCGCCGGCCTCCCTGCCCCGCAG'
            'TCTCCCGCCCTTTCCCAGACACCGGAGGATCCCAGACCAG'
            'CCCAGTATTCCAGACGGAGAGGGTGCCCGAGTCCTCCGCA'
            'GGGCCCCGCCTCCCTAGGCGCTCGGCGCTGGAGACCTAGG'
            'GGAGCGCGCGAGAGAAACGCTGAGCTCCACCCGGCTGGGG'
            'CGGAGGCTGAGCACTATACCCGGGACTCGAACCCCGCGCC'
            'CTCGTCCCCCAGCTTTAGGCATCCCACGGCGAGCCCTGGC'
            'CCAACGCCCCCGCGCCGCGTCTTCCCCCGCACCCCCGCCC'
            'GGAGCCGAGTCCCGGGTTGCCAGCGCGCCCATCCCGCCCC'
            'CAACAGGTCTCGGGCAGCGGGGAGCGCGCTGCCTGTCCCG'
            'CAGGCGCTGGGAGGCGAGCCGGAGCGAAGCCGAGGGCCGG'
            'CCGCGCGCGCGGGAGGGGTCGCCGAGGGCTGCGAGCCCGG'
            'CCACTTACGCGCGGGCTCTGGGGAACCCCATGGAGAGGAG'
            'CACGTCCAGCGCCGATCCATGCTTGATGGTGCCGGGGCGC'
            'TGTTGGCGGTTCCTCCGGGGGGTGACTTTGCTGTACAGCT'
            'CCTCTCTCGCAGCCATGCCGAGCGGACTGGGGTGGCCGTA'
            'CTGAGCCATGGCTCAGCGGCCAGCCATCGCCAGGGAAGGA'
            'GCCCGAGCCGGATGGTTCGGGAAGGAGGAGGGAGTCGGGC'
            'TGGGCCAGGGGAGGGGGCTCGGGGACCCAGAGCCAGGCAG'
            'GAGGCGGCGGCAGCGGAGGCGGCGCTCCGAGCTTCCCCAG'
            'TGTCGGGACTCTGA'
        )
        orfs = calculate_orf(seq)
        for orf in orfs:
            assert seq[orf.start - 1 : orf.start + 2] == 'ATG'
        orfs = sorted(orfs)
        assert len(orfs) == 2
        assert orfs[0] == Interval(1, 894)
        assert orfs[1] == Interval(590, 724)

        seq = (
            'AAGGAGAGAAAATGGCGTCCACGGATTACAGTACCTATAGCCAAGCTGCAGCGCAGCAGGGCTACAGTGCTTACACCGCCCAGCCCACTCAAGGATATGC'
            'ACAGACCACCCAGGCATATGGGCAACAAAGCTATGGAACCTATGGACAGCCCACTGATGTCAGCTATACCCAGGCTCAGACCACTGCAACCTATGGGCAG'
            'ACCGCCTATGCAACTTCTTATGGACAGCCTCCCACTGGTTATACTACTCCAACTGCCCCCCAGGCATACAGCCAGCCTGTCCAGGGGTATGGCACTGGTG'
            'CTTATGATACCACCACTGCTACAGTCACCACCACCCAGGCCTCCTATGCAGCTCAGTCTGCATATGGCACTCAGCCTGCTTATCCAGCCTATGGGCAGCA'
            'GCCAGCAGCCACTGCACCTACAAGCTATTCCTCTACACAGCCGACTAGTTATGATCAGAGCAGTTACTCTCAGCAGAACACCTATGGGCAACCGAGCAGC'
            'TATGGACAGCAGAGTAGCTATGGTCAACAAAGCAGCTATGGGCAGCAGCCTCCCACTAGTTACCCACCCCAAACTGGATCCTACAGCCAAGCTCCAAGTC'
            'AATATAGCCAACAGAGCAGCAGCTACGGGCAGCAGAGTTCATTCCGACAGGACCACCCCAGTAGCATGGGTGTTTATGGGCAGGAGTCTGGAGGATTTTC'
            'CGGACCAGGAGAGAACCGGAGCATGAGTGGCCCTGATAACCGGGGCAGGGGAAGAGGGGGATTTGATCGTGGAGGCATGAGCAGAGGTGGGCGGGGAGGA'
            'GGACGCGGTGGAATGGGCAGCGCTGGAGAGCGAGGTGGCTTCAATAAGCCTGGTGGACCCATGGATGAAGGACCAGATCTTGATCTAGGCCCACCTGTAG'
            'ATCCAGATGAAGACTCTGACAACAGTGCAATTTATGTACAAGGATTAAATGACAGTGTGACTCTAGATGATCTGGCAGACTTCTTTAAGCAGTGTGGGGT'
            'TGTTAAGATGAACAAGAGAACTGGGCAACCCATGATCCACATCTACCTGGACAAGGAAACAGGAAAGCCCAAAGGCGATGCCACAGTGTCCTATGAAGAC'
            'CCACCCACTGCCAAGGCTGCCGTGGAATGGTTTGATGGGAAAGATTTTCAAGGGAGCAAACTTAAAGTCTCCCTTGCTCGGAAGAAGCCTCCAATGAACA'
            'GTATGCGGGGTGGTCTGCCACCCCGTGAGGGCAGAGGCATGCCACCACCACTCCGTGGAGGTCCAGGAGGCCCAGGAGGTCCTGGGGGACCCATGGGTCG'
            'CATGGGAGGCCGTGGAGGAGATAGAGGAGGCTTCCCTCCAAGAGGACCCCGGGGTTCCCGAGGGAACCCCTCTGGAGGAGGAAACGTCCAGCACCGAGCT'
            'GGAGACTGGCAGTGTCCCAATCCGGGTTGTGGAAACCAGAACTTCGCCTGGAGAACAGAGTGCAACCAGTGTAAGGCCCCAAAGCCTGAAGGCTTCCTCC'
            'CGCCACCCTTTCCGCCCCCGGGTGGTGATCGTGGCAGAGGTGGCCCTGGTGGCATGCGGGGAGGAAGAGGTGGCCTCATGGATCGTGGTGGTCCCGGTGG'
            'AATGTTCAGAGGTGGCCGTGGTGGAGACAGAGGTGGCTTCCGTGGTGGCCGGGGCATGGACCGAGGTGGCTTTGGTGGAGGAAGACGAGGTGGCCCTGGG'
            'GGGCCCCCTGGACCTTTGATGGAACAGATGGGAGGAAGAAGAGGAGGACGTGGAGGACCTGGAAAAATGGATAAAGGCGAGCACCGTCAGGAGCGCAGAG'
            'ATCGGCCCTACTAGATGCAGAGACCCCGCAGAGCTGCATTGACTACCAGATTTATTTTTTAAACCAGAAAATGTTTTAAATTTATAATTCCATATTTATA'
            'ATGTTGGCCACAACATTATGATTATTCCTTGTCTGTACTTTAGTATTTTTCACCATTTGTGAAGAAACATTAAAACAAGTTAAATGGTA'
        )

        orfs = calculate_orf(seq)
        for orf in orfs:
            assert seq[orf.start - 1 : orf.start + 2] == 'ATG'


class TestAnnotateEvents:
    def test_annotate_events(self):
        reference_annotations = load_annotations(get_data('mock_reference_annotations.full.json'))
        b1 = Breakpoint('fakereference9', 658, orient=ORIENT.RIGHT, strand=STRAND.POS)
        b2 = Breakpoint('fakereference9', 10237, orient=ORIENT.RIGHT, strand=STRAND.NEG)
        bpp = BreakpointPair(
            b1,
            b2,
            stranded=True,
            opposing_strands=True,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
        )
        annotations = annotate_events(
            [bpp], reference_genome=REFERENCE_GENOME, annotations=reference_annotations, filters=[]
        )
        assert len(annotations) == 4
        assert annotations[0].transcript1.get_strand() == STRAND.POS
        assert annotations[0].transcript2.get_strand() == STRAND.NEG
        assert annotations[0].transcript1.name == 'ENST00000375851'
        assert annotations[0].transcript2.name is None
        for ann in annotations:
            print(ann.transcript1, ann.transcript2)
        annotations = annotate_events(
            [bpp], reference_genome=REFERENCE_GENOME, annotations=reference_annotations
        )
        assert len(annotations) == 2
