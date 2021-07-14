import argparse
from functools import partial

import pytest
from mavis.annotate.genomic import Gene, PreTranscript, Transcript
from mavis.bam import cigar as _cigar
from mavis.bam.cache import BamCache
from mavis.bam.read import SamRead
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, STRAND
from mavis.interval import Interval
from mavis.schemas import DEFAULTS
from mavis.validate.base import Evidence
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

from . import MockBamFileHandle, MockObject, MockRead, mock_read_pair

REFERENCE_GENOME = None


@pytest.fixture
def distance_setup():
    n = argparse.Namespace()
    n.transcript = PreTranscript(
        [(1001, 1100), (1501, 1600), (2001, 2100), (2201, 2300)], strand='+'
    )
    for patt in n.transcript.generate_splicing_patterns():
        n.transcript.transcripts.append(Transcript(n.transcript, patt))
    n.trans_evidence = MockObject(
        annotations={},
        read_length=100,
        max_expected_fragment_size=550,
        call_error=11,
        overlapping_transcripts={n.transcript},
    )
    setattr(
        n.trans_evidence,
        '_select_transcripts',
        lambda *pos: n.trans_evidence.overlapping_transcripts,
    )
    setattr(
        n.trans_evidence,
        'distance',
        partial(TranscriptomeEvidence.distance, n.trans_evidence),
    )
    return n


class TestDistance:
    def test_exonic(self, distance_setup):
        assert distance_setup.trans_evidence.distance(1001, 1550) == Interval(149)

    def test_intergenic_exonic(self, distance_setup):
        dist = distance_setup.trans_evidence.distance(101, 1550)
        assert dist == Interval(1049, 1049)

    def test_intergenic_intergenic(self, distance_setup):
        dist = distance_setup.trans_evidence.distance(101, 300)
        assert dist == Interval(199)

    def test_aligned_intronic(self, distance_setup):
        dist = distance_setup.trans_evidence.distance(1102, 1499)
        assert dist == Interval(5)

    def test_indel_at_exon_boundary(self, distance_setup):
        assert distance_setup.trans_evidence.distance(1101, 1501) == Interval(2)

    def test_no_annotations(self, distance_setup):
        dist = distance_setup.trans_evidence.distance(101, 300, [])
        assert dist == Interval(199)

    def test_intergenic_intronic(self, distance_setup):
        dist = distance_setup.trans_evidence.distance(101, 1400)
        assert dist == Interval(1101)

    def test_empty_intron(self, distance_setup):
        t2 = PreTranscript([(1001, 1100), (1501, 1600), (2001, 2200), (2201, 2300)], strand='+')
        for patt in t2.generate_splicing_patterns():
            t2.transcripts.append(Transcript(t2, patt))
        print(t2)
        print(distance_setup.trans_evidence.overlapping_transcripts)
        distance_setup.trans_evidence.overlapping_transcripts.add(t2)
        dist = distance_setup.trans_evidence.distance(1001, 2301)
        assert dist == Interval(400, 400)


class TestTransStandardize:
    def test_shift_overaligned(self):
        # qwertyuiopas---kkkkk------dfghjklzxcvbnm
        # ..........      ................
        gene = Gene('1', 1, 1000, strand='+')
        transcript = PreTranscript(exons=[(1, 12), (20, 28)], gene=gene, strand='+')
        for spl_patt in transcript.generate_splicing_patterns():
            transcript.transcripts.append(Transcript(transcript, spl_patt))
        gene.transcripts.append(transcript)
        read = SamRead(
            reference_name='1',
            reference_start=0,
            cigar=_cigar.convert_string_to_cigar('14=7D12='),
            query_sequence='qwertyuiopasdfghjklzxcvbnm',
        )
        evidence = TranscriptomeEvidence(
            annotations={},
            reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnm')},
            bam_cache=MockObject(get_read_reference_name=lambda r: r.reference_name),
            break1=Breakpoint('1', 1, orient='L', strand='+'),
            break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220,
        )
        evidence.overlapping_transcripts.add(transcript)
        new_read = evidence.standardize_read(read)
        assert new_read.cigar == _cigar.convert_string_to_cigar('12=7N14=')

    def test_shift_overaligned_left(self):
        # qwertyuiopasdf---kkkkkdf------ghjklzxcvbnm
        # ..........      ................
        gene = Gene('1', 1, 1000, strand='+')
        transcript = PreTranscript(exons=[(1, 14), (22, 28)], gene=gene, strand='+')
        for spl_patt in transcript.generate_splicing_patterns():
            transcript.transcripts.append(Transcript(transcript, spl_patt))
        gene.transcripts.append(transcript)
        read = SamRead(
            reference_name='1',
            reference_start=0,
            cigar=_cigar.convert_string_to_cigar('12=7D14='),
            query_sequence='qwertyuiopasdfghjklzxcvbnm',
        )
        evidence = TranscriptomeEvidence(
            annotations={},
            reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnmsbcdefhi')},
            bam_cache=MockObject(get_read_reference_name=lambda r: r.reference_name),
            break1=Breakpoint('1', 1, orient='L', strand='+'),
            break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220,
        )
        evidence.overlapping_transcripts.add(transcript)
        new_read = evidence.standardize_read(read)
        assert new_read.cigar == _cigar.convert_string_to_cigar('14=7N12=')

    def test_shift_no_transcripts(self):
        read = SamRead(
            reference_name='1',
            reference_start=0,
            cigar=_cigar.convert_string_to_cigar('14=7D18='),
            query_sequence='qwertyuiopasdfdfghjklzxcvbnm',
        )
        evidence = TranscriptomeEvidence(
            annotations={},
            reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnm')},
            bam_cache=None,
            break1=Breakpoint('1', 1, orient='L', strand='+'),
            break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220,
        )
        new_cigar = evidence.exon_boundary_shift_cigar(read)
        assert new_cigar == _cigar.convert_string_to_cigar('14=7D18=')


@pytest.fixture
def read_length():
    return 50


@pytest.fixture
def trans_evidence(read_length):
    return TranscriptomeEvidence(
        {},  # fake the annotations
        Breakpoint('1', 1051, 1051, 'L'),
        Breakpoint('1', 1551, 1551, 'R'),
        None,
        None,  # bam_cache and reference_genome
        opposing_strands=False,
        read_length=read_length,
        stdev_fragment_size=100,
        median_fragment_size=100,
        config={'validate.stdev_count_abnormal': 1},
    )


@pytest.fixture
def genomic_evidence(read_length):
    return GenomeEvidence(
        Breakpoint('1', 1051, 1051, 'L'),
        Breakpoint('1', 1551, 1551, 'R'),
        None,
        None,  # bam_cache and reference_genome
        opposing_strands=False,
        read_length=read_length,
        stdev_fragment_size=100,
        median_fragment_size=100,
        config={'validate.stdev_count_abnormal': 1},
    )


class TestComputeFragmentSizes:
    def test_genomic_vs_trans_no_annotations(self, genomic_evidence, read_length, trans_evidence):
        # should be identical
        read, mate = mock_read_pair(
            MockRead('name', '1', 1051 - read_length + 1, 1051, is_reverse=False),
            MockRead('name', '1', 2300, 2300 + read_length - 1, is_reverse=True),
        )
        assert genomic_evidence.compute_fragment_size(
            read, mate
        ) == trans_evidence.compute_fragment_size(read, mate)

    def test_reverse_reads(self, genomic_evidence, trans_evidence):
        read, mate = mock_read_pair(
            MockRead('name', '1', 1001, 1100, is_reverse=False),
            MockRead('name', '1', 2201, 2301, is_reverse=True),
        )
        assert genomic_evidence.compute_fragment_size(read, mate) == Interval(1300)
        assert genomic_evidence.compute_fragment_size(mate, read) == Interval(1300)
        assert trans_evidence.compute_fragment_size(read, mate) == Interval(1300)
        assert trans_evidence.compute_fragment_size(mate, read) == Interval(1300)


@pytest.fixture
def traverse_setup():
    n = argparse.Namespace()
    n.transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.POS)
    for patt in n.transcript.generate_splicing_patterns():
        n.transcript.transcripts.append(Transcript(n.transcript, patt))

    n.trans_evidence = MockObject(
        annotations={},
        read_length=100,
        max_expected_fragment_size=550,
        call_error=11,
        overlapping_transcripts={n.transcript},
    )
    setattr(
        n.trans_evidence,
        '_select_transcripts',
        lambda *pos: n.trans_evidence.overlapping_transcripts,
    )
    setattr(
        n.trans_evidence,
        'traverse',
        partial(TranscriptomeEvidence.traverse, n.trans_evidence),
    )
    return n


class TestTraverse:
    def test_left_before_transcript(self, traverse_setup):
        exp_pos = Evidence.traverse(900, 500 - 1, ORIENT.LEFT)
        assert traverse_setup.trans_evidence.traverse(900, 500 - 1, ORIENT.LEFT) == exp_pos

    def test_left_after_transcript(self, traverse_setup):
        exp_pos = Evidence.traverse(2200, 100, ORIENT.LEFT)
        assert traverse_setup.trans_evidence.traverse(2200, 100, ORIENT.LEFT) == exp_pos

    def test_left_at_end(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1900, 500, ORIENT.LEFT)
        assert gpos == Interval(900)

    def test_left_within_transcript_exonic(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1750, 200 - 1, ORIENT.LEFT)
        assert gpos == Interval(1051)

    def test_left_within_exon(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1750, 20 - 1, ORIENT.LEFT)
        assert gpos.start == 1731
        assert gpos.end == 1731

    def test_left_within_transcript_intronic(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1600, 150 - 1, ORIENT.LEFT)
        assert gpos == Interval(1451)

    def test_right_before_transcript(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(500, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(599)

    def test_right_before_transcript2(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(901, 500 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1900)

    def test_right_after_transcript(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(2201, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(2300)

    def test_right_within_transcript(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1351, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1750)

    def test_right_within_exon(self, traverse_setup):
        gpos = traverse_setup.trans_evidence.traverse(1351, 10 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1360)


@pytest.fixture
def tranverse_trans_rev_setup():
    n = argparse.Namespace()
    n.transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.NEG)
    for patt in n.transcript.generate_splicing_patterns():
        n.transcript.transcripts.append(Transcript(n.transcript, patt))

    n.trans_evidence = MockObject(
        annotations={},
        read_length=100,
        max_expected_fragment_size=550,
        call_error=11,
        overlapping_transcripts={n.transcript},
    )
    setattr(
        n.trans_evidence,
        '_select_transcripts',
        lambda *pos: n.trans_evidence.overlapping_transcripts,
    )
    setattr(
        n.trans_evidence,
        'traverse',
        partial(TranscriptomeEvidence.traverse, n.trans_evidence),
    )
    return n


class TestTraverseTransRev:
    def test_left_before_transcript(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(900, 500 - 1, ORIENT.LEFT)
        assert gpos == Interval(401)
        assert GenomeEvidence.traverse(900, 500 - 1, ORIENT.LEFT) == gpos

    def test_left_after_transcript(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(2200, 100, ORIENT.LEFT)
        assert GenomeEvidence.traverse(2200, 100, ORIENT.LEFT) == gpos
        assert gpos == Interval(2100)

    def test_left_after_transcript2(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1900, 500 - 1, ORIENT.LEFT)
        assert gpos == Interval(901)

    def test_left_within_transcript_exonic(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1750, 200 - 1, ORIENT.LEFT)
        assert gpos == Interval(1051)

    def test_left_within_exon(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1750, 20 - 1, ORIENT.LEFT)
        assert gpos.start == 1731
        assert gpos.end == 1731

    def test_left_within_transcript_intronic(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1600, 150 - 1, ORIENT.LEFT)
        assert gpos == Interval(1451)

    def test_right_before_transcript(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(500, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(599)

    def test_right_before_transcript2(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(901, 500 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1900)

    def test_right_after_transcript(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(2201, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(2300)

    def test_right_within_transcript(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1351, 100 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1750)

    def test_right_within_exon(self, tranverse_trans_rev_setup):
        gpos = tranverse_trans_rev_setup.trans_evidence.traverse(1351, 10 - 1, ORIENT.RIGHT)
        assert gpos == Interval(1360)


@pytest.fixture
def trans_window_setup():
    n = argparse.Namespace()
    gene = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
    n.pre_transcript = PreTranscript(
        gene=gene, exons=[(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)]
    )
    gene.unspliced_transcripts.append(n.pre_transcript)
    for spl in n.pre_transcript.generate_splicing_patterns():
        n.pre_transcript.transcripts.append(Transcript(n.pre_transcript, spl))
    n.annotations = {gene.chr: [gene]}
    n.genome_evidence = MockObject(
        annotations={},
        read_length=100,
        max_expected_fragment_size=550,
        config={**DEFAULTS, 'validate.call_error': 11},
    )
    n.trans_evidence = MockObject(
        annotations={},
        read_length=100,
        max_expected_fragment_size=550,
        overlapping_transcripts={n.pre_transcript},
        config={**DEFAULTS, 'validate.call_error': 11},
    )
    setattr(
        n.trans_evidence,
        '_select_transcripts',
        lambda *pos: n.trans_evidence.overlapping_transcripts,
    )
    setattr(
        n.trans_evidence,
        'traverse',
        partial(TranscriptomeEvidence.traverse, n.trans_evidence),
    )
    return n


def transcriptome_window(ev, breakpoint, transcripts=None):
    if transcripts:
        ev.overlapping_transcripts.update(transcripts)
    return TranscriptomeEvidence.generate_window(ev, breakpoint)


class TestTranscriptomeEvidenceWindow:
    def test_before_start(self, trans_window_setup):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b
        ) == GenomeEvidence.generate_window(trans_window_setup.genome_evidence, b)

        b = Breakpoint(chr='1', start=500, orient=ORIENT.RIGHT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b
        ) == GenomeEvidence.generate_window(trans_window_setup.genome_evidence, b)

    def test_after_end(self, trans_window_setup):
        b = Breakpoint(chr='1', start=6000, orient=ORIENT.RIGHT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b
        ) == GenomeEvidence.generate_window(trans_window_setup.genome_evidence, b)

    def test_exonic_long_exon(self, trans_window_setup):
        b = Breakpoint(chr='1', start=3200, orient=ORIENT.RIGHT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b
        ) == GenomeEvidence.generate_window(trans_window_setup.genome_evidence, b)

    def test_intronic_long_exon(self, trans_window_setup):
        b = Breakpoint(chr='1', start=2970, orient=ORIENT.RIGHT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b
        ) == GenomeEvidence.generate_window(trans_window_setup.genome_evidence, b)

    def test_intronic_long_intron(self, trans_window_setup):
        b = Breakpoint(chr='1', start=1800, orient=ORIENT.RIGHT)
        assert transcriptome_window(trans_window_setup.trans_evidence, b) == Interval(1490, 2360)

    def test_intronic_short_exon_right(self, trans_window_setup):
        b = Breakpoint(chr='1', start=1690, orient=ORIENT.RIGHT)
        assert transcriptome_window(trans_window_setup.trans_evidence, b) == Interval(1580, 3500)

    def test_intronic_short_exon_left(self, trans_window_setup):
        b = Breakpoint(chr='1', start=2200, orient=ORIENT.LEFT)
        assert transcriptome_window(trans_window_setup.trans_evidence, b) == Interval(1440, 2310)

    def test_multiple_transcripts(self, trans_window_setup):
        #  [(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = trans_window_setup.annotations['1'][0]
        t2 = PreTranscript(gene=gene, exons=[(1001, 1100), (1200, 1300), (2100, 2200)])
        for patt in t2.generate_splicing_patterns():
            t2.transcripts.append(Transcript(t2, patt))
        gene.transcripts.append(t2)
        # 989 - 2561
        # 989 - 3411
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b, [trans_window_setup.pre_transcript, t2]
        ) == Interval(1040, 3160)

    def test_many_small_exons(self, trans_window_setup):
        g = Gene('fake', 17271277, 17279592, strand='+')
        pre_transcript = PreTranscript(
            gene=g,
            exons=[
                (17271277, 17271984),
                (17272649, 17272709),
                (17275586, 17275681),
                (17275769, 17275930),
                (17276692, 17276817),
                (17277168, 17277388),  # 220
                (17277845, 17277888),  # 44
                (17278293, 17278378),  # 86
                (17279229, 17279592),  # 364
            ],
        )
        g.transcripts.append(pre_transcript)
        for patt in pre_transcript.generate_splicing_patterns():
            pre_transcript.transcripts.append(Transcript(pre_transcript, patt))
        b = Breakpoint(chr='fake', start=17279591, orient=ORIENT.LEFT)
        assert transcriptome_window(
            trans_window_setup.trans_evidence, b, [pre_transcript]
        ) == Interval(17277321, 17279701)


class TestNetSizeTrans:
    def test_net_zero(self):
        transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.POS)
        for patt in transcript.generate_splicing_patterns():
            transcript.transcripts.append(Transcript(transcript, patt))
        trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={transcript},
        )
        setattr(
            trans_evidence,
            '_select_transcripts',
            lambda *pos: trans_evidence.overlapping_transcripts,
        )
        setattr(
            trans_evidence,
            'distance',
            partial(TranscriptomeEvidence.distance, trans_evidence),
        )

        bpp = BreakpointPair(
            Breakpoint('1', 1099, orient=ORIENT.LEFT),
            Breakpoint('1', 1302, orient=ORIENT.RIGHT),
            untemplated_seq='TT',
        )
        dist = partial(TranscriptomeEvidence.distance, trans_evidence)
        assert bpp.net_size() == Interval(-200)
        assert bpp.net_size(dist) == Interval(0)


class TestGenomeEvidenceWindow:
    def test_orient_ns(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        window = GenomeEvidence.generate_window(
            MockObject(
                read_length=100,
                max_expected_fragment_size=550,
                config={**DEFAULTS, 'validate.call_error': 11},
            ),
            bpp,
        )
        assert window.start == 440
        assert window.end == 1560
        assert len(window) == 1121

    def test_orient_left(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        window = GenomeEvidence.generate_window(
            MockObject(
                read_length=100,
                max_expected_fragment_size=550,
                config={**DEFAULTS, 'validate.call_error': 11},
            ),
            bpp,
        )
        assert window.start == 440
        assert window.end == 1110
        assert len(window) == 671

    def test_orient_right(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        window = GenomeEvidence.generate_window(
            MockObject(
                read_length=100,
                max_expected_fragment_size=550,
                config={**DEFAULTS, 'validate.call_error': 11},
            ),
            bpp,
        )
        assert window.start == 890
        assert window.end == 1560
        assert len(window) == 671

    def test_window_accessors(self):
        ge = GenomeEvidence(
            Breakpoint('1', 1500, orient=ORIENT.LEFT),
            Breakpoint('1', 6001, orient=ORIENT.RIGHT),
            None,
            None,  # bam_cache and reference_genome
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=500,
            median_fragment_size=100,
            config={'validate.stdev_count_abnormal': 1, 'validate.call_error': 0},
        )
        assert ge.outer_window1.start == 901
        assert ge.outer_window1.end == 1649
        assert ge.outer_window2.end == 6600
        assert ge.outer_window2.start == 5852

        assert ge.inner_window1.start == 1351
        assert ge.inner_window1.end == 1649
        assert ge.inner_window2.end == 6150
        assert ge.inner_window2.start == 5852


@pytest.fixture
def flanking_ge(read_length):
    return GenomeEvidence(
        Breakpoint('1', 1500, orient=ORIENT.LEFT),
        Breakpoint('1', 6001, orient=ORIENT.RIGHT),
        BamCache(MockBamFileHandle({'1': 0})),
        None,  # reference_genome
        opposing_strands=False,
        read_length=150,
        stdev_fragment_size=500,
        median_fragment_size=100,
        config={'validate.stdev_count_abnormal': 1, 'validate.call_error': 0},
    )
    # outer windows (901, 1649)  (5852, 6600)
    # inner windows (1351, 1649)  (5852, 6150)


class TestGenomeEvidenceAddReads:
    def test_collect_flanking_pair_error_unmapped_read(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        read.is_unmapped = True
        with pytest.raises(ValueError):
            flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_mate_unmapped(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        mate.is_unmapped = True
        with pytest.raises(ValueError):
            flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_query_names_dont_match(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test1', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        with pytest.raises(ValueError):
            flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_template_lengths_dont_match(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False, template_length=50),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        mate.template_length = 55
        with pytest.raises(ValueError):
            flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_read_low_mq(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        read.mapping_quality = 0
        assert not flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_mate_low_mq(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        mate.mapping_quality = 0
        assert not flanking_ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_interchromosomal(self, flanking_ge):
        read, mate = mock_read_pair(
            MockRead('test', 1, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True),
        )
        assert not flanking_ge.collect_flanking_pair(read, mate)
