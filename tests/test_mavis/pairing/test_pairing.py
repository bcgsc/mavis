import unittest

import pytest

from mavis.annotate.genomic import PreTranscript
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CALL_METHOD, COLUMNS, ORIENT, PROTOCOL, STRAND, SVTYPE
from mavis.pairing import pairing


@pytest.fixture
def genomic_event1():
    return BreakpointPair(
        Breakpoint('1', 1),
        Breakpoint('1', 10),
        opposing_strands=True,
        **{
            COLUMNS.event_type: SVTYPE.DEL,
            COLUMNS.call_method: CALL_METHOD.CONTIG,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.protocol: PROTOCOL.GENOME,
        },
    )


@pytest.fixture
def genomic_event2():
    return BreakpointPair(
        Breakpoint('1', 1),
        Breakpoint('1', 10),
        opposing_strands=True,
        **{
            COLUMNS.event_type: SVTYPE.DEL,
            COLUMNS.call_method: CALL_METHOD.CONTIG,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.protocol: PROTOCOL.GENOME,
        },
    )


@pytest.fixture
def unspliced_transcript1():
    return PreTranscript(exons=[(1, 100), (301, 400), (501, 600)], strand=STRAND.POS, name='t1')


@pytest.fixture
def unspliced_transcript2():
    return PreTranscript(
        exons=[(1001, 1100), (1301, 1400), (1501, 1600)], strand=STRAND.POS, name='t2'
    )


@pytest.fixture
def transcripts(unspliced_transcript1, unspliced_transcript2):
    return {
        unspliced_transcript1.name: unspliced_transcript1,
        unspliced_transcript2.name: unspliced_transcript2,
    }


@pytest.fixture
def distances():
    return {CALL_METHOD.CONTIG: 0, CALL_METHOD.FLANK: 0, CALL_METHOD.SPLIT: 10}


class TestPairing:
    def test_genome_protocol_diff_chrom(self, genomic_event1, genomic_event2, transcripts):
        genomic_event2.break1.chr = '2'
        assert not pairing.equivalent(genomic_event1, genomic_event2, transcripts)

    def test_genome_protocol_diff_orient(self, genomic_event1, genomic_event2, transcripts):
        genomic_event2.break1.orient = ORIENT.LEFT
        genomic_event1.break1.orient = ORIENT.RIGHT
        assert not pairing.equivalent(genomic_event1, genomic_event2, transcripts)

    def test_genome_protocol_diff_strand(self, genomic_event1, genomic_event2, transcripts):
        genomic_event2.break1.strand = STRAND.POS
        genomic_event1.break1.strand = STRAND.NEG
        assert not pairing.equivalent(genomic_event1, genomic_event2, transcripts)

    def test_genome_protocol_diff_event_type(self, genomic_event1, genomic_event2, transcripts):
        genomic_event2.data[COLUMNS.event_type] = SVTYPE.DEL
        genomic_event1.data[COLUMNS.event_type] = SVTYPE.INS
        assert not pairing.equivalent(genomic_event1, genomic_event2, transcripts)

    def test_genome_protocol_ns_orient(self, genomic_event1, genomic_event2, transcripts):
        genomic_event2.break1.orient = ORIENT.LEFT
        genomic_event1.break2.orient = ORIENT.RIGHT
        assert pairing.equivalent(genomic_event1, genomic_event2, transcripts)

    def test_genome_protocol_by_contig(
        self, genomic_event1, genomic_event2, transcripts, distances
    ):
        genomic_event1.call_method = CALL_METHOD.CONTIG
        genomic_event2.call_method = CALL_METHOD.CONTIG
        distances[CALL_METHOD.CONTIG] = 0
        distances[CALL_METHOD.SPLIT] = 10
        assert pairing.equivalent(genomic_event1, genomic_event2, distances=distances)

        genomic_event1.break1.start = 2
        genomic_event1.break1.end = 20
        assert not pairing.equivalent(genomic_event1, genomic_event2, distances=distances)

    def test_genome_protocol_by_split(self, genomic_event1, genomic_event2, transcripts, distances):
        genomic_event1.call_method = CALL_METHOD.SPLIT
        genomic_event2.call_method = CALL_METHOD.SPLIT
        assert pairing.equivalent(genomic_event1, genomic_event2, distances=distances)
        distances[CALL_METHOD.FLANK] = 100
        distances[CALL_METHOD.SPLIT] = 10
        genomic_event1.break1.start = 11
        genomic_event1.break1.end = 20
        assert not pairing.equivalent(genomic_event1, genomic_event2, distances=distances)

    def test_genome_protocol_by_flanking(
        self, genomic_event1, genomic_event2, transcripts, distances
    ):
        genomic_event1.call_method = CALL_METHOD.FLANK
        genomic_event2.call_method = CALL_METHOD.FLANK
        assert pairing.equivalent(genomic_event1, genomic_event2, distances=distances)
        distances[CALL_METHOD.FLANK] = 10
        distances[CALL_METHOD.SPLIT] = 100
        genomic_event1.break1.start = 11
        genomic_event1.break1.end = 20
        assert not pairing.equivalent(genomic_event1, genomic_event2, distances=distances)

    def test_mixed_protocol_fusions_same_sequence(
        self, genomic_event1, genomic_event2, transcripts
    ):
        genome_ev = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 10,
            },
        )
        trans_ev = BreakpointPair(
            Breakpoint('1', 50),
            Breakpoint('1', 60),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.TRANS,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 10,
            },
        )
        assert not pairing.equivalent(genome_ev, trans_ev, transcripts)
        genome_ev.data[COLUMNS.fusion_sequence_fasta_id] = 'a'
        trans_ev.data[COLUMNS.fusion_sequence_fasta_id] = 'a'
        assert pairing.inferred_equivalent(genome_ev, trans_ev, transcripts)

    def test_mixed_protocol_fusions_same_sequence_diff_translation(
        self, genomic_event1, genomic_event2, transcripts
    ):
        genome_ev = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: 'a',
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 10,
            },
        )
        trans_ev = BreakpointPair(
            Breakpoint('1', 50),
            Breakpoint('1', 60),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: 'a',
                COLUMNS.protocol: PROTOCOL.TRANS,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 50,
            },
        )
        assert not pairing.inferred_equivalent(genome_ev, trans_ev, transcripts)

    def test_mixed_protocol_fusions_different_sequence(
        self, genomic_event1, genomic_event2, transcripts
    ):
        genome_ev = BreakpointPair(
            Breakpoint('1', 1),
            Breakpoint('1', 10),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: 'a',
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 10,
            },
        )
        trans_ev = BreakpointPair(
            Breakpoint('1', 50),
            Breakpoint('1', 60),
            opposing_strands=True,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: 'b',
                COLUMNS.protocol: PROTOCOL.TRANS,
                COLUMNS.transcript1: None,
                COLUMNS.transcript2: None,
                COLUMNS.fusion_cdna_coding_start: 1,
                COLUMNS.fusion_cdna_coding_end: 10,
            },
        )
        assert not pairing.inferred_equivalent(genome_ev, trans_ev, transcripts)

    def test_mixed_protocol_one_predicted_one_match(
        self, genomic_event1, genomic_event2, transcripts, unspliced_transcript1
    ):
        genome_ev = BreakpointPair(
            Breakpoint('1', 350, orient=ORIENT.LEFT),
            Breakpoint('1', 400, orient=ORIENT.RIGHT),
            opposing_strands=False,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.transcript1: unspliced_transcript1.name,
                COLUMNS.transcript2: None,
            },
        )
        trans_ev = BreakpointPair(
            Breakpoint('1', 350, orient=ORIENT.LEFT),
            Breakpoint('1', 400, orient=ORIENT.RIGHT),
            opposing_strands=False,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.TRANS,
                COLUMNS.transcript1: unspliced_transcript1.name,
                COLUMNS.transcript2: None,
            },
        )
        assert pairing.equivalent(genome_ev, trans_ev, transcripts)
        assert pairing.equivalent(trans_ev, genome_ev, transcripts)

        genome_ev.data[COLUMNS.transcript2] = unspliced_transcript1.name
        genome_ev.data[COLUMNS.transcript1] = None
        trans_ev.data[COLUMNS.transcript2] = unspliced_transcript1.name
        trans_ev.data[COLUMNS.transcript1] = None
        assert pairing.inferred_equivalent(genome_ev, trans_ev, transcripts)
        assert pairing.inferred_equivalent(trans_ev, genome_ev, transcripts)

    def test_mixed_protocol_one_predicted_one_mismatch(
        self, genomic_event1, genomic_event2, transcripts, unspliced_transcript1
    ):
        genome_ev = BreakpointPair(
            Breakpoint('1', 350, orient=ORIENT.LEFT),
            Breakpoint('1', 400, orient=ORIENT.RIGHT),
            opposing_strands=False,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.GENOME,
                COLUMNS.transcript1: unspliced_transcript1.name,
                COLUMNS.transcript2: None,
            },
        )
        trans_ev = BreakpointPair(
            Breakpoint('1', 350, orient=ORIENT.LEFT),
            Breakpoint('1', 400, orient=ORIENT.RIGHT),
            opposing_strands=False,
            **{
                COLUMNS.event_type: SVTYPE.DEL,
                COLUMNS.call_method: CALL_METHOD.CONTIG,
                COLUMNS.fusion_sequence_fasta_id: None,
                COLUMNS.protocol: PROTOCOL.TRANS,
                COLUMNS.transcript1: unspliced_transcript1.name,
                COLUMNS.transcript2: None,
            },
        )
        assert pairing.equivalent(genome_ev, trans_ev, transcripts)
        assert pairing.equivalent(trans_ev, genome_ev, transcripts)

        genome_ev.data[COLUMNS.transcript2] = unspliced_transcript1.name
        genome_ev.data[COLUMNS.transcript1] = None
        trans_ev.data[COLUMNS.transcript2] = unspliced_transcript1.name
        trans_ev.data[COLUMNS.transcript1] = None
        assert pairing.inferred_equivalent(genome_ev, trans_ev, transcripts)
        assert pairing.inferred_equivalent(trans_ev, genome_ev, transcripts)

    def test_mixed_protocol_both_predicted(self):

        raise unittest.SkipTest('TODO')

    def test_mixed_protocol_neither_predicted_one_match(self):
        raise unittest.SkipTest('TODO')

    def test_mixed_protocol_neither_predicted_no_match(self):
        raise unittest.SkipTest('TODO')

    def test_mixed_protocol_neither_predicted_both_match(self):
        raise unittest.SkipTest('TODO')

    def test_transcriptome_protocol(self):
        raise unittest.SkipTest('TODO')


@pytest.fixture
def positive_transcript():
    return PreTranscript([(101, 200), (301, 400), (501, 600)], strand=STRAND.POS)


@pytest.fixture
def negative_transcript():
    return PreTranscript([(101, 200), (301, 400), (501, 600)], strand=STRAND.NEG)


class TestBreakpointPrediction:
    def test_exonic_five_prime(self, positive_transcript):
        b = Breakpoint('1', 350, orient=ORIENT.LEFT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 2
        assert breaks[0].start == 200
        assert breaks[1] == b

    def test_exonic_five_prime_first_exon(self, positive_transcript):
        b = Breakpoint('1', 150, orient=ORIENT.LEFT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 1
        assert breaks[0] == b

    def test_exonic_three_prime(self, positive_transcript):
        b = Breakpoint('1', 350, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 2
        assert breaks[1].start == 501
        assert breaks[0] == b

    def test_exonic_three_prime_last_exon(self, positive_transcript):
        b = Breakpoint('1', 550, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 1
        assert breaks[0] == b

    def test_intronic_five_prime(self, positive_transcript):
        b = Breakpoint('1', 450, orient=ORIENT.LEFT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 1
        assert breaks[0].start == 400

    def test_intronic_three_prime(self, positive_transcript):
        b = Breakpoint('1', 250, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, positive_transcript)
        assert len(breaks) == 1
        assert breaks[0].start == 301

    def test_outside_transcript(self, positive_transcript):
        b = Breakpoint('1', 100, orient=ORIENT.RIGHT)
        with pytest.raises(AssertionError):
            pairing.predict_transcriptome_breakpoint(b, positive_transcript)

    # for neg transcripts
    def test_exonic_three_prime_neg(self, negative_transcript):
        b = Breakpoint('1', 350, orient=ORIENT.LEFT, strand=STRAND.NEG)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 2
        assert breaks[0].start == 200
        assert breaks[1] == b

    def test_intronic_three_prime_neg(self, negative_transcript):
        b = Breakpoint('1', 450, orient=ORIENT.LEFT, strand=STRAND.NEG)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 1
        assert breaks[0].start == 400

    def test_exonic_five_prime_neg_first_exon(self, negative_transcript):
        b = Breakpoint('1', 150, orient=ORIENT.LEFT)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 1
        assert breaks[0] == b

    def test_exonic_three_prime_neg_first_exon(self, negative_transcript):
        b = Breakpoint('1', 150, orient=ORIENT.LEFT)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 1
        assert breaks[0] == b

    def test_exonic_five_prime_neg(self, negative_transcript):
        b = Breakpoint('1', 350, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 2
        assert breaks[1].start == 501
        assert breaks[0] == b

    def test_exonic_five_prime_neg_last_exon(self, negative_transcript):
        b = Breakpoint('1', 550, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 1
        assert breaks[0] == b

    def test_intronic_five_prime_neg(self, negative_transcript):
        b = Breakpoint('1', 250, orient=ORIENT.RIGHT)
        breaks = pairing.predict_transcriptome_breakpoint(b, negative_transcript)
        assert len(breaks) == 1
        assert breaks[0].start == 301


class TestEquivalent:
    def test_useq_uncertainty(self):
        event1 = BreakpointPair(
            Breakpoint('1', 157540650, orient='L'),
            Breakpoint('1', 157540877, orient='R'),
            event_type='deletion',
            call_method='contig',
            untemplated_seq='GCCTGGCCGCA',
        )
        event2 = BreakpointPair(
            Breakpoint('1', 157540661, orient='L'),
            Breakpoint('1', 157540877, orient='R'),
            event_type='deletion',
            call_method='spanning reads',
        )
        assert pairing.equivalent(event1, event2)

    def test_useq_uncertainty2(self):
        event1 = BreakpointPair(
            Breakpoint('1', 32, orient='L'),
            Breakpoint('1', 61, orient='R'),
            event_type='deletion',
            call_method='contig',
            untemplated_seq='A',
        )
        event2 = BreakpointPair(
            Breakpoint('1', 24, orient='L'),
            Breakpoint('1', 61, orient='R'),
            event_type='deletion',
            call_method='contig',
            untemplated_seq='TTTTTTTTT',
        )
        assert pairing.equivalent(event1, event2)
