from unittest import mock

import mavis.bam.cigar as _cigar
from mavis.annotate.file_io import load_reference_genome
from mavis.bam.cache import BamCache
from mavis.bam.read import SamRead
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CIGAR, ORIENT, STRAND, reverse_complement
from mavis.interval import Interval
from mavis.validate import align
from mavis.validate.assemble import Contig
from mavis.validate.evidence import GenomeEvidence
from mavis_config import DEFAULTS

from ..util import blat_only, bwa_only, get_data
from .mock import MockLongString, MockObject, MockRead

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))


class TestCallReadEvents:
    def test_hardclipping(self):
        read = SamRead(reference_name='15')
        read.reference_start = 71491944
        read.cigar = _cigar.convert_string_to_cigar('12=1D25=113H')
        read.query_sequence = 'GTGTGTGGTGTGGGGTGTGTGGTGTGTGTGGTGTGTG'
        read.is_reverse = True

        expected_bpp = BreakpointPair(
            Breakpoint('15', 71491956, orient='L', strand='-'),
            Breakpoint('15', 71491958, orient='R', strand='-'),
            untemplated_seq='',
        )
        events = align.call_read_events(read, is_stranded=True)
        assert len(events) == 1
        assert events[0].break1 == expected_bpp.break1
        assert events[0].break2 == expected_bpp.break2


class TestAlign:
    @blat_only
    def test_blat_contigs(self):
        ev = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
            config={
                'validate.stdev_count_abnormal': 2,
                'validate.min_splits_reads_resolution': 1,
                'validate.min_flanking_pairs_resolution': 1,
            },
        )
        ev.contigs = [
            Contig(
                'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAG'
                'TCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTG'
                'TTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT',
                0,
            )
        ]
        print(ev.contigs[0].seq)
        seq = align.align_sequences(
            {'seq': ev.contigs[0].seq},
            BAM_CACHE,
            REFERENCE_GENOME,
            aligner_reference=get_data('mock_reference_genome.2bit'),
            aligner='blat',
        )
        print(seq)
        align.select_contig_alignments(ev, seq)
        print(ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        assert alignment.read1.reference_id == 1
        assert alignment.read2.reference_id == 1
        assert align.query_coverage_interval(alignment.read1) == Interval(125, 244)
        assert align.query_coverage_interval(alignment.read2) == Interval(117, 244)
        assert alignment.read1.reference_start == 1114
        assert alignment.read2.reference_start == 2187
        assert alignment.read1.cigar == [(CIGAR.S, 125), (CIGAR.EQ, 120)]
        assert alignment.read2.cigar == [(CIGAR.S, 117), (CIGAR.EQ, 128)]

    @bwa_only
    def test_bwa_contigs(self):
        ev = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
            config={
                'validate.stdev_count_abnormal': 2,
                'validate.min_splits_reads_resolution': 1,
                'validate.min_flanking_pairs_resolution': 1,
            },
        )
        ev.contigs = [
            Contig(
                'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAG'
                'TCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTG'
                'TTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT',
                0,
            )
        ]
        print(ev.contigs[0].seq)
        seq = align.align_sequences(
            {'seq': ev.contigs[0].seq},
            BAM_CACHE,
            REFERENCE_GENOME,
            aligner_reference=get_data('mock_reference_genome.fa'),
            aligner='bwa mem',
            aligner_output_file='mem.out',
            aligner_fa_input_file='mem.in.fa',
        )
        align.select_contig_alignments(ev, seq)
        print(ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        assert alignment.read2.query_sequence == reverse_complement(alignment.read1.query_sequence)
        assert alignment.read1.reference_name == 'reference3'
        assert alignment.read2.reference_name == 'reference3'
        assert alignment.read1.reference_id == 1
        assert alignment.read2.reference_id == 1
        assert align.query_coverage_interval(alignment.read1) == Interval(125, 244)
        assert align.query_coverage_interval(alignment.read2) == Interval(117, 244)
        assert alignment.read1.reference_start == 1114
        assert alignment.read2.reference_start == 2187
        assert alignment.read1.cigar == [(CIGAR.S, 125), (CIGAR.EQ, 120)]
        assert alignment.read2.cigar == [(CIGAR.S, 117), (CIGAR.EQ, 128)]

    @blat_only
    def test_blat_contigs_deletion(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
        )
        ev.contigs = [
            Contig(
                'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT'
                'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT',
                0,
            )
        ]
        seq = align.align_sequences(
            {'seq': ev.contigs[0].seq},
            BAM_CACHE,
            REFERENCE_GENOME,
            aligner_reference=get_data('mock_reference_genome.2bit'),
            aligner='blat',
        )
        for query, reads in seq.items():
            print('>>>', query)
            for read in reads:
                print(repr(read))
        align.select_contig_alignments(ev, seq)
        alignments = list(ev.contigs[0].alignments)
        print('alignments:')
        for aln in alignments:
            print(aln, repr(aln.read1), repr(aln.read2))
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment.read2 is None
        assert alignment.read1.reference_id == 0
        assert not alignment.read1.is_reverse
        assert align.query_coverage_interval(alignment.read1) == Interval(0, 175)
        assert alignment.read1.reference_start == 1612
        assert alignment.read1.cigar == [(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)]

    @blat_only
    def test_blat_contigs_deletion_revcomp(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
        )
        seq = (
            'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT'
            'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT'
        )
        ev.contigs = [Contig(reverse_complement(seq), 0)]
        align.select_contig_alignments(
            ev,
            align.align_sequences(
                {'seq': ev.contigs[0].seq},
                BAM_CACHE,
                REFERENCE_GENOME,
                aligner_reference=get_data('mock_reference_genome.2bit'),
                aligner='blat',
            ),
        )
        print('alignments:', ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        print(alignment)
        assert alignment.read2 is None
        assert alignment.read1.reference_id == 0
        assert alignment.read1.is_reverse
        assert alignment.read1.query_sequence == seq
        assert align.query_coverage_interval(alignment.read1) == Interval(0, 175)
        assert alignment.read1.reference_start == 1612
        assert alignment.read1.cigar == [(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)]


class TestBreakpointContigRemappedDepth:
    def test_break_left_deletion(self):
        contig = Contig(' ' * 60, None)
        contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=10))
        contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=20))
        contig.add_mapped_sequence(MockObject(reference_start=50, reference_end=60))

        b = Breakpoint('10', 1030, 1030, orient=ORIENT.LEFT)
        read = MockRead(
            cigarstring='35M10D5I20M',
            reference_start=999,
            reference_name='10',
        )
        align.DiscontinuousAlignment.breakpoint_contig_remapped_depth(b, contig, read)


class TestSplitEvents:
    def test_read_with_exons(self):
        contig = MockRead(
            query_sequence='CTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGGATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTACGGGGTGACCGTTTGGGAGTTGATGACCTTTGGATCCAA',
            cigarstring='68M678D50M15D34M6472D185M10240D158M891D74M8I5883D29M',
            reference_name='7',
            reference_start=55241669,
        )
        assert len(align.call_read_events(contig)) == 6


class TestCallBreakpointPair:
    def test_single_one_event(self):
        r = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='10M3I7D10M',
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG',
        )
        bpps = align.call_read_events(r)
        assert len(bpps) == 1
        bpp = bpps[0]
        assert bpp.opposing_strands is False
        assert bpp.break1.start == 10
        assert bpp.break1.end == 10
        assert bpp.break2.start == 18
        assert bpp.break2.end == 18
        assert bpp.untemplated_seq == 'GGG'

    def test_ins_and_del(self):
        r = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='10M3I5M7D5M',
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG',
        )
        # only report the major del event for now
        bpps = align.call_read_events(r)
        assert len(bpps) == 2
        bpp = bpps[0]
        assert bpp.opposing_strands is False
        assert bpp.break1.start == 10
        assert bpp.break1.end == 10
        assert bpp.break2.start == 11
        assert bpp.break2.end == 11
        assert bpp.untemplated_seq == 'GGG'
        bpp = bpps[1]
        assert bpp.opposing_strands is False
        assert bpp.break1.start == 15
        assert bpp.break1.end == 15
        assert bpp.break2.start == 23
        assert bpp.break2.end == 23

    def test_single_insertion(self):
        r = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='10M8I5M',
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG',
        )
        bpp = align.call_read_events(r)[0]
        assert bpp.opposing_strands is False
        assert bpp.break1.start == 10
        assert bpp.break1.end == 10
        assert bpp.break2.start == 11
        assert bpp.break2.end == 11
        assert bpp.untemplated_seq == 'GGGTAGCT'

    def test_single_duplication(self):
        r = MockRead(
            query_name='seq1',
            reference_name='gene3',
            reference_start=27155,
            cigarstring='65M6I95D21M17S',
            query_sequence='TAGTTGGATCTCTGTGCTGACTGACTGACAGACAGACTTTAGTGTCTGTGTGCTGACTGACAGACAGACTTTAGTGTCTGTGTGCTGACT'
            'GACAGACTCTAGTAGTGTC',
        )
        bpp = align.call_read_events(r)[0]
        assert bpp.break1.start == 27220
        assert bpp.break2.start == 27316
        assert bpp.untemplated_seq == 'AGACTT'

    def test_single_duplication_with_leading_untemp(self):
        r = MockRead(
            query_sequence=(
                'CTCCCACCAGGAGCTCGTCCTCACCACGTCCTGCACCAGCACCTCCAGCTCCCGCAGCAGCGCCTCGCCCCCACGGTGCGCGCTCCGCGCCGGTTCC'
                'ATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATCGGCTCCGTGGGTTCCATGGACT'
                'CTGTGGGCTCGGGCCCGACGCGCACGGAGGACTGGAGGACTGGGGCGTGTGTCTGCGGTGCAGGCGAGGCGGGGCGGGC'
            ),
            query_name='duplication_with_untemp',
            reference_name='reference17',
            reference_start=1882,
            cigarstring='126=54I93=',
            is_reverse=False,
        )
        bpp = align.call_read_events(r)[0]
        assert bpp.untemplated_seq == 'AGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATCGGCTCCGT'
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT

    def test_single_duplication_with_no_untemp(self):
        r = MockRead(
            query_sequence=(
                'GGATGATTTACCTTGGGTAATGAAACTCAGATTTTGCTGTTGTTTTTGTTCGATTTTGCTGTTGTTTTTGTTCCAAAGTGTTTTATACTGATAAAGCAACC'
                'CCGGTTTAGCATTGCCATTGGTAA'
            ),
            query_name='duplication_with_untemp',
            reference_name='reference3',
            reference_start=1497,
            cigarstring='51=22I52=',
            is_reverse=False,
        )
        # repeat: GATTTTGCTGTTGTTTTTGTTC
        bpp = align.convert_to_duplication(align.call_read_events(r)[0], REFERENCE_GENOME)
        assert bpp.untemplated_seq == ''
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.LEFT
        assert 1548 == bpp.break2.start
        assert 1527 == bpp.break1.start

    def test_single_duplication_with_trailing_untemp(self):
        r = MockRead(
            query_sequence=(
                'GGATGATTTACCTTGGGTAATGAAACTCA'
                'GATTTTGCTGTTGTTTTTGTTC'
                'GATTTTGCTGTTGTTTTTGTTC'
                'GTCAA'
                'CAAAGTGTTTTATACTGATAAAGCAACCCCGGTTTAGCATTGCCATTGGTAA'
            ),
            query_name='duplication_with_untemp',
            reference_name='reference3',
            reference_start=1497,
            cigarstring='51=27I52=',
            is_reverse=False,
        )
        # repeat: GATTTTGCTGTTGTTTTTGTTC
        print(r)
        print(REFERENCE_GENOME['reference3'][1497 : 1497 + 51])
        print(REFERENCE_GENOME['reference3'][1548 - 21 : 1548 + 1])
        bpp = align.call_read_events(r)[0]
        print(bpp)
        bpp = align.convert_to_duplication(bpp, REFERENCE_GENOME)
        print(bpp)
        assert bpp.untemplated_seq == 'GTCAA'
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.LEFT
        assert 1548 == bpp.break2.start
        assert 1527 == bpp.break1.start

    def test_read_pair_indel(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT 1-30     1-?
        # r1  AAATTTCCCgggaattccggatcgatcgat 1-9      1-9
        # r2  aaatttcccgggaattccggaTCGATCGAT 22-30    100-108
        # i   ---------GGGAATTCCGGA--------- 10-21    n/a
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='9M21S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='21S9M',
            query_sequence=seq,
            is_reverse=False,
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.POS
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.untemplated_seq == 'GGGAATTCCGGA'
        assert bpp.break1.start == 9
        assert bpp.break2.start == 100
        assert bpp.break1.seq == 'AAATTTCCC'
        assert bpp.break2.seq == 'TCGATCGAT'

    def test_read_pair_deletion(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccggaTCGATCGAT
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='21M9S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='21S9M',
            query_sequence=seq,
            is_reverse=False,
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.POS
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.untemplated_seq == ''
        assert bpp.break1.start == 21
        assert bpp.break2.start == 100

    def test_read_pair_translocation(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccggaTCGATCGAT
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='2',
            reference_start=0,
            cigarstring='21M9S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='21S9M',
            query_sequence=seq,
            is_reverse=False,
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.POS
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '2'
        assert bpp.untemplated_seq == ''

    def test_read_pair_deletion_overlapping_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccGGATCGATCGAT

        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='21M9S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='18S12M',
            query_sequence=seq,
            is_reverse=False,
        )
        assert r1.reference_end == 21
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.POS
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.untemplated_seq == ''
        assert bpp.break1.start == 21
        assert bpp.break2.start == 103
        assert bpp.break1.seq == 'AAATTTCCCGGGAATTCCGGA'
        assert bpp.break2.seq == 'TCGATCGAT'

    def test_read_pair_inversion_overlapping_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ------------------GGA---------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='21M9S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='12M18S',
            query_sequence=reverse_complement(seq),
            is_reverse=True,
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.NEG
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.untemplated_seq == ''
        assert bpp.break1.start == 21
        assert bpp.break2.start == 108
        assert bpp.break1.seq == 'AAATTTCCCGGGAATTCCGGA'
        assert bpp.break2.seq == reverse_complement('TCGATCGAT')

    def test_read_pair_large_inversion_overlapping_query_coverage(self):
        s = 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'

        read1 = MockRead(
            reference_name='3',
            reference_start=1114,
            cigarstring='125S120=',
            query_sequence=s,
            is_reverse=False,
        )
        read2 = MockRead(
            reference_name='3',
            reference_start=2187,
            cigarstring='117S8=1D120M',
            query_sequence=reverse_complement(s),
            is_reverse=True,
        )
        bpp = align.call_paired_read_event(read1, read2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.NEG
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.untemplated_seq == ''
        assert bpp.break1.start == 1115
        assert bpp.break2.start == 2188 + 3
        print(bpp.break1.seq)
        print(bpp.break2.seq)
        assert (
            bpp.break1.seq
            == 'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'
        )
        assert (
            bpp.break2.seq
            == 'GCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCAAATTCTGTGTTTACAGGGCTTTCATGCTCAG'
        )

    def test_read_pair_inversion_gap_in_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTccggatcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ----------------CC------------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_name='1',
            reference_start=0,
            cigarstring='16M14S',
            query_sequence=seq,
            is_reverse=False,
        )

        r2 = MockRead(
            reference_name='1',
            reference_start=99,
            cigarstring='12M18S',
            query_sequence=reverse_complement(seq),
            is_reverse=True,
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        assert bpp.break1.strand == STRAND.POS
        assert bpp.break2.strand == STRAND.NEG
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.untemplated_seq == 'CC'
        assert bpp.break1.start == 16
        assert bpp.break2.start == 111
        assert bpp.break1.seq == 'AAATTTCCCGGGAATT'
        assert bpp.break2.seq == reverse_complement('GGATCGATCGAT')


class TestConvertToDuplication:
    def test_insertion_to_duplication(self):
        # BPP(Breakpoint(3:60204611L), Breakpoint(3:60204612R), opposing=False, seq='CATACATACATACATACATACATACATACATA')
        # insertion contig [seq2] contig_alignment_score: 0.99, contig_alignment_mq: Interval(255, 255)
        # (3:60132614[seq2]140=71788D69=32I86=, None))
        bpp = BreakpointPair(
            Breakpoint('3', 60204611, orient='L'),
            Breakpoint('3', 60204612, orient='R'),
            untemplated_seq='CATACATACATACATACATACATACATACATA',
            opposing_strands=False,
        )
        reference_genome = {
            '3': MockObject(
                seq=MockLongString(
                    'CAGGGTCTGAGCTCTTAACTCTATACTGCCTACATACATACATACATACATACATATATACATACATATATAAATT',
                    offset=60204555,
                )
            )
        }
        print(reference_genome['3'].seq[60204588 : 60204588 + 8], 'CATACATA')
        setattr(bpp, 'read1', MockObject(query_sequence='', query_name=None))
        setattr(bpp, 'read2', None)
        event = align.convert_to_duplication(bpp, reference_genome)
        print(event)
        assert event.break1.orient == ORIENT.RIGHT
        assert event.break1.start == 60204588
        assert event.break2.orient == ORIENT.LEFT
        assert event.break2.start == 60204611
        # CATACATACATACATACATACATACATACATA
        # ........................********
        assert event.untemplated_seq == 'CATACATA'

    def test_single_bp_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('3', 121, orient='L'),
            Breakpoint('3', 122, orient='R'),
            untemplated_seq='T',
            opposing_strands=False,
        )
        reference_genome = {
            '3': MockObject(seq=MockLongString('ATCGAGCTACGGATCTTTTTTCGATCGATCAATA', offset=100))
        }
        print(reference_genome['3'].seq[120 - 10 : 121])
        setattr(bpp, 'read1', MockObject(query_sequence='', query_name=None))
        setattr(bpp, 'read2', None)
        event = align.convert_to_duplication(bpp, reference_genome)
        print(event)
        assert event.break1.orient == ORIENT.RIGHT
        assert event.break1.start == 121
        assert event.break2.orient == ORIENT.LEFT
        assert event.break2.start == 121
        assert event.untemplated_seq == ''


class TestSelectContigAlignments:
    def test_inversion_and_deletion(self):
        s = 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'
        evidence = MockObject(
            interchromosomal=False,
            opposing_strands=True,
            break1=MockObject(orient=ORIENT.RIGHT, chr='3'),
            break2=MockObject(orient=ORIENT.RIGHT, chr='3'),
            contigs=[MockObject(seq=s, alignments=set())],
            standardize_read=lambda x: x,
            config={
                **DEFAULTS,
                'validate.contig_aln_merge_inner_anchor': 5,
                'validate.contig_aln_min_query_consumption': 0.9,
            },
            outer_window1=Interval(1000, 1200),
            outer_window2=Interval(2000, 2200),
            LR=False,
            LL=False,
            RR=True,
            RL=False,
            reference_genome=None,
            bam_cache=mock.Mock(stranded=False),
        )
        read1 = MockRead(
            reference_start=1114,
            cigarstring='125S120=',
            query_sequence=s,
            is_reverse=False,
            reference_name='3',
            alignment_rank=0,
        )
        read2 = MockRead(
            reference_start=2187,
            cigarstring='117S7=1X120=',
            query_sequence=reverse_complement(s),
            is_reverse=True,
            reference_name='3',
            alignment_rank=1,
        )
        raw_alignments = {s: [read1, read2]}
        align.select_contig_alignments(evidence, raw_alignments)
        alignments = list(evidence.contigs[0].alignments)
        assert len(alignments) == 2


class TestGetAlignerVersion:
    def test_get_blat_36x2(self):
        content = 'blat - Standalone BLAT v. 36x2 fast sequence search command line tool\n'
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            assert align.get_aligner_version(align.SUPPORTED_ALIGNER.BLAT) == '36x2'

    def test_get_blat_36(self):
        content = "blat - Standalone BLAT v. 36 fast sequence search command line tool"
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            assert align.get_aligner_version(align.SUPPORTED_ALIGNER.BLAT) == '36'

    def test_get_bwa_0_7_15(self):
        content = (
            "\nProgram: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.15-r1140"
        )
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            assert align.get_aligner_version(align.SUPPORTED_ALIGNER.BWA_MEM) == '0.7.15-r1140'

    def test_get_bwa_0_7_12(self):
        content = (
            "\nProgram: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.12-r1039"
        )
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            assert align.get_aligner_version(align.SUPPORTED_ALIGNER.BWA_MEM) == '0.7.12-r1039'
