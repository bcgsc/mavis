import shutil
import unittest
from unittest import mock

from mavis import align
from mavis.annotate.file_io import load_reference_genome
from mavis.assemble import Contig
from mavis.bam.cache import BamCache
import mavis.bam.cigar as _cigar
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CIGAR, ORIENT, reverse_complement, STRAND, SVTYPE
from mavis.interval import Interval
from mavis.validate.evidence import GenomeEvidence
from mavis.validate.constants import DEFAULTS
from mavis.bam.read import SamRead

from . import MockBamFileHandle, MockObject, MockLongString, MockRead
from ..util import get_data

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))


class TestCallReadEvents(unittest.TestCase):

    def test_hardclipping(self):
        read = SamRead(reference_name='15')
        read.reference_start = 71491944
        read.cigar = _cigar.convert_string_to_cigar('12=1D25=113H')
        read.query_sequence = 'GTGTGTGGTGTGGGGTGTGTGGTGTGTGTGGTGTGTG'
        read.is_reverse = True

        expected_bpp = BreakpointPair(
            Breakpoint('15', 71491956, orient='L', strand='-'),
            Breakpoint('15', 71491958, orient='R', strand='-'),
            untemplated_seq='')
        events = align.call_read_events(read, is_stranded=True)
        self.assertEqual(1, len(events))
        self.assertEqual(expected_bpp.break1, events[0].break1)
        self.assertEqual(expected_bpp.break2, events[0].break2)


class TestAlign(unittest.TestCase):
    def setUp(self):
        self.cache = BamCache(MockBamFileHandle({'Y': 23, 'fake': 0, 'reference3': 3}))

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
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
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1
        )
        ev.contigs = [
            Contig(
                'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAG'
                'TCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTG'
                'TTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT', 0)
        ]
        print(ev.contigs[0].seq)
        seq = align.align_sequences({'seq': ev.contigs[0].seq}, BAM_CACHE, REFERENCE_GENOME, aligner_reference=get_data('mock_reference_genome.2bit'), aligner='blat')
        print(seq)
        align.select_contig_alignments(ev, seq)
        print(ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        self.assertEqual(1, alignment.read1.reference_id)
        self.assertEqual(1, alignment.read2.reference_id)
        self.assertEqual(Interval(125, 244), align.query_coverage_interval(alignment.read1))
        self.assertEqual(Interval(117, 244), align.query_coverage_interval(alignment.read2))
        self.assertEqual(1114, alignment.read1.reference_start)
        self.assertEqual(2187, alignment.read2.reference_start)
        self.assertEqual([(CIGAR.S, 125), (CIGAR.EQ, 120)], alignment.read1.cigar)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.EQ, 128)], alignment.read2.cigar)

    @unittest.skipIf(not shutil.which('bwa'), 'missing the command')
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
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1
        )
        ev.contigs = [
            Contig(
                'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAG'
                'TCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTG'
                'TTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT', 0)
        ]
        print(ev.contigs[0].seq)
        seq = align.align_sequences(
            {'seq': ev.contigs[0].seq}, BAM_CACHE, REFERENCE_GENOME,
            aligner_reference=get_data('mock_reference_genome.fa'),
            aligner='bwa mem',
            aligner_output_file='mem.out',
            aligner_fa_input_file='mem.in.fa'
        )
        align.select_contig_alignments(ev, seq)
        print(ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        self.assertEqual(reverse_complement(alignment.read1.query_sequence), alignment.read2.query_sequence)
        self.assertEqual('reference3', alignment.read1.reference_name)
        self.assertEqual('reference3', alignment.read2.reference_name)
        self.assertEqual(1, alignment.read1.reference_id)
        self.assertEqual(1, alignment.read2.reference_id)
        self.assertEqual(Interval(125, 244), align.query_coverage_interval(alignment.read1))
        self.assertEqual(Interval(117, 244), align.query_coverage_interval(alignment.read2))
        self.assertEqual(1114, alignment.read1.reference_start)
        self.assertEqual(2187, alignment.read2.reference_start)
        self.assertEqual([(CIGAR.S, 125), (CIGAR.EQ, 120)], alignment.read1.cigar)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.EQ, 128)], alignment.read2.cigar)

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_deletion(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100
        )
        ev.contigs = [
            Contig(
                'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT'
                'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT', 0)
        ]
        seq = align.align_sequences({'seq': ev.contigs[0].seq}, BAM_CACHE, REFERENCE_GENOME, aligner_reference=get_data('mock_reference_genome.2bit'), aligner='blat')
        for query, reads in seq.items():
            print('>>>', query)
            for read in reads:
                print(repr(read))
        align.select_contig_alignments(ev, seq)
        alignments = list(ev.contigs[0].alignments)
        print('alignments:')
        for aln in alignments:
            print(aln, repr(aln.read1), repr(aln.read2))
        self.assertEqual(1, len(alignments))
        alignment = alignments[0]
        self.assertTrue(alignment.read2 is None)
        self.assertEqual(0, alignment.read1.reference_id)
        self.assertTrue(not alignment.read1.is_reverse)
        self.assertEqual(Interval(0, 175), align.query_coverage_interval(alignment.read1))
        self.assertEqual(1612, alignment.read1.reference_start)
        self.assertEqual([(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)], alignment.read1.cigar)

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_inversion(self):
        raise unittest.SkipTest('TODO')

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_deletion_revcomp(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=BAM_CACHE,
            reference_genome=REFERENCE_GENOME,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100
        )
        seq = 'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT' \
              'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT'
        ev.contigs = [Contig(reverse_complement(seq), 0)]
        align.select_contig_alignments(ev, align.align_sequences({'seq': ev.contigs[0].seq}, BAM_CACHE, REFERENCE_GENOME, aligner_reference=get_data('mock_reference_genome.2bit'), aligner='blat'))
        print('alignments:', ev.contigs[0].alignments)
        alignment = list(ev.contigs[0].alignments)[0]
        print(alignment)
        self.assertTrue(alignment.read2 is None)
        self.assertEqual(0, alignment.read1.reference_id)
        self.assertTrue(alignment.read1.is_reverse)
        self.assertEqual(seq, alignment.read1.query_sequence)
        self.assertEqual(Interval(0, 175), align.query_coverage_interval(alignment.read1))
        self.assertEqual(1612, alignment.read1.reference_start)
        self.assertEqual([(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)], alignment.read1.cigar)


class TestBreakpointContigRemappedDepth(unittest.TestCase):
    def setUp(self):
        self.contig = Contig(' ' * 60, None)
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=10))
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=20))
        self.contig.add_mapped_sequence(MockObject(reference_start=50, reference_end=60))

    def test_break_left_deletion(self):
        b = Breakpoint('10', 1030, 1030, orient=ORIENT.LEFT)
        read = MockRead(
            cigar=_cigar.convert_string_to_cigar('35M10D5I20M'),
            reference_start=999,
            reference_name='10'
        )
        align.SplitAlignment.breakpoint_contig_remapped_depth(b, self.contig, read)


class TestSplitEvents(unittest.TestCase):
    def test_read_with_exons(self):
        contig = MockRead(
            query_sequence='CTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGGATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTACGGGGTGACCGTTTGGGAGTTGATGACCTTTGGATCCAA',
            cigar=_cigar.convert_string_to_cigar('68M678D50M15D34M6472D185M10240D158M891D74M8I5883D29M'),
            reference_name='7',
            reference_id=6,
            reference_start=55241669
        )
        self.assertEqual(6, len(align.call_read_events(contig)))


class TestCallBreakpointPair(unittest.TestCase):

    def test_single_one_event(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        bpps = align.call_read_events(r)
        self.assertEqual(1, len(bpps))
        bpp = bpps[0]
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(18, bpp.break2.start)
        self.assertEqual(18, bpp.break2.end)
        self.assertEqual('GGG', bpp.untemplated_seq)

    def test_ins_and_del(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.M, 5), (CIGAR.D, 7), (CIGAR.M, 5)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        # only report the major del event for now
        bpps = align.call_read_events(r)
        self.assertEqual(2, len(bpps))
        bpp = bpps[0]
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(11, bpp.break2.start)
        self.assertEqual(11, bpp.break2.end)
        self.assertEqual('GGG', bpp.untemplated_seq)
        bpp = bpps[1]
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(15, bpp.break1.start)
        self.assertEqual(15, bpp.break1.end)
        self.assertEqual(23, bpp.break2.start)
        self.assertEqual(23, bpp.break2.end)

    def test_single_insertion(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 8), (CIGAR.M, 5)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        bpp = align.call_read_events(r)[0]
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(11, bpp.break2.start)
        self.assertEqual(11, bpp.break2.end)
        self.assertEqual('GGGTAGCT', bpp.untemplated_seq)

    def test_single_duplication(self):
        r = MockRead(
            name='seq1',
            reference_name='gene3',
            reference_start=27155,
            cigar=[(CIGAR.M, 65), (CIGAR.I, 6), (CIGAR.D, 95), (CIGAR.M, 21), (CIGAR.S, 17)],
            query_sequence='TAGTTGGATCTCTGTGCTGACTGACTGACAGACAGACTTTAGTGTCTGTGTGCTGACTGACAGACAGACTTTAGTGTCTGTGTGCTGACT'
                           'GACAGACTCTAGTAGTGTC'
        )
        bpp = align.call_read_events(r)[0]
        self.assertEqual(27220, bpp.break1.start)
        self.assertEqual(27316, bpp.break2.start)
        self.assertEqual('AGACTT', bpp.untemplated_seq)

    def test_single_duplication_with_leading_untemp(self):
        r = MockRead(
            query_sequence=(
                'CTCCCACCAGGAGCTCGTCCTCACCACGTCCTGCACCAGCACCTCCAGCTCCCGCAGCAGCGCCTCGCCCCCACGGTGCGCGCTCCGCGCCGGTTCC'
                'ATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATCGGCTCCGTGGGTTCCATGGACT'
                'CTGTGGGCTCGGGCCCGACGCGCACGGAGGACTGGAGGACTGGGGCGTGTGTCTGCGGTGCAGGCGAGGCGGGGCGGGC'),
            query_name='duplication_with_untemp',
            reference_id=16,
            reference_name='reference17',
            reference_start=1882,
            cigar=[(CIGAR.EQ, 126), (CIGAR.I, 54), (CIGAR.EQ, 93)],
            is_reverse=False)
        bpp = align.call_read_events(r)[0]
        self.assertEqual('AGGTTCCATGGGCTCCGTAGGTTCCATGGGCTCCGTAGGTTCCATCGGCTCCGT', bpp.untemplated_seq)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)

    def test_single_duplication_with_no_untemp(self):
        r = MockRead(
            query_sequence=(
                'GGATGATTTACCTTGGGTAATGAAACTCAGATTTTGCTGTTGTTTTTGTTCGATTTTGCTGTTGTTTTTGTTCCAAAGTGTTTTATACTGATAAAGCAACC'
                'CCGGTTTAGCATTGCCATTGGTAA'),
            query_name='duplication_with_untemp',
            reference_id=2,
            reference_name='reference3',
            reference_start=1497,
            cigar=[(CIGAR.EQ, 51), (CIGAR.I, 22), (CIGAR.EQ, 52)],
            is_reverse=False)
        # repeat: GATTTTGCTGTTGTTTTTGTTC
        bpp = align.convert_to_duplication(align.call_read_events(r)[0], REFERENCE_GENOME)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(bpp.break2.start, 1548)
        self.assertEqual(bpp.break1.start, 1527)

    def test_single_duplication_with_trailing_untemp(self):
        r = MockRead(
            query_sequence=(
                'GGATGATTTACCTTGGGTAATGAAACTCA'
                'GATTTTGCTGTTGTTTTTGTTC'
                'GATTTTGCTGTTGTTTTTGTTC' 'GTCAA'
                'CAAAGTGTTTTATACTGATAAAGCAACCCCGGTTTAGCATTGCCATTGGTAA'),
            query_name='duplication_with_untemp',
            reference_id=2,
            reference_name='reference3',
            reference_start=1497,
            cigar=[(CIGAR.EQ, 51), (CIGAR.I, 27), (CIGAR.EQ, 52)],
            is_reverse=False)
        # repeat: GATTTTGCTGTTGTTTTTGTTC
        print(r)
        print(REFERENCE_GENOME['reference3'][1497:1497 + 51])
        print(REFERENCE_GENOME['reference3'][1548 - 21:1548 + 1])
        bpp = align.call_read_events(r)[0]
        print(bpp)
        bpp = align.convert_to_duplication(bpp, REFERENCE_GENOME)
        print(bpp)
        self.assertEqual('GTCAA', bpp.untemplated_seq)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(bpp.break2.start, 1548)
        self.assertEqual(bpp.break1.start, 1527)

    def test_read_pair_indel(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT 1-30     1-?
        # r1  AAATTTCCCgggaattccggatcgatcgat 1-9      1-9
        # r2  aaatttcccgggaattccggaTCGATCGAT 22-30    100-108
        # i   ---------GGGAATTCCGGA--------- 10-21    n/a
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 9), (CIGAR.S, 21)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 21), (CIGAR.M, 9)],
            query_sequence=seq,
            is_reverse=False
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('GGGAATTCCGGA', bpp.untemplated_seq)
        self.assertEqual(9, bpp.break1.start)
        self.assertEqual(100, bpp.break2.start)
        self.assertEqual('AAATTTCCC', bpp.break1.seq)
        self.assertEqual('TCGATCGAT', bpp.break2.seq)

    def test_read_pair_deletion(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccggaTCGATCGAT
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 21), (CIGAR.M, 9)],
            query_sequence=seq,
            is_reverse=False
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(100, bpp.break2.start)

    def test_read_pair_translocation(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccggaTCGATCGAT
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='2',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 21), (CIGAR.M, 9)],
            query_sequence=seq,
            is_reverse=False
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('1', bpp.break1.chr)
        self.assertEqual('2', bpp.break2.chr)
        self.assertEqual('', bpp.untemplated_seq)

    def test_read_pair_deletion_overlapping_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat
        # r2  aaatttcccgggaattccGGATCGATCGAT

        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.S, 18), (CIGAR.M, 12)],
            query_sequence=seq,
            is_reverse=False
        )
        self.assertEqual(21, r1.reference_end)
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(103, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATTCCGGA', bpp.break1.seq)
        self.assertEqual('TCGATCGAT', bpp.break2.seq)

    def test_read_pair_inversion_overlapping_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTCCGGAtcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ------------------GGA---------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 21), (CIGAR.S, 9)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.M, 12), (CIGAR.S, 18)],
            query_sequence=reverse_complement(seq),
            is_reverse=True
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(108, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATTCCGGA', bpp.break1.seq)
        self.assertEqual(reverse_complement('TCGATCGAT'), bpp.break2.seq)

    def test_read_pair_large_inversion_overlapping_query_coverage(self):
        s = 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'

        read1 = MockRead(
            reference_id=3, reference_start=1114, cigar=[(CIGAR.S, 125), (CIGAR.EQ, 120)], query_sequence=s,
            is_reverse=False
        )
        read2 = MockRead(
            reference_id=3, reference_start=2187, cigar=[(CIGAR.S, 117), (CIGAR.EQ, 8), (CIGAR.D, 1), (CIGAR.M, 120)],
            query_sequence=reverse_complement(s), is_reverse=True
        )
        bpp = align.call_paired_read_event(read1, read2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(1115, bpp.break1.start)
        self.assertEqual(2188 + 3, bpp.break2.start)
        print(bpp.break1.seq)
        print(bpp.break2.seq)
        self.assertEqual(
            'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAG'
            'GGTTTTCATTTCTGTATGTTAAT', bpp.break1.seq)
        self.assertEqual(
            'GCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATCCA'
            'AATTCTGTGTTTACAGGGCTTTCATGCTCAG', bpp.break2.seq)

    def test_read_pair_inversion_gap_in_query_coverage(self):
        # seq AAATTTCCCGGGAATTCCGGATCGATCGAT
        # r1  AAATTTCCCGGGAATTccggatcgatcgat +
        # r2c aaatttcccgggaattccGGATCGATCGAT -
        # i   ----------------CC------------
        # r2  ATCTATCGATCCggaattcccgggaaattt 100+12 = 111 - 3 = 108
        seq = 'AAATTTCCCGGGAATTCCGGATCGATCGAT'  # 30
        r1 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 16), (CIGAR.S, 14)],
            query_sequence=seq,
            is_reverse=False
        )

        r2 = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=99,
            cigar=[(CIGAR.M, 12), (CIGAR.S, 18)],
            query_sequence=reverse_complement(seq),
            is_reverse=True
        )
        bpp = align.call_paired_read_event(r1, r2, is_stranded=True)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('CC', bpp.untemplated_seq)
        self.assertEqual(16, bpp.break1.start)
        self.assertEqual(111, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATT', bpp.break1.seq)
        self.assertEqual(reverse_complement('GGATCGATCGAT'), bpp.break2.seq)


class TestConvertToDuplication(unittest.TestCase):

    def test_insertion_to_duplication(self):
        # BPP(Breakpoint(3:60204611L), Breakpoint(3:60204612R), opposing=False, seq='CATACATACATACATACATACATACATACATA')
        # insertion contig [seq2] contig_alignment_score: 0.99, contig_alignment_mq: Interval(255, 255)
        # (3:60132614[seq2]140=71788D69=32I86=, None))
        bpp = BreakpointPair(
            Breakpoint('3', 60204611, orient='L'), Breakpoint('3', 60204612, orient='R'),
            untemplated_seq='CATACATACATACATACATACATACATACATA',
            opposing_strands=False
        )
        reference_genome = {'3': MockObject(
            seq=MockLongString('CAGGGTCTGAGCTCTTAACTCTATACTGCCTACATACATACATACATACATACATATATACATACATATATAAATT', offset=60204555))}
        print(reference_genome['3'].seq[60204588:60204588 + 8], 'CATACATA')
        setattr(bpp, 'read1', MockObject(query_sequence='', query_name=None))
        setattr(bpp, 'read2', None)
        event = align.convert_to_duplication(bpp, reference_genome)
        print(event)
        self.assertEqual(ORIENT.RIGHT, event.break1.orient)
        self.assertEqual(60204588, event.break1.start)
        self.assertEqual(ORIENT.LEFT, event.break2.orient)
        self.assertEqual(60204611, event.break2.start)
        # CATACATACATACATACATACATACATACATA
        # ........................********
        self.assertEqual('CATACATA', event.untemplated_seq)

    def test_single_bp_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('3', 121, orient='L'), Breakpoint('3', 122, orient='R'),
            untemplated_seq='T',
            opposing_strands=False
        )
        reference_genome = {'3': MockObject(
            seq=MockLongString('ATCGAGCTACGGATCTTTTTTCGATCGATCAATA', offset=100))}
        print(reference_genome['3'].seq[120 - 10:121])
        setattr(bpp, 'read1', MockObject(query_sequence='', query_name=None))
        setattr(bpp, 'read2', None)
        event = align.convert_to_duplication(bpp, reference_genome)
        print(event)
        self.assertEqual(ORIENT.RIGHT, event.break1.orient)
        self.assertEqual(121, event.break1.start)
        self.assertEqual(ORIENT.LEFT, event.break2.orient)
        self.assertEqual(121, event.break2.start)
        self.assertEqual('', event.untemplated_seq)


class TestSelectContigAlignments(unittest.TestCase):
    def test_inversion_and_deletion(self):
        s = 'CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT'
        evidence = MockObject(
            interchromosomal=False,
            opposing_strands=True,
            break1=MockObject(orient=ORIENT.RIGHT, chr='3'),
            break2=MockObject(orient=ORIENT.RIGHT, chr='3'),
            contigs=[MockObject(seq=s, alignments=set())],
            standardize_read=lambda x: x,
            contig_aln_max_event_size=DEFAULTS.contig_aln_max_event_size,
            contig_aln_merge_inner_anchor=5,
            contig_aln_merge_outer_anchor=DEFAULTS.contig_aln_merge_outer_anchor,
            contig_aln_min_query_consumption=0.9,
            contig_aln_min_extend_overlap=DEFAULTS.contig_aln_min_extend_overlap,
            contig_aln_min_anchor_size=DEFAULTS.contig_aln_min_anchor_size,
            contig_aln_min_score=DEFAULTS.contig_aln_min_score,
            outer_window1=Interval(1000, 1200),
            outer_window2=Interval(2000, 2200),
            reference_genome=None,
            bam_cache=mock.Mock(stranded=False)
        )
        read1 = SamRead(
            reference_id=3, reference_start=1114, cigar=[(CIGAR.S, 125), (CIGAR.EQ, 120)], query_sequence=s,
            is_reverse=False, reference_name='3', alignment_rank=0
        )
        read2 = SamRead(
            reference_id=3, reference_start=2187, cigar=[(CIGAR.S, 117), (CIGAR.EQ, 8), (CIGAR.D, 1), (CIGAR.EQ, 120)],
            query_sequence=reverse_complement(s), is_reverse=True, reference_name='3', alignment_rank=1
        )
        raw_alignments = {s: [read1, read2]}
        align.select_contig_alignments(evidence, raw_alignments)
        alignments = list(evidence.contigs[0].alignments)
        self.assertEqual(2, len(alignments))


class TestGetAlignerVersion(unittest.TestCase):

    def test_get_blat_36x2(self):
        content = 'blat - Standalone BLAT v. 36x2 fast sequence search command line tool\n'
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            self.assertEqual('36x2', align.get_aligner_version(align.SUPPORTED_ALIGNER.BLAT))

    def test_get_blat_36(self):
        content = "blat - Standalone BLAT v. 36 fast sequence search command line tool"
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            self.assertEqual('36', align.get_aligner_version(align.SUPPORTED_ALIGNER.BLAT))

    def test_get_bwa_0_7_15(self):
        content = "\nProgram: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.15-r1140"
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            self.assertEqual('0.7.15-r1140', align.get_aligner_version(align.SUPPORTED_ALIGNER.BWA_MEM))

    def test_get_bwa_0_7_12(self):
        content = "\nProgram: bwa (alignment via Burrows-Wheeler transformation)\nVersion: 0.7.12-r1039"
        with mock.patch('subprocess.getoutput', mock.Mock(return_value=content)):
            self.assertEqual('0.7.12-r1039', align.get_aligner_version(align.SUPPORTED_ALIGNER.BWA_MEM))
