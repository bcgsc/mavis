import unittest
from mavis.constants import CIGAR, ORIENT, STRAND, reverse_complement
from mavis.breakpoint import BreakpointPair, Breakpoint
from mavis.annotate import load_reference_genome
from . import MockRead
from . import REFERENCE_GENOME_FILE

REFERENCE_GENOME = None
REF_CHR = 'fake'


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME[REF_CHR].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


class TestCallBreakpointPair(unittest.TestCase):

    def test_single_one_event(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(18, bpp.break2.start)
        self.assertEqual(18, bpp.break2.end)
        self.assertEqual('GGG', bpp.untemplated_seq)

    def test_single_delins(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.M, 5), (CIGAR.D, 7), (CIGAR.M, 5)],
            query_sequence='ACTGAATCGTGGGTAGCTGCTAG'
        )
        # only report the major del event for now
        bpp = BreakpointPair.call_breakpoint_pair(r)
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
        bpp = BreakpointPair.call_breakpoint_pair(r)
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
        bpp = BreakpointPair.call_breakpoint_pair(r, REFERENCE_GENOME=REFERENCE_GENOME)
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
        bpp = BreakpointPair.call_breakpoint_pair(r, REFERENCE_GENOME=REFERENCE_GENOME)
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
        bpp = BreakpointPair.call_breakpoint_pair(r, REFERENCE_GENOME=REFERENCE_GENOME)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(bpp.break2.start, 1548)
        self.assertEqual(bpp.break1.start, 1527)

    def test_single_duplication_with_trailing_untemp(self):
        r = MockRead(
            query_sequence=(
                'GGATGATTTACCTTGGGTAATGAAACTCAGATTTTGCTGTTGTTTTTGTTC'
                'GATTTTGCTGTTGTTTTTGTTC' 'GTCAA'
                'CAAAGTGTTTTATACTGATAAAGCAACCCCGGTTTAGCATTGCCATTGGTAA'),
            query_name='duplication_with_untemp',
            reference_id=2,
            reference_name='reference3',
            reference_start=1497,
            cigar=[(CIGAR.EQ, 51), (CIGAR.I, 27), (CIGAR.EQ, 52)],
            is_reverse=False)
        # repeat: GATTTTGCTGTTGTTTTTGTTC
        bpp = BreakpointPair.call_breakpoint_pair(r, REFERENCE_GENOME=REFERENCE_GENOME)
        self.assertEqual('GTCAA', bpp.untemplated_seq)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(bpp.break2.start, 1548)
        self.assertEqual(bpp.break1.start, 1527)

    def test_single_multi_events(self):
        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.I, 3), (CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
                           'CTC'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(20, bpp.break1.start)
        self.assertEqual(20, bpp.break1.end)
        self.assertEqual(28, bpp.break2.start)
        self.assertEqual(28, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual('ACTGAATCGTGGGTAGCTGCTAG', bpp.break1.seq)
        self.assertEqual('GGGCCTGCTC', bpp.break2.seq)

        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.D, 3), (CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(23, bpp.break1.start)
        self.assertEqual(23, bpp.break1.end)
        self.assertEqual(31, bpp.break2.start)
        self.assertEqual(31, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_seq)

        r = MockRead(
            reference_id=0,
            reference_name='1',
            reference_start=0,
            cigar=[(CIGAR.M, 10), (CIGAR.D, 7), (CIGAR.M, 10), (CIGAR.D, 3), (CIGAR.M, 10)],
            query_sequence='ACTGAATCGT'
                           'GGGTAGCTGC'
                           'TAGGGGCCTG'
        )
        bpp = BreakpointPair.call_breakpoint_pair(r)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(10, bpp.break1.start)
        self.assertEqual(10, bpp.break1.end)
        self.assertEqual(18, bpp.break2.start)
        self.assertEqual(18, bpp.break2.end)
        self.assertEqual('', bpp.untemplated_seq)

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
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
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
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual('', bpp.untemplated_seq)
        self.assertEqual(21, bpp.break1.start)
        self.assertEqual(100, bpp.break2.start)

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
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
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
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
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
        # first part of the inversion
        pslx_row = {
            'block_count': 1,
            'tstarts': [1114],
            'block_sizes': [120],
            'qname': 'seq1',
            'tname': 'reference3',
            'qstarts': [125],
            'strand': '+',
            'qseq_full': s,
            'score': 1,
            'qseqs': [
                'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGG'
                'TTTTCATTTCTGTATGTTAAT'],
            'tseqs': [
                'TCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGG'
                'TTTTCATTTCTGTATGTTAAT']
        }
        read1 = MockRead(
            reference_id=3, reference_start=1114, cigar=[(CIGAR.S, 125), (CIGAR.EQ, 120)], query_sequence=s,
            is_reverse=False
        )
        read2 = MockRead(
            reference_id=3, reference_start=2187, cigar=[(CIGAR.S, 117), (CIGAR.EQ, 128)],
            query_sequence=reverse_complement(s), is_reverse=True
        )
        bpp = BreakpointPair.call_breakpoint_pair(read1, read2)
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
        bpp = BreakpointPair.call_breakpoint_pair(r1, r2)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.NEG, bpp.break2.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual('CC', bpp.untemplated_seq)
        self.assertEqual(16, bpp.break1.start)
        self.assertEqual(111, bpp.break2.start)
        self.assertEqual('AAATTTCCCGGGAATT', bpp.break1.seq)
        self.assertEqual(reverse_complement('GGATCGATCGAT'), bpp.break2.seq)


class TestBreakpointSequenceHomology(unittest.TestCase):

    def test_left_pos_right_pos(self):
        b1 = Breakpoint(REF_CHR, 157, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1788, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CAATGC', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

        b1 = Breakpoint(REF_CHR, 589, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 704, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('TTAA', 'ATAGC'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_left_pos_left_neg(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.NEG, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_left_neg_left_pos(self):
        # CCC|AAA ------------ TTT|GGG
        # CCC                      CCC
        #     TTT              TTT
        b1 = Breakpoint(REF_CHR, 1459, strand=STRAND.NEG, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 2914, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CCC', 'TTT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_right_pos_right_neg(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_right_neg_right_pos(self):
        # CCC|AAA ------------ TTT|GGG
        # GGG                      GGG
        #     AAA              AAA
        b1 = Breakpoint(REF_CHR, 1460, strand=STRAND.NEG, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 2915, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('AAA', 'GGG'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_close_del(self):
        # ....TT|TT....
        b1 = Breakpoint(REF_CHR, 1001, strand=STRAND.POS, orient=ORIENT.LEFT)
        b2 = Breakpoint(REF_CHR, 1002, strand=STRAND.POS, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('', ''), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_close_dup(self):
        # ....GATACATTTCTTCTTGAAAA...
        # -------------<=============
        # ===============>-----------
        # -------------CT-CT--------- first break homology
        # ------------T--T----------- second break homology
        b1 = Breakpoint(REF_CHR, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        self.assertEqual(('CT', 'TT'), bpp.breakpoint_sequence_homology(REFERENCE_GENOME))

    def test_non_specific_error(self):
        b1 = Breakpoint(REF_CHR, 740, 745, strand=STRAND.POS, orient=ORIENT.RIGHT)
        b2 = Breakpoint(REF_CHR, 747, strand=STRAND.POS, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2)
        with self.assertRaises(AttributeError):
            bpp.breakpoint_sequence_homology(REFERENCE_GENOME)
