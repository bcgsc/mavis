from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.validate import *
from structural_variant.annotate import load_reference_genome, Gene, usTranscript, Transcript
from structural_variant.constants import ORIENT, READ_PAIR_TYPE
from structural_variant.read_tools import BamCache
from tests import MockRead
import unittest
from tests import REFERENCE_GENOME_FILE, BAM_INPUT

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(BAM_INPUT)
    global READS
    READS = {}
    for read in BAM_CACHE.fetch('reference3', 1, 8000):
        if read.qname not in READS:
            READS[read.qname] = [None, None]
        if read.is_supplementary:
            continue
        if read.is_read1:
            READS[read.qname][0] = read
        else:
            READS[read.qname][1] = read
    # add a check to determine if it is the expected bam file


class TestEvidenceWindow(unittest.TestCase):
    def setUp(self):
        self.annotations = {}
        gene = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        t1 = usTranscript(gene=gene, exons=[(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        gene.unspliced_transcripts.append(t1)
        for spl in t1.generate_splicing_patterns():
            t = Transcript(t1, spl)
            t1.transcripts.append(t)
        self.annotations[gene.chr] = [gene]

    def test_generate_window_orient_ns(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(440, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_window_orient_left(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(440, w[0])
        self.assertEqual(1110, w[1])
        self.assertEqual(671, len(w))

    def test_generate_window_orient_right(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        w = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(890, w[0])
        self.assertEqual(1560, w[1])

    def test_generate_transcriptome_window_before_start(self):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

        b = Breakpoint(chr='1', start=500, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_after_end(self):
        b = Breakpoint(chr='1', start=5000, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_exonic_long_exon(self):
        b = Breakpoint(chr='1', start=3200, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_exon(self):
        b = Breakpoint(chr='1', start=2970, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_long_intron(self):
        b = Breakpoint(chr='1', start=2000, orient=ORIENT.RIGHT)
        w1 = Evidence.generate_window(
            b, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        w2 = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(w1, w2)

    def test_generate_transcriptome_window_intronic_short_exon_right(self):
        b = Breakpoint(chr='1', start=1690, orient=ORIENT.RIGHT)
        w = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1580, 3500), w)

    def test_generate_transcriptome_window_intronic_short_exon_left(self):
        b = Breakpoint(chr='1', start=2200, orient=ORIENT.LEFT)
        w = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1440, 2310), w)

    def test_generate_transcriptome_window_multiple_transcripts(self):
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = self.annotations['1'][0]
        t2 = usTranscript(gene=gene, exons=[(1000, 1100), (1200, 1300), (2100, 2200)])
        gene.transcripts.append(t2)

        w = Evidence.generate_transcriptome_window(
            b, self.annotations, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1040, 3160), w)

    def test_generate_transcriptome_window(self):
        g = Gene('fake', 17271277, 17279592)
        ust = usTranscript(
            gene=g,
            exons=[
                (17271277, 17271984),
                (17272649, 17272709),
                (17275586, 17275681),
                (17275769, 17275930),
                (17276692, 17276817),
                (17277168, 17277388),
                (17277845, 17277888),
                (17278293, 17278378),
                (17279229, 17279592)
            ])
        g.transcripts.append(ust)
        ref = {'fake': [g]}
        b = Breakpoint(chr='fake', start=17279591, orient=ORIENT.LEFT)
        w1 = Evidence.generate_transcriptome_window(
            b, ref, read_length=100, median_insert_size=250, call_error=10, stdev_insert_size=50, stdev_count_abnormal=2)
        self.assertEqual(Interval(1, 2), w1)


class TestEvidenceGathering(unittest.TestCase):
    def setUp(self):
    # test loading of evidence for event found on reference3 1114 2187
        self.ev1 = Evidence(
            BreakpointPair(
                Breakpoint('reference3',1114,orient=ORIENT.RIGHT),
                Breakpoint('reference3',2187,orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            BAM_CACHE, REFERENCE_GENOME,
            read_length=125,
            stdev_insert_size=100,
            median_insert_size=380,
            stdev_count_abnormal=3,
            min_flanking_reads_resolution=3
        )

    def test_add_split_read(self):
        ev1_sr = MockRead(query_name='HISEQX1_11:3:1105:15351:25130:split',
                           reference_id=1, cigar=[(4,68),(7,82)], reference_start=1114,
                           reference_end=1154, query_alignment_start=110,
                           query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
                          query_alignment_end=150, flag=113,
                          next_reference_id=1, next_reference_start=2341)
        self.ev1.add_split_read(ev1_sr,True)
        self.assertEqual(ev1_sr, list(self.ev1.split_reads[0])[0])


    def test_add_split_read_errors(self):
        #wrong cigar string
        ev1_sr = MockRead(query_name='HISEQX1_11:4:1203:3062:55280:split',
                          reference_id=1, cigar=[(7,110),(7,40)], reference_start=1114,
                          reference_end=1154, query_alignment_start=110,
                          query_sequence='CTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATG',
                          query_alignment_end=150, flag=371,
                          next_reference_id=1, next_reference_start=2550)
        with self.assertRaises(UserWarning):
            self.ev1.add_split_read(ev1_sr,True)


    def test_add_flanking_read(self):
        ev1_fr = MockRead(query_name='HISEQX_11:3:2206:22140:26976:flanking',
                          reference_id=1, reference_start=2214,
                          reference_end=2364, flag=113,
                          next_reference_id=1, next_reference_start=1120)
        self.ev1.add_flanking_read(ev1_fr)
        self.assertEqual(ev1_fr, list(self.ev1.flanking_reads[1])[0])

    def test_add_flanking_read_errors(self):
        #read pairs overlap by 1 but in windows
        ev1_fr = MockRead(reference_id=1, reference_start=1903,
                          reference_end=2053, flag=113,
                          next_reference_id=1, next_reference_start=2052)
        with self.assertRaises(UserWarning):
            self.ev1.add_flanking_read(ev1_fr)

#    @unittest.skip("demonstrating skipping")
    def test_load_evidence(self):
        self.ev1.load_evidence()
        self.assertEqual(4,len(self.ev1.split_reads[0]))
        self.assertEqual(8, len(self.ev1.flanking_reads[0]))

#    @unittest.skip("demonstrating skipping")
    def test_assemble_split_reads(self):
        sr1 = MockRead(query_name='HISEQX1_11:3:1105:15351:25130:split',
                       query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
                       flag=113)
        sr2 = MockRead(query_sequence='GTCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTT', flag=121)
        sr3 = MockRead(query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG', flag=113)
        sr5 = MockRead(query_sequence='CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAAT', flag=113)
        sr6 = MockRead(query_sequence='GCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAA', flag=113)
        sr7 = MockRead(query_sequence='TGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATA', flag=113)
        sr8 = MockRead(query_sequence='CCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTA', flag=113)
        sr9 = MockRead(query_sequence='TGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGC', flag=113)
        sr10 = MockRead(query_sequence='ACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTG', flag=113)
        sr11  = MockRead(query_sequence='TTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACT', flag=113)
        sr12 = MockRead(query_sequence='GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAG', flag=113)
        sr13 = MockRead(query_sequence='TTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATG', flag=113)
        sr14 = MockRead(query_sequence='TTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCT', flag=113)
        sr15 = MockRead(query_sequence='GTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACA', flag=113)
        sr16 = MockRead(query_sequence='CCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATAT', flag=113)
        sr17 = MockRead(query_sequence='CGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGG', flag=113)
        sr18 = MockRead(query_sequence='GGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAA', flag=113)
        sr19 = MockRead(query_sequence='TGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCA', flag=113)
        sr20 = MockRead(query_sequence='CATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTT', flag=113)
        sr21 = MockRead(query_sequence='TCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTT', flag=113)
        sr22 = MockRead(query_sequence='TTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTG', flag=113)
        sr23 = MockRead(query_sequence='TGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTAT', flag=113)
        sr24 = MockRead(query_sequence='CTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTT', flag=113)
        sr25 = MockRead(query_sequence='AGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT', flag=113)
        self.ev1.split_reads=([],[sr1,sr3,sr7,sr9,sr12,sr15,sr19,sr24]) #subset needed to make a contig
#        self.ev1.split_reads=([],[sr1,sr3,sr5,sr6,sr7,sr8,sr9,sr10,sr11,sr12,sr13,sr14,sr15,sr16,sr17,sr18,sr19,sr20,sr21,sr22,sr23,sr24]) #full set of reads produces different contig from subset.
        self.ev1.half_mapped=([],[sr2])
        self.ev1.assemble_split_reads()
        print(self.ev1.contigs)
        self.assertEqual('CAACAATATGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATC', self.ev1.contigs[0].seq)

class TestEventCall(unittest.TestCase):
    def setUp(self):
        self.ev1 = Evidence(
            BreakpointPair(
                Breakpoint('reference3',1114,orient=ORIENT.RIGHT),
                Breakpoint('reference3',2187,orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            BAM_CACHE, REFERENCE_GENOME,
            read_length=125,
            stdev_insert_size=100,
            median_insert_size=380,
            stdev_count_abnormal=3,
            min_flanking_reads_resolution=3
        )
        self.ev = EventCall(
            Breakpoint('reference3',1114,orient=ORIENT.RIGHT),
            Breakpoint('reference3',2187,orient=ORIENT.RIGHT),
            self.ev1,
            [SVTYPE.INV],
            call_method='split reads')

    def test_count_flanking_support_empty(self):
        c = self.ev.count_flanking_support()
        self.assertEqual((0,0,0), c)

    def test_count_flanking_support(self):
        self.ev.evidence.flanking_reads[0].add(MockRead(query_name="test1", template_length=500))
        self.ev.evidence.flanking_reads[0].add(MockRead(query_name="test2", template_length=560))
        c = self.ev.count_flanking_support()
        self.assertEqual((2,530,30), c)

    def test_count_split_read_support_empty(self):
        c = self.ev.count_split_read_support()
        self.assertEqual((0,0,0,0,0), c)


    def test_count_split_read_support(self):
        self.ev.evidence.split_reads[0].add(MockRead(query_name="test1", cigar=[(CIGAR.S,110),(CIGAR.EQ,40)],
                                                     reference_start = 1114, reference_end = 1150 ))
        self.ev.evidence.split_reads[0].add(MockRead(query_name="test2", cigar=[(CIGAR.EQ,30),(CIGAR.S,120)],
                                                     reference_start = 1108, reference_end = 1115))
        self.ev.evidence.split_reads[0].add(MockRead(query_name="test3", cigar=[(CIGAR.S,30),(CIGAR.EQ,120)],
                                                     reference_start = 1114, reference_end=1154,
                                                     tags=[(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN,1)]))
        self.ev.evidence.split_reads[1].add(MockRead(query_name="test4", cigar=[(CIGAR.EQ,30),(CIGAR.S,120)],
                                                     reference_start = 2187))
        self.ev.evidence.split_reads[1].add(MockRead(query_name="test5", cigar=[(CIGAR.S,30),(CIGAR.EQ,120)],
                                                     reference_start = 2187))
        self.ev.evidence.split_reads[1].add(MockRead(query_name="test1", cigar=[(CIGAR.S,30),(CIGAR.EQ,120)],
                                                     reference_start = 2187, reference_end=2307,
                                                     tags=[(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN,1)]))
        c = self.ev.count_split_read_support()
        self.assertEqual((2,1,2,1,1), c)


class TestEvidence(unittest.TestCase):
    def test__call_by_flanking_reads_intra(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 100, orient=ORIENT.LEFT),
                Breakpoint('fake', 500, orient=ORIENT.RIGHT),
                opposing_strands=False
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_flanking_reads_resolution=1
        )
        ev.flanking_reads[0].add(MockRead(reference_start=20, reference_end=60, next_reference_start=600))
        ev.flanking_reads[0].add(MockRead(reference_start=40, reference_end=80, next_reference_start=650))
        break1, break2 = Evidence._call_by_flanking_reads(ev, SVTYPE.DEL)
        self.assertEqual(80, break1.start)
        self.assertEqual(209, break1.end)
        self.assertEqual(461, break2.start)
        self.assertEqual(600, break2.end)

    def test__call_by_flanking_reads_intra_tighten_on_overlap(self):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 100, orient=ORIENT.LEFT),
                Breakpoint('fake', 200, orient=ORIENT.RIGHT),
                opposing_strands=False
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_flanking_reads_resolution=1
        )
        ev.flanking_reads[0].add(
            MockRead(reference_start=20, reference_end=60, next_reference_start=150, template_length=200))
        ev.flanking_reads[0].add(
            MockRead(reference_start=40, reference_end=80, next_reference_start=200, template_length=200))
        break1, break2 = Evidence._call_by_flanking_reads(ev, SVTYPE.DEL)
        self.assertEqual(80, break1.start)
        self.assertEqual(149, break1.end)
        self.assertEqual(81, break2.start)
        self.assertEqual(150, break2.end)

    def test__call_by_flanking_reads_close_to_zero(self):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 100, orient=ORIENT.RIGHT),
                Breakpoint('fake', 500, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_flanking_reads_resolution=1
        )
        ev.flanking_reads[0].add(MockRead(reference_start=20, reference_end=60, next_reference_start=250))
        ev.flanking_reads[0].add(MockRead(reference_start=40, reference_end=80, next_reference_start=300))
        break1, break2 = Evidence._call_by_flanking_reads(ev, SVTYPE.DEL)

        self.assertEqual(0, break1.start)
        self.assertEqual(20, break1.end)
        self.assertEqual(111, break2.start)
        self.assertEqual(250, break2.end)

    def test__call_by_supporting_reads_split_read(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 50,150, orient=ORIENT.RIGHT),
                Breakpoint('fake', 450,550, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1
        )
        ev.split_reads[0].add(MockRead(query_name='t1',reference_start=100, reference_end=120, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.split_reads[1].add(MockRead(query_name='t1',reference_start=500, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.split_reads[0].add(MockRead(query_name='t2',reference_start=100, reference_end=120, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.split_reads[1].add(MockRead(query_name='t2',reference_start=500, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))

        break1, break2 = Evidence._call_by_supporting_reads(ev, SVTYPE.INV)[0]

        self.assertEqual(100, break1.start)
        self.assertEqual(100, break1.end)
        self.assertEqual(500, break2.start)
        self.assertEqual(500, break2.end)

    def test__call_by_supporting_reads_split_read2(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 50,150, orient=ORIENT.RIGHT),
                Breakpoint('fake', 450,550, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1
        )
        ev.split_reads[0].add(MockRead(query_name='t1',reference_start=100, reference_end=120, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.split_reads[1].add(MockRead(query_name='t1',reference_start=500, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))

        break1, break2 = Evidence._call_by_supporting_reads(ev, SVTYPE.INV)[0]

        self.assertEqual(100, break1.start)
        self.assertEqual(100, break1.end)
        self.assertEqual(500, break2.start)
        self.assertEqual(500, break2.end)

    def test__call_by_supporting_reads_flanking_split_read(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 50,150, orient=ORIENT.RIGHT),
                Breakpoint('fake', 450,550, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_reads_resolution=1
           )
        ev.split_reads[0].add(MockRead(query_name='t1',reference_start=100, reference_end=120, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.flanking_reads[1].add(MockRead(query_name='t2',reference_start=505, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)], next_reference_start = 150))
        break1, break2 = Evidence._call_by_supporting_reads(ev, SVTYPE.INV)[0]
        print(break1.start,break1.end,break2.start,break2.end)
        self.assertEqual(100, break1.start)
        self.assertEqual(100, break1.end)
        self.assertEqual(355, break2.start)
        self.assertEqual(505, break2.end)

    def test__call_by_supporting_reads_split_flanking_read(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 50,150, orient=ORIENT.RIGHT),
                Breakpoint('fake', 450,550, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_reads_resolution=1
        )
        ev.split_reads[1].add(MockRead(query_name='t1',reference_start=500, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)]))
        ev.flanking_reads[0].add(MockRead(query_name='t2',reference_start=120, reference_end=140, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)], next_reference_start = 520))
        break1, break2 = Evidence._call_by_supporting_reads(ev, SVTYPE.INV)[0]
        print(break1.start,break1.end,break2.start,break2.end)
        self.assertEqual(0, break1.start)
        self.assertEqual(120, break1.end)
        self.assertEqual(500, break2.start)
        self.assertEqual(500, break2.end)

    def test__call_by_supporting_reads_flanking_read(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 50,150, orient=ORIENT.RIGHT),
                Breakpoint('fake', 450,550, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_reads_resolution=1
        )
        ev.flanking_reads[1].add(MockRead(query_name='t1',reference_start=500, reference_end=520, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)], next_reference_start= 150))
        ev.flanking_reads[0].add(MockRead(query_name='t2',reference_start=120, reference_end=140, cigar = [(CIGAR.S,20),(CIGAR.EQ,20)], next_reference_start = 520))
        break1, break2 = Evidence._call_by_supporting_reads(ev, SVTYPE.INV)[0]
        print(break1.start,break1.end,break2.start,break2.end)
        self.assertEqual(0, break1.start)
        self.assertEqual(120, break1.end)
        self.assertEqual(350, break2.start)
        self.assertEqual(500, break2.end)
        pass

    def test__call_by_supporting_reads_empty(self):
        pass

    def test__call_by_supporting_reads(self):
        pass

    def test_expected_insert_size_range(self):
        ev = Evidence(
            BreakpointPair(
                Breakpoint('fake', 100, orient=ORIENT.RIGHT),
                Breakpoint('fake', 500, orient=ORIENT.RIGHT),
                opposing_strands=True
            ),
            None, None,
            read_length=40,
            stdev_insert_size=25,
            median_insert_size=100,
            stdev_count_abnormal=2
        )
        i = ev.expected_insert_size_range()
        self.assertEqual(50, i.start)
        self.assertEqual(150, i.end)

    def test_read_pair_type_LR(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.LR, Evidence.read_pair_type(r))

    def test_read_pair_type_LR_reverse(self):
        r = MockRead(
            reference_id=1,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.LR, Evidence.read_pair_type(r))

    def test_read_pair_type_LL(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.LL, Evidence.read_pair_type(r))

    def test_read_pair_type_RR(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.RR, Evidence.read_pair_type(r))

    def test_read_pair_type_RL(self):
        r = MockRead(
            reference_id=0,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=True,
            mate_is_reverse=False
        )
        self.assertEqual(READ_PAIR_TYPE.RL, Evidence.read_pair_type(r))

    def test_read_pair_type_RL_reverse(self):
        r = MockRead(
            reference_id=1,
            next_reference_id=0,
            reference_start=1,
            next_reference_start=2,
            is_reverse=False,
            mate_is_reverse=True
        )
        self.assertEqual(READ_PAIR_TYPE.RL, Evidence.read_pair_type(r))

class MockEvidence:
    def __init__(self, ref=None):
        self.HUMAN_REFERENCE_GENOME = ref


if __name__ == "__main__":
    unittest.main()
