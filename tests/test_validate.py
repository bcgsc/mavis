from structural_variant.breakpoint import Breakpoint
from structural_variant.annotate import load_reference_genome, Gene, usTranscript, Transcript
from structural_variant.constants import ORIENT, STRAND, CIGAR, PYSAM_READ_FLAGS, SVTYPE, CALL_METHOD
from structural_variant.interval import Interval
from structural_variant.bam.cache import BamCache
from tests import MockRead
import unittest
from tests import REFERENCE_GENOME_FILE, BAM_INPUT, FULL_BAM_INPUT

from structural_variant.validate.evidence import GenomeEvidence, TranscriptomeEvidence
import structural_variant.validate.call as call
from structural_variant.validate.call import EventCall

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(BAM_INPUT)
    global FULL_BAM_CACHE
    FULL_BAM_CACHE = BamCache(FULL_BAM_INPUT)
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


class TestGenomeEvidenceWindow(unittest.TestCase):

    def test_orient_ns(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        w = GenomeEvidence._generate_window(
            b, read_length=100, max_expected_fragment_size=550, call_error=11)
        self.assertEqual(440, w[0])
        self.assertEqual(1560, w[1])
        self.assertEqual(1121, len(w))

    def test_orient_left(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        w = GenomeEvidence._generate_window(
            b, read_length=100, call_error=11, max_expected_fragment_size=550)
        self.assertEqual(440, w[0])
        self.assertEqual(1110, w[1])
        self.assertEqual(671, len(w))

    def test_orient_right(self):
        b = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        w = GenomeEvidence._generate_window(
            b, read_length=100, call_error=11, max_expected_fragment_size=550)
        self.assertEqual(890, w[0])
        self.assertEqual(1560, w[1])
        self.assertEqual(671, len(w))


class TestTranscriptomeEvidenceWindow(unittest.TestCase):

    def setUp(self):
        self.annotations = {}
        gene = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        self.ust = usTranscript(gene=gene, exons=[(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        gene.unspliced_transcripts.append(self.ust)
        for spl in self.ust.generate_splicing_patterns():
            t = Transcript(self.ust, spl)
            self.ust.transcripts.append(t)
        self.annotations[gene.chr] = [gene]
        self.read_length = 100
        self.max_expected_fragment_size = 550
        self.call_error = 11

    def transcriptome_window(self, breakpoint, annotations=None):
        return TranscriptomeEvidence._generate_window(
            breakpoint, [self.ust] if annotations is None else annotations,
            read_length=self.read_length,
            call_error=self.call_error,
            max_expected_fragment_size=self.max_expected_fragment_size
        )

    def genome_window(self, breakpoint):
        return GenomeEvidence._generate_window(
            breakpoint,
            read_length=self.read_length,
            call_error=self.call_error,
            max_expected_fragment_size=self.max_expected_fragment_size
        )

    def test_before_start(self):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

        b = Breakpoint(chr='1', start=500, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_after_end(self):
        b = Breakpoint(chr='1', start=5000, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_exonic_long_exon(self):
        b = Breakpoint(chr='1', start=3200, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_intronic_long_exon(self):
        b = Breakpoint(chr='1', start=2970, orient=ORIENT.RIGHT)
        self.assertEqual(Interval(1440, 3561), self.transcriptome_window(b))

    def test_intronic_long_intron(self):
        b = Breakpoint(chr='1', start=2000, orient=ORIENT.RIGHT)
        self.assertEqual(Interval(1440, 3561), self.transcriptome_window(b))

    def test_intronic_short_exon_right(self):
        b = Breakpoint(chr='1', start=1690, orient=ORIENT.RIGHT)
        self.assertEqual(Interval(1090, 3511), self.transcriptome_window(b))

    def test_intronic_short_exon_left(self):
        b = Breakpoint(chr='1', start=2200, orient=ORIENT.LEFT)
        self.assertEqual(Interval(690, 3111), self.transcriptome_window(b))

    def test_multiple_transcripts(self):
        #  [(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = self.annotations['1'][0]
        t2 = usTranscript(gene=gene, exons=[(1001, 1100), (1200, 1300), (2100, 2200)])
        gene.transcripts.append(t2)
        # 989 - 2561
        # 989 - 3411
        self.assertEqual(Interval(990, 3411), self.transcriptome_window(b, [self.ust, t2]))

    def test_many_small_exons(self):
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
        self.assertEqual(Interval(17277321, 17279702), self.transcriptome_window(b, [ust]))


#@unittest.skip('skip because slow')
class TestFullEvidenceGathering(unittest.TestCase):
    # need to make the assertions more specific by checking the actual names of the reads found in each bin
    # rather than just the counts.
    def genome_evidence(self, break1, break2, opposing_strands):
        return GenomeEvidence(
            break1, break2, FULL_BAM_CACHE, REFERENCE_GENOME,
            opposing_strands=opposing_strands,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            stdev_count_abnormal=3,
            min_flanking_pairs_resolution=3,
            max_sc_preceeding_anchor=3
        )

    def count_original_reads(self, reads):
        count = 0
        for read in reads:
            if not read.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                count += 1
                print(read.query_name)
            elif not read.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                count += 1
        return count

    def test_load_evidence_translocation(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 520, orient=ORIENT.RIGHT),
            Breakpoint('reference19', 964, orient=ORIENT.LEFT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(14, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(20, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(21, len(ev1.flanking_pairs))

        # second example
        ev1 = self.genome_evidence(
            Breakpoint('reference2', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference4', 2000, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(21, self.count_original_reads(ev1.split_reads[0]))
        # one of the reads that appears to look good in the bam is too low quality % match
        self.assertEqual(40, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(57, len(ev1.flanking_pairs))

    def test_load_evidence_inversion(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True
        )

        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(54, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(20, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(104, len(ev1.flanking_pairs))
        
        # second example 
        ev1 = self.genome_evidence(
            Breakpoint('reference7', 15000, orient=ORIENT.RIGHT),
            Breakpoint('reference7', 19000, orient=ORIENT.RIGHT),
            opposing_strands=True
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(15, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(27, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(52, len(ev1.flanking_pairs))

    def test_load_evidence_duplication(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference7', 5000, orient=ORIENT.RIGHT),
            Breakpoint('reference7', 11000, orient=ORIENT.LEFT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(34, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(12, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(65, len(ev1.flanking_pairs))

    def test_load_evidence_deletion(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference20', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference20', 6000, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(22, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(14, self.count_original_reads(ev1.split_reads[1]))
        
        # second example
        ev1 = self.genome_evidence(
            Breakpoint('referenceX', 2000, orient=ORIENT.LEFT),
            Breakpoint('referenceX', 6000, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(4, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(10, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(27, len(ev1.flanking_pairs))
        
        # third example
        ev1 = self.genome_evidence(
            Breakpoint('referenceX', 10000, orient=ORIENT.LEFT),
            Breakpoint('referenceX', 14000, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(8, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(9, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(26, len(ev1.flanking_pairs))

    def test_load_evidence_low_qual_deletion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference19', 5000, 5620, orient=ORIENT.LEFT),
            Breakpoint('reference19', 8620, 9240, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(0, len(ev1.split_reads[0]))
        self.assertEqual(0, len(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.flanking_pairs))


class TestEvidenceGathering(unittest.TestCase):

    def setUp(self):
        # test loading of evidence for event found on reference3 1114 2187
        self.ev1 = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            BAM_CACHE, REFERENCE_GENOME,
            opposing_strands=True,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            stdev_count_abnormal=3,
            min_flanking_pairs_resolution=3
        )

    def test_add_split_read(self):
        ev1_sr = MockRead(query_name='HISEQX1_11:3:1105:15351:25130:split',
                          reference_id=1, cigar=[(4, 68), (7, 82)], reference_start=1114,
                          reference_end=1154, query_alignment_start=110,
                          query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
                          query_alignment_end=150, flag=113,
                          next_reference_id=1, next_reference_start=2341)
        self.ev1.add_split_read(ev1_sr, True)
        self.assertEqual(ev1_sr, list(self.ev1.split_reads[0])[0])

    def test_add_split_read_failure(self):
        # wrong cigar string
        ev1_sr = MockRead(query_name='HISEQX1_11:4:1203:3062:55280:split',
                          reference_id=1, cigar=[(7, 110), (7, 40)], reference_start=1114,
                          reference_end=1154, query_alignment_start=110,
                          query_sequence='CTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATG',
                          query_alignment_end=150, flag=371,
                          next_reference_id=1, next_reference_start=2550)
        self.assertFalse(self.ev1.add_split_read(ev1_sr, True))

    def test_add_flanking_pair(self):
        self.ev1.add_flanking_pair(
            MockRead(
                reference_id=1, reference_start=2214, reference_end=2364, is_reverse=True,
                next_reference_id=1, next_reference_start=1120, mate_is_reverse=True
            ),
            MockRead(
                reference_id=1, reference_start=1120, reference_end=2364, is_reverse=True,
                next_reference_id=1, next_reference_start=1120, mate_is_reverse=True,
                is_read1=False
            )
        )
        self.assertEqual(1, len(self.ev1.flanking_pairs))

    def test_add_flanking_pair_not_overlapping_evidence_window(self):
        # first read in pair does not overlap the first evidence window
        # therefore this should return False and not add to the flanking_pairs
        read = MockRead(
            reference_id=1, reference_start=1903,
            reference_end=2053,
            is_read1=True, is_reverse=True, mate_is_reverse=True, is_paired=True,
            next_reference_id=1, next_reference_start=2052
        )
        mate = MockRead(
            reference_id=1, reference_start=2052,
            reference_end=2053,
            is_read1=False, is_reverse=True, mate_is_reverse=True, is_paired=True,
            next_reference_id=1, next_reference_start=1903
        )
        self.assertFalse(self.ev1.add_flanking_pair(read, mate))
        self.assertEqual(0, len(self.ev1.flanking_pairs))

#    @unittest.skip("demonstrating skipping")
    def test_load_evidence(self):
        print(self.ev1)
        self.ev1.load_evidence()
        print(self.ev1.spanning_reads)
        self.assertEqual(
            2,
            len([r for r in self.ev1.split_reads[0] if not r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT)]))
        self.assertEqual(7, len(self.ev1.flanking_pairs))
        self.assertEqual(
            2,
            len([r for r in self.ev1.split_reads[1] if not r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT)]))

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
        sr11 = MockRead(query_sequence='TTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACT', flag=113)
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
        self.ev1.split_reads = ([], [sr1, sr3, sr7, sr9, sr12, sr15, sr19, sr24])  # subset needed to make a contig
#        self.ev1.split_reads=([],[sr1,sr3,sr5,sr6,sr7,sr8,sr9,sr10,sr11,sr12,sr13,sr14,sr15,sr16,sr17,sr18,sr19,sr20,sr21,sr22,sr23,sr24]) #full set of reads produces different contig from subset.
        # full contig with more read support should be
        # CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT
        self.ev1.half_mapped = ([], [sr2])
        self.ev1.assemble_split_reads()
        print(self.ev1.contigs)
        self.assertEqual(
            'CAACAATATGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATC', self.ev1.contigs[0].sequence)


class TestEventCall(unittest.TestCase):

    def setUp(self):
        self.ev1 = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            BAM_CACHE, REFERENCE_GENOME,
            opposing_strands=True,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            stdev_count_abnormal=3,
            min_flanking_pairs_resolution=3
        )
        self.ev = EventCall(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            source_evidence=self.ev1,
            event_type=SVTYPE.INV,
            call_method=CALL_METHOD.SPLIT
        )

    def test_count_flanking_support_empty(self):
        c = self.ev.count_flanking_support()
        self.assertEqual(3, len(c))
        self.assertEqual((set(), 0, 0), c)

    def test_count_flanking_support(self):
        # 1114 ++
        # 2187 ++
        self.ev.source_evidence.flanking_pairs.add((
            MockRead(
                query_name="test1",
                reference_id=3, next_reference_id=3,
                template_length=500, 
                reference_start=1150, 
                reference_end=1200, 
                next_reference_start=2200,
                is_read1=True,
                is_reverse=True, mate_is_reverse=True
            ), 
            MockRead(
                query_name="test1", 
                reference_id=3, next_reference_id=3,
                template_length=-500, 
                reference_start=2200, 
                reference_end=2250, 
                next_reference_start=1150,
                is_read1=False,
                is_reverse=True, mate_is_reverse=True
            )
        ))
        self.ev.source_evidence.flanking_pairs.add((
            MockRead(
                query_name="test2", 
                reference_id=3, next_reference_id=3,
                template_length=560, 
                reference_start=1150, 
                reference_end=1200, 
                next_reference_start=2200,
                is_read1=True,
                is_reverse=True, mate_is_reverse=True
            ), 
            MockRead(
                query_name="test2",
                reference_id=3, next_reference_id=3,
                template_length=-560, 
                reference_start=2200, 
                reference_end=2250, 
                next_reference_start=1150,
                is_read1=False,
                is_reverse=True, mate_is_reverse=True
            )
        ))
        reads, median, stdev = self.ev.count_flanking_support()
        self.assertEqual(2, len(reads))
        self.assertEqual(530, median)
        self.assertEqual(30, stdev)

    def test_count_split_read_support_empty(self):
        c = self.ev.count_split_read_support()
        self.assertEqual(0, sum([len(x) for x in c]))
        self.assertEqual(4, len(c))

    def test_count_split_read_support(self):
        self.ev.source_evidence.split_reads[0].add(
            MockRead(
                query_name="test1", cigar=[(CIGAR.S, 110), (CIGAR.EQ, 40)],
                reference_start=1114, reference_end=1150
            ))
        self.ev.source_evidence.split_reads[0].add(
            MockRead(
                query_name="test2", cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)],
                reference_start=1108, reference_end=1115
            ))
        self.ev.source_evidence.split_reads[0].add(
            MockRead(
                query_name="test3", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=1114, reference_end=1154,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)]
            ))
        self.ev.source_evidence.split_reads[1].add(
            MockRead(
                query_name="test4", cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)], reference_start=2187
            ))
        self.ev.source_evidence.split_reads[1].add(
            MockRead(
                query_name="test5", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)], reference_start=2187
            ))
        self.ev.source_evidence.split_reads[1].add(
            MockRead(
                query_name="test1", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=2187, reference_end=2307,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)]
            ))
        f, ftgt, s, stgt = self.ev.count_split_read_support()
        self.assertEqual(2, len(f))
        self.assertEqual(1, len(ftgt))
        self.assertEqual(2, len(s))
        self.assertEqual(1, len(stgt))
        self.assertEqual(1, len(f & s))


class TestCallBySupportingReads(unittest.TestCase):

    def setUp(self):
        self.ev = GenomeEvidence(
            Breakpoint('fake', 50, 150, orient=ORIENT.RIGHT),
            Breakpoint('fake', 450, 550, orient=ORIENT.RIGHT),
            None, None,
            opposing_strands=True,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1
        )

    def test_empty(self):
        with self.assertRaises(UserWarning):
            break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

    def test_call_both_by_split_read(self):
        self.ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[0].add(
            MockRead(query_name='t2', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t2', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )

        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, break1.start)
        self.assertEqual(101, break1.end)
        self.assertEqual(501, break2.start)
        self.assertEqual(501, break2.end)

    def test_call_both_by_split_read_low_resolution(self):
        self.ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )

        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, break1.start)
        self.assertEqual(101, break1.end)
        self.assertEqual(501, break2.start)
        self.assertEqual(501, break2.end)

    def test_mixed_split_then_flanking(self):
        self.ev.split_reads[0].add(
            MockRead(
                query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)]
            )
        )
        self.ev.flanking_pairs.add((
            MockRead(query_name='t2', reference_start=150, reference_end=150, next_reference_start=505),
            MockRead(query_name='t2', reference_start=505, reference_end=520, next_reference_start=150)
        ))
        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, break1.start)
        self.assertEqual(101, break1.end)
        self.assertEqual(451, break2.start)
        self.assertEqual(506, break2.end)

    def test_split_flanking_read(self):
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.flanking_pairs.add((
            MockRead(query_name='t2', reference_start=120, reference_end=140, next_reference_start=520),
            MockRead(query_name='t2', reference_start=520, reference_end=520, next_reference_start=120)
        ))
        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(71, break1.start)
        self.assertEqual(121, break1.end)
        self.assertEqual(501, break2.start)
        self.assertEqual(501, break2.end)

    def test_both_by_flanking_pairs(self):
        self.ev.flanking_pairs.add((
            MockRead(
                query_name='t1', reference_start=150, reference_end=150, next_reference_start=500
            ),
            MockRead(
                query_name='t1', reference_start=500, reference_end=520, next_reference_start=150
            )
        ))
        self.ev.flanking_pairs.add((
            MockRead(
                query_name='t2', reference_start=120, reference_end=140, next_reference_start=520
            ),
            MockRead(
                query_name='t2', reference_start=520, reference_end=520, next_reference_start=120
            )
        ))
        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]
        # 120-149  ..... 500-519
        # max frag = 150 - 80 = 70
        self.assertEqual(82, break1.start)
        self.assertEqual(121, break1.end)
        self.assertEqual(452, break2.start)  # 70 - 21 = 49
        self.assertEqual(501, break2.end)

    def test_call_both_by_split_reads_multiple_calls(self):
        self.ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[0].add(
            MockRead(query_name='t2', reference_start=110, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t2', reference_start=520, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )

        evs = call._call_by_supporting_reads(self.ev, SVTYPE.INV)
        self.assertEqual(4, len(evs))


class TestCallByFlankingReads(unittest.TestCase):
    def setUp(self):
        self.ev_LR = GenomeEvidence(
            Breakpoint('fake', 100, orient=ORIENT.LEFT),
            Breakpoint('fake', 200, orient=ORIENT.RIGHT),
            None, None,
            opposing_strands=False,
            read_length=25,
            stdev_fragment_size=25,
            median_fragment_size=100,
            stdev_count_abnormal=2,
            min_flanking_pairs_resolution=1
        )

    def test_call_both_intrachromosomal_LR(self):
        # --LLL-100------------500-RRR-------
        # max fragment size: 100 + 2 * 25 = 150
        # max distance = 150 - read_length = 100
        # coverage ranges: 40->80    600->700
        self.assertEqual(150, self.ev_LR.max_expected_fragment_size)
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=19, reference_end=60, next_reference_start=599),
            MockRead(reference_start=599, reference_end=650, next_reference_start=19)
        ))
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=39, reference_end=80, next_reference_start=649),
            MockRead(reference_start=649, reference_end=675, next_reference_start=39)
        ))
        # add a pair that will be ignored
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=39, reference_end=50, next_reference_start=91),
            MockRead(reference_start=91, reference_end=110, next_reference_start=39)
        ))
        break1, break2 = call._call_by_flanking_pairs(self.ev_LR, SVTYPE.DEL)
        self.assertEqual(80, break1.start)
        self.assertEqual(80 + 39, break1.end)
        self.assertEqual(600 - 24, break2.start)
        self.assertEqual(600, break2.end)

    def test_call_both_intrachromosomal_LR_coverage_overlaps_range(self):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=21, reference_end=60, next_reference_start=80),
            MockRead(reference_start=80, reference_end=120, next_reference_start=21)
        ))
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=41, reference_end=80, next_reference_start=110),
            MockRead(reference_start=110, reference_end=140, next_reference_start=41)
        ))
        # pair to skip
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=39, reference_end=80, next_reference_start=649),
            MockRead(reference_start=649, reference_end=675, next_reference_start=39)
        ))
        break1, break2 = call._call_by_flanking_pairs(self.ev_LR, SVTYPE.INS)
        self.assertEqual(80, break1.start)
        self.assertEqual(80, break1.end) # 119
        self.assertEqual(81, break2.start)
        self.assertEqual(81, break2.end)

    def test_intrachromosomal_flanking_coverage_overlap_error(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=19, reference_end=60, next_reference_start=599),
            MockRead(reference_start=599, reference_end=650, next_reference_start=19)
        ))
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=620, reference_end=80, next_reference_start=780),
            MockRead(reference_start=780, reference_end=820, next_reference_start=620)
        ))
        with self.assertRaises(AssertionError):
            call._call_by_flanking_pairs(self.ev_LR, SVTYPE.DEL)

    def test_coverage_larger_than_max_expected_variance_error(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=19, reference_end=60, next_reference_start=599),
            MockRead(reference_start=599, reference_end=650, next_reference_start=19)
        ))
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=301, reference_end=350, next_reference_start=780),
            MockRead(reference_start=780, reference_end=820, next_reference_start=301)
        ))
        with self.assertRaises(AssertionError):
            call._call_by_flanking_pairs(self.ev_LR, SVTYPE.DEL)

    def test_call_both_close_to_zero(self):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        ev = GenomeEvidence(
            Breakpoint('fake', 100, orient=ORIENT.RIGHT),
            Breakpoint('fake', 500, orient=ORIENT.RIGHT),
            None, None,
            opposing_strands=True,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=180,
            stdev_count_abnormal=2,
            min_flanking_pairs_resolution=1
        )
        ev.flanking_pairs.add((
            MockRead(reference_start=19, reference_end=60, next_reference_start=149),
            MockRead(reference_start=149, reference_end=150, next_reference_start=19)
        ))
        ev.flanking_pairs.add((
            MockRead(reference_start=39, reference_end=80, next_reference_start=199),
            MockRead(reference_start=199, reference_end=200, next_reference_start=39)
        ))
        break1, break2 = call._call_by_flanking_pairs(ev, SVTYPE.INV)

        self.assertEqual(1, break1.start)
        self.assertEqual(20, break1.end)
        self.assertEqual(81, break2.start)
        self.assertEqual(150, break2.end)

    def test_call_first_with_second_given_incompatible_error(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=100, reference_end=120, next_reference_start=200),
            MockRead(reference_start=200, reference_end=220, next_reference_start=100)
        ))
        with self.assertRaises(AssertionError):
            break1, break2 = call._call_by_flanking_pairs(
                self.ev_LR, SVTYPE.INV,
                second_breakpoint_called=Breakpoint(self.ev_LR.break2.chr, 110, orient=ORIENT.RIGHT)
            )

    def test_call_first_with_second_given_and_overlap(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=100, reference_end=120, next_reference_start=200),
            MockRead(reference_start=200, reference_end=220, next_reference_start=100)
        ))
        b2 = Breakpoint(self.ev_LR.break2.chr, 119, 150, orient=ORIENT.RIGHT)
        break1, break2 = call._call_by_flanking_pairs(
            self.ev_LR, SVTYPE.INV,
            second_breakpoint_called=b2
        )
        self.assertEqual(b2, break2)
        self.assertEqual(120, break1.start)
        self.assertEqual(149, break1.end)

    def test_call_second_with_first_given_incompatible_error(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=100, reference_end=120, next_reference_start=200),
            MockRead(reference_start=200, reference_end=220, next_reference_start=100)
        ))
        with self.assertRaises(AssertionError):
            break1, break2 = call._call_by_flanking_pairs(
                self.ev_LR, SVTYPE.INV,
                first_breakpoint_called=Breakpoint(self.ev_LR.break2.chr, 210, orient=ORIENT.LEFT)
            )

    def test_call_second_with_first_given_and_overlap(self):
        self.ev_LR.flanking_pairs.add((
            MockRead(reference_start=100, reference_end=120, next_reference_start=200),
            MockRead(reference_start=200, reference_end=220, next_reference_start=100)
        ))
        b1 = Breakpoint(self.ev_LR.break2.chr, 185, orient=ORIENT.LEFT)
        break1, break2 = call._call_by_flanking_pairs(
            self.ev_LR, SVTYPE.INV,
            first_breakpoint_called=b1
        )
        self.assertEqual(b1, break1)
        self.assertEqual(186, break2.start)
        self.assertEqual(201, break2.end)


class MockEvidence:

    def __init__(self, ref=None):
        self.HUMAN_REFERENCE_GENOME = ref


if __name__ == "__main__":
    unittest.main()
