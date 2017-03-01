from mavis.breakpoint import Breakpoint
from mavis.annotate import load_reference_genome, Gene, usTranscript, Transcript
from mavis.constants import ORIENT, STRAND, PYSAM_READ_FLAGS
from mavis.interval import Interval
from mavis.bam.cache import BamCache
from . import MockRead, mock_read_pair
import unittest
import itertools
from . import REFERENCE_GENOME_FILE, BAM_INPUT, FULL_BAM_INPUT
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

REFERENCE_GENOME = None
RUN_FULL = int(os.environ.get('RUN_FULL', 0))


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


class TestTraverseExonicDistance(unittest.TestCase):

    def setUp(self):
        self.ust1 = usTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.POS)

    def test_left_before_transcript(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(900, 500, ORIENT.LEFT, [self.ust1])
        self.assertEqual(401, gpos.start)
        self.assertEqual(401, gpos.end)

    def test_left_after_transcript(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(2200, 100, ORIENT.LEFT, [self.ust1])
        self.assertEqual(2101, gpos.start)
        self.assertEqual(2101, gpos.end)

        gpos = TranscriptomeEvidence.traverse_exonic_distance(1900, 500, ORIENT.LEFT, [self.ust1])
        self.assertEqual(901, gpos.start)
        self.assertEqual(1401, gpos.end)

    def test_left_within_transcript_exonic(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(1750, 200, ORIENT.LEFT, [self.ust1])
        self.assertEqual(1051, gpos.start)
        self.assertEqual(1551, gpos.end)

    def test_left_within_exon(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(1750, 20, ORIENT.LEFT, [self.ust1])
        self.assertEqual(1731, gpos.start)
        self.assertEqual(1731, gpos.end)

    def test_left_within_transcript_intronic(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(1600, 150, ORIENT.LEFT, [self.ust1])
        self.assertEqual(1051, gpos.start)
        self.assertEqual(1451, gpos.end)

    def test_right_before_transcript(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(500, 100, ORIENT.RIGHT, [self.ust1])
        self.assertEqual(Interval(599), gpos)

        gpos = TranscriptomeEvidence.traverse_exonic_distance(901, 500, ORIENT.RIGHT, [self.ust1])
        self.assertEqual(Interval(1400, 1900), gpos)

    def test_right_after_transcript(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(2201, 100, ORIENT.RIGHT, [self.ust1])
        self.assertEqual(Interval(2300), gpos)

    def test_right_within_transcript(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(1351, 100, ORIENT.RIGHT, [self.ust1])
        self.assertEqual(Interval(1450, 1750), gpos)

    def test_right_within_exon(self):
        gpos = TranscriptomeEvidence.traverse_exonic_distance(1351, 10, ORIENT.RIGHT, [self.ust1])
        self.assertEqual(Interval(1360), gpos)


class TestComputeExonicDistance(unittest.TestCase):
    pass


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

@unittest.skipIf(not RUN_FULL, 'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestFullEvidenceGathering(unittest.TestCase):
    # need to make the assertions more specific by checking the actual names of the reads found in each bin
    # rather than just the counts.
    def genome_evidence(self, break1, break2, opposing_strands):
        ge = GenomeEvidence(
            break1, break2, FULL_BAM_CACHE, REFERENCE_GENOME,
            opposing_strands=opposing_strands,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            stdev_count_abnormal=3,
            min_flanking_pairs_resolution=3,
            max_sc_preceeding_anchor=3
        )
        print(ge.break1.chr, ge.outer_windows[0])
        print(ge.break1.chr, ge.inner_windows[0])
        print(ge.break2.chr, ge.outer_windows[1])
        print(ge.break2.chr, ge.inner_windows[1])
        return ge

    def count_original_reads(self, reads):
        count = 0
        for read in sorted(reads, key=lambda x: x.query_name):
            if not read.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                count += 1
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
        self.assertEqual(35, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(11, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(64, len(ev1.flanking_pairs))

    def test_load_evidence_deletion1(self):
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
    
    def test_load_evidence_deletion2(self):
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
    
    def test_load_evidence_deletion3(self):
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
    
    def test_load_evidence_deletion4(self):
        # forth example
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 3609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 3818, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        self.assertEqual(20, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(17, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(40, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion1(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference11', 6000, orient=ORIENT.LEFT),
            Breakpoint('reference11', 6003, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))

        self.assertEqual(5, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(3, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(20, len(ev1.spanning_reads))
        self.assertEqual(6, len(ev1.flanking_pairs))
    
    def test_load_evidence_small_deletion2(self):
        # second example
        ev1 = self.genome_evidence(
            Breakpoint('reference11', 10000, orient=ORIENT.LEFT),
            Breakpoint('reference11', 10030, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read, mate in ev1.flanking_pairs:
            print(read.query_name)

        self.assertEqual(27, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(52, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(19, len(ev1.spanning_reads))
        self.assertEqual(7, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test1(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 2001, orient=ORIENT.LEFT),
            Breakpoint('reference12', 2120, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read, mate in ev1.flanking_pairs:
            print(read.query_name)

        self.assertEqual(18, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(16, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(22, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test2(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 3609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 3818, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        self.assertEqual(20, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(17, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(40, len(set(ev1.flanking_pairs)))

    def test_load_evidence_small_deletion_test3(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 8609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 8927, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(27, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(5, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(53, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test4(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 12609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 13123, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        self.assertEqual(33, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(6, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(77, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test5(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 17109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 17899, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        self.assertEqual(19, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(11, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(48, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test6(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 22109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 24330, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        self.assertEqual(18, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(13, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(53, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test7(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 28109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 31827, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        self.assertEqual(39, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(13, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(49, len(ev1.flanking_pairs))

    def test_load_evidence_small_deletion_test8(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 36109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 42159, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        self.assertEqual(59, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(8, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(59, len(ev1.flanking_pairs))


    @unittest.skip('skip because too complex')
    def test_load_evidence_complex_deletion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 6001, orient=ORIENT.LEFT),
            Breakpoint('reference12', 6016, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(76, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(83, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(1, len(ev1.spanning_reads))
        self.assertEqual(2, len(ev1.flanking_pairs))

    @unittest.skip('skip because high coverage')
    def test_load_evidence_small_insertion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference1', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference1', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(17, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(17, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(48, len(ev1.spanning_reads))
        self.assertEqual(4, len(ev1.flanking_pairs))

    @unittest.skip('skip because too high coverage')
    def test_load_evidence_small_insertion_high_coverage(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference9', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference9', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(37, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(52, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(37, len(ev1.spanning_reads))
        self.assertEqual(9, len(ev1.flanking_pairs))

        ev1 = self.genome_evidence(
            Breakpoint('reference16', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference16', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(27, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(52, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(19, len(ev1.spanning_reads))
        self.assertEqual(9, len(ev1.flanking_pairs))

    def test_load_evidence_small_duplication(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 10000, orient=ORIENT.RIGHT),
            Breakpoint('reference12', 10021, orient=ORIENT.LEFT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(29, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(51, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(0, len(ev1.spanning_reads))
        self.assertEqual(0, len(ev1.flanking_pairs))

        #Example 2
        ev1 = self.genome_evidence(
            Breakpoint('reference17', 1974, orient=ORIENT.RIGHT),
            Breakpoint('reference17', 2020, orient=ORIENT.LEFT),
            opposing_strands=False
        )
        ev1.load_evidence()

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        self.assertEqual(25, self.count_original_reads(ev1.split_reads[0]))
        self.assertEqual(56, self.count_original_reads(ev1.split_reads[1]))
        self.assertEqual(3, len(ev1.spanning_reads))
        self.assertEqual(0, len(ev1.flanking_pairs))

    def test_load_evidence_low_qual_deletion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference19', 4847, 4847, orient=ORIENT.LEFT),
            Breakpoint('reference19', 5219, 5219, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        ev1.load_evidence()
        print(len(ev1.spanning_reads))
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
        pair = mock_read_pair(
            MockRead(reference_id=1, reference_start=1903, reference_end=2053, is_reverse=True),
            MockRead(reference_id=1, reference_start=2052, reference_end=2053, is_reverse=True)
        )
        self.assertFalse(self.ev1.add_flanking_pair(*pair))
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
        self.ev1.assemble_contig()
        print(self.ev1.contigs)
        self.assertEqual(
            'CAACAATATGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATC', self.ev1.contigs[0].seq)



class MockEvidence:

    def __init__(self, ref=None):
        self.HUMAN_REFERENCE_GENOME = ref


if __name__ == "__main__":
    unittest.main()
