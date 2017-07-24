from mavis.breakpoint import Breakpoint
from mavis.annotate import Gene, usTranscript, Transcript
from mavis.constants import ORIENT, STRAND
from mavis.interval import Interval
from mavis.bam.cache import BamCache
from . import MockRead, mock_read_pair
import unittest
from . import MockBamFileHandle
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

REFERENCE_GENOME = None


class TestComputeExonicDistance(unittest.TestCase):
    def setUp(self):
        self.t1 = usTranscript([(1001, 1100), (1501, 1600), (2001, 2100), (2201, 2300)], strand='+')

    def test_intergenic_exonic(self):
        d = TranscriptomeEvidence.compute_exonic_distance(101, 1550, [self.t1])
        self.assertEqual(Interval(1050, 1450), d)

    def test_intergenic_intergenic(self):
        d = TranscriptomeEvidence.compute_exonic_distance(101, 300, [self.t1])
        self.assertEqual(Interval(200), d)

    def test_no_annotations(self):
        d = TranscriptomeEvidence.compute_exonic_distance(101, 300, [])
        self.assertEqual(Interval(200), d)

    def test_intergenic_intronic(self):
        d = TranscriptomeEvidence.compute_exonic_distance(101, 1400, [self.t1])
        self.assertEqual(Interval(1000, 1300), d)

    def test_empty_intron(self):
        t2 = usTranscript([(1001, 1100), (1501, 1600), (2001, 2200), (2201, 2300)], strand='+')
        d = TranscriptomeEvidence.compute_exonic_distance(1001, 2300, [self.t1, t2])
        self.assertEqual(Interval(400, 1300), d)


class TestComputeFragmentSizes(unittest.TestCase):
    def setUp(self):
        b1 = Breakpoint('1', 1051, 1051, 'L')
        b2 = Breakpoint('1', 1551, 1551, 'R')
        self.read_length = 50
        self.trans_ev = TranscriptomeEvidence(
            {},  # fake the annotations
            b1, b2,
            None, None,  # bam_cache and reference_genome
            opposing_strands=False,
            read_length=self.read_length,
            stdev_fragment_size=100,
            median_fragment_size=100,
            stdev_count_abnormal=1,
        )
        self.genomic_ev = GenomeEvidence(
            b1, b2,
            None, None,  # bam_cache and reference_genome
            opposing_strands=False,
            read_length=self.read_length,
            stdev_fragment_size=100,
            median_fragment_size=100,
            stdev_count_abnormal=1
        )

    def test_genomic_vs_trans_no_annotations(self):
        # should be identical
        read, mate = mock_read_pair(
            MockRead('name', '1', 1051 - self.read_length + 1, 1051, is_reverse=False),
            MockRead('name', '1', 2300, 2300 + self.read_length - 1, is_reverse=True)
        )
        self.assertEqual(
            self.trans_ev.compute_fragment_size(read, mate),
            self.genomic_ev.compute_fragment_size(read, mate)
        )

    def test_reverse_reads(self):
        read, mate = mock_read_pair(
            MockRead('name', '1', 1001, 1100, is_reverse=False),
            MockRead('name', '1', 2201, 2301, is_reverse=True)
        )
        self.assertEqual(Interval(1300), self.genomic_ev.compute_fragment_size(read, mate))
        self.assertEqual(Interval(1300), self.genomic_ev.compute_fragment_size(mate, read))
        self.assertEqual(Interval(1300), self.trans_ev.compute_fragment_size(read, mate))
        self.assertEqual(Interval(1300), self.trans_ev.compute_fragment_size(mate, read))


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
        self.assertEqual(Interval(17277321, 17279701), self.transcriptome_window(b, [ust]))


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

    def test_window_accessors(self):
        ge = GenomeEvidence(
            Breakpoint('1', 1500, orient=ORIENT.LEFT),
            Breakpoint('1', 6001, orient=ORIENT.RIGHT),
            None, None,  # bam_cache and reference_genome
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=500,
            median_fragment_size=100,
            call_error=0,
            stdev_count_abnormal=1
        )
        self.assertEqual(901, ge.outer_window1.start)
        self.assertEqual(1649, ge.outer_window1.end)
        self.assertEqual(6600, ge.outer_window2.end)
        self.assertEqual(5852, ge.outer_window2.start)

        self.assertEqual(1351, ge.inner_window1.start)
        self.assertEqual(1649, ge.inner_window1.end)
        self.assertEqual(6150, ge.inner_window2.end)
        self.assertEqual(5852, ge.inner_window2.start)


class TestGenomeEvidenceAddReads(unittest.TestCase):

    def setUp(self):
        self.ge = GenomeEvidence(
            Breakpoint('1', 1500, orient=ORIENT.LEFT),
            Breakpoint('1', 6001, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle({'1': 0})), None,  # reference_genome
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=500,
            median_fragment_size=100,
            call_error=0,
            stdev_count_abnormal=1
        )
        # outer windows (901, 1649)  (5852, 6600)
        # inner windows (1351, 1649)  (5852, 6150)

    def test_collect_flanking_pair_error_unmapped_read(self):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        read.is_unmapped = True
        with self.assertRaises(ValueError):
            self.ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_mate_unmapped(self):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        mate.is_unmapped = True
        with self.assertRaises(ValueError):
            self.ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_query_names_dont_match(self):
        read, mate = mock_read_pair(
            MockRead('test1', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        with self.assertRaises(ValueError):
            self.ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_error_template_lengths_dont_match(self):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False, template_length=50),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        mate.template_length = 55
        with self.assertRaises(ValueError):
            self.ge.collect_flanking_pair(read, mate)

    def test_collect_flanking_pair_read_low_mq(self):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        read.mapping_quality = 0
        self.assertFalse(self.ge.collect_flanking_pair(read, mate))

    def test_collect_flanking_pair_mate_low_mq(self):
        read, mate = mock_read_pair(
            MockRead('test', 0, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        mate.mapping_quality = 0
        self.assertFalse(self.ge.collect_flanking_pair(read, mate))

    def test_collect_flanking_pair_interchromosomal(self):
        read, mate = mock_read_pair(
            MockRead('test', 1, 900, 1000, is_reverse=False),
            MockRead('test', 0, 6000, 6099, is_reverse=True)
        )
        self.assertFalse(self.ge.collect_flanking_pair(read, mate))
