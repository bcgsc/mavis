import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from mavis.breakpoint import Breakpoint
from mavis.annotate import load_reference_genome, Gene, usTranscript, Transcript
from mavis.constants import ORIENT, STRAND, CIGAR, PYSAM_READ_FLAGS, SVTYPE, CALL_METHOD
from mavis.interval import Interval
from mavis.bam.cache import BamCache
from tests import MockRead, mock_read_pair
import unittest
from tests import REFERENCE_GENOME_FILE, BAM_INPUT, FULL_BAM_INPUT, MockBamFileHandle
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence
import mavis.validate.call as call
from mavis.validate.call import EventCall

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
            {}, # fake the annotations
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
