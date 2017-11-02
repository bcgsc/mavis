import shutil
import unittest

from mavis.align import align_contigs, query_coverage_interval, SplitAlignment
from mavis.annotate.file_io import load_reference_genome
from mavis.assemble import Contig
from mavis.bam.cache import BamCache
import mavis.bam.cigar as cigar_tools
from mavis.breakpoint import Breakpoint
from mavis.constants import CIGAR, ORIENT, reverse_complement
from mavis.interval import Interval
from mavis.validate.evidence import GenomeEvidence

from . import BAM_INPUT, MockBamFileHandle, MockObject, MockRead, REFERENCE_GENOME_FILE, REFERENCE_GENOME_FILE_2BIT


REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(REFERENCE_GENOME_FILE)
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(BAM_INPUT)


class TestAlign(unittest.TestCase):
    def setUp(self):
        self.cache = BamCache(MockBamFileHandle({'Y': 23, 'fake': 0, 'reference3': 3}))

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs(self):
        ev = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True,
            bam_cache=None,
            reference_genome=None,
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
        align_contigs([ev], BAM_CACHE, REFERENCE_GENOME, aligner_reference=REFERENCE_GENOME_FILE_2BIT,
                      aligner='blat')
        read1, read2 = ev.contigs[0].alignments[0]
        self.assertEqual(reverse_complement(read1.query_sequence), read2.query_sequence)
        self.assertEqual(1, read1.reference_id)
        self.assertEqual(1, read2.reference_id)
        self.assertEqual(Interval(125, 244), query_coverage_interval(read1))
        self.assertEqual(Interval(117, 244), query_coverage_interval(read2))
        self.assertEqual(1114, read1.reference_start)
        self.assertEqual(2187, read2.reference_start)
        self.assertEqual([(CIGAR.S, 125), (CIGAR.EQ, 120)], read1.cigar)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.EQ, 128)], read2.cigar)

    @unittest.skipIf(not shutil.which('bwa'), 'missing the command')
    def test_bwa_contigs(self):
        ev = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True,
            bam_cache=None,
            reference_genome=None,
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
        align_contigs(
            [ev], BAM_CACHE, REFERENCE_GENOME,
            aligner_reference=REFERENCE_GENOME_FILE,
            aligner='bwa mem',
            aligner_output_file='mem.out',
            aligner_fa_input_file='mem.in.fa'
        )
        print(ev.contigs[0].alignments)
        read1, read2 = ev.contigs[0].alignments[0]
        self.assertEqual(reverse_complement(read1.query_sequence), read2.query_sequence)
        self.assertEqual('reference3', read1.reference_name)
        self.assertEqual('reference3', read2.reference_name)
        self.assertEqual(1, read1.reference_id)
        self.assertEqual(1, read2.reference_id)
        self.assertEqual(Interval(125, 244), query_coverage_interval(read1))
        self.assertEqual(Interval(117, 244), query_coverage_interval(read2))
        self.assertEqual(1114, read1.reference_start)
        self.assertEqual(2187, read2.reference_start)
        self.assertEqual([(CIGAR.S, 125), (CIGAR.EQ, 120)], read1.cigar)
        self.assertEqual([(CIGAR.S, 117), (CIGAR.EQ, 128)], read2.cigar)

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_deletion(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=None,
            reference_genome=None,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100
        )
        ev.contigs = [
            Contig(
                'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT'
                'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT', 0)
        ]
        align_contigs([ev], BAM_CACHE, REFERENCE_GENOME, aligner_reference=REFERENCE_GENOME_FILE_2BIT, aligner='blat')
        read1, read2 = ev.contigs[0].alignments[0]
        self.assertTrue(read2 is None)
        self.assertEqual(0, read1.reference_id)
        self.assertTrue(not read1.is_reverse)
        self.assertEqual(Interval(0, 175), query_coverage_interval(read1))
        self.assertEqual(1612, read1.reference_start)
        self.assertEqual([(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)], read1.cigar)

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_inversion(self):
        raise unittest.SkipTest('TODO')

    @unittest.skipIf(not shutil.which('blat'), 'missing the blat command')
    def test_blat_contigs_deletion_revcomp(self):
        ev = GenomeEvidence(
            Breakpoint('fake', 1714, orient=ORIENT.LEFT),
            Breakpoint('fake', 2968, orient=ORIENT.RIGHT),
            opposing_strands=False,
            bam_cache=None,
            reference_genome=None,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100
        )
        seq = 'GGTATATATTTCTCAGATAAAAGATATTTTCCCTTTTATCTTTCCCTAAGCTCACACTACATATATTGCATTTATCTTATATCTGCTTTAAAACCTATTTAT' \
              'TATGTCATTTAAATATCTAGAAAAGTTATGACTTCACCAGGTATGAAAAATATAAAAAGAACTCTGTCAAGAAT'
        ev.contigs = [Contig(reverse_complement(seq), 0)]
        align_contigs([ev], BAM_CACHE, REFERENCE_GENOME, aligner_reference=REFERENCE_GENOME_FILE_2BIT, aligner='blat')
        read1, read2 = ev.contigs[0].alignments[0]
        self.assertTrue(read2 is None)
        self.assertEqual(0, read1.reference_id)
        self.assertTrue(read1.is_reverse)
        self.assertEqual(seq, read1.query_sequence)
        self.assertEqual(Interval(0, 175), query_coverage_interval(read1))
        self.assertEqual(1612, read1.reference_start)
        self.assertEqual([(CIGAR.EQ, 102), (CIGAR.D, 1253), (CIGAR.EQ, 74)], read1.cigar)


class TestBreakpointContigRemappedDepth(unittest.TestCase):
    def setUp(self):
        self.contig = Contig(' ' * 60, None)
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=10))
        self.contig.add_mapped_sequence(MockObject(reference_start=0, reference_end=20))
        self.contig.add_mapped_sequence(MockObject(reference_start=50, reference_end=60))

    def test_break_left_deletion(self):
        b = Breakpoint('10', 1030, 1030, orient=ORIENT.LEFT)
        read = MockRead(
            cigar=cigar_tools.convert_string_to_cigar('35M10D5I20M'),
            reference_start=999,
            reference_name='10'
        )
        SplitAlignment.breakpoint_contig_remapped_depth(b, self.contig, read)
