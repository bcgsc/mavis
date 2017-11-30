from functools import partial
import unittest

from mavis.annotate.genomic import Gene, Transcript, PreTranscript
from mavis.bam.cache import BamCache
from mavis.bam.read import SamRead
from mavis.bam import cigar as _cigar
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CIGAR, ORIENT, STRAND
from mavis.interval import Interval
from mavis.validate.constants import DEFAULTS
from mavis.validate.base import Evidence
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

from . import mock_read_pair, MockBamFileHandle, MockRead, MockObject


REFERENCE_GENOME = None


class TestDistance(unittest.TestCase):
    def setUp(self):
        self.transcript = PreTranscript([(1001, 1100), (1501, 1600), (2001, 2100), (2201, 2300)], strand='+')
        for patt in self.transcript.generate_splicing_patterns():
            self.transcript.transcripts.append(Transcript(self.transcript, patt))
        self.trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={self.transcript}
        )
        setattr(self.trans_evidence, '_select_transcripts', lambda *pos: self.trans_evidence.overlapping_transcripts)
        setattr(self.trans_evidence, 'distance', partial(TranscriptomeEvidence.distance, self.trans_evidence))

    def test_exonic(self):
        self.assertEqual(Interval(149), self.trans_evidence.distance(1001, 1550))

    def test_intergenic_exonic(self):
        dist = self.trans_evidence.distance(101, 1550)
        self.assertEqual(Interval(1049, 1049), dist)

    def test_intergenic_intergenic(self):
        dist = self.trans_evidence.distance(101, 300)
        self.assertEqual(Interval(199), dist)

    def test_aligned_intronic(self):
        dist = self.trans_evidence.distance(1102, 1499)
        self.assertEqual(Interval(5), dist)

    def test_indel_at_exon_boundary(self):
        self.assertEqual(Interval(2), self.trans_evidence.distance(1101, 1501))

    def test_no_annotations(self):
        dist = self.trans_evidence.distance(101, 300, [])
        self.assertEqual(Interval(199), dist)

    def test_intergenic_intronic(self):
        dist = self.trans_evidence.distance(101, 1400)
        self.assertEqual(Interval(1101), dist)

    def test_empty_intron(self):
        t2 = PreTranscript([(1001, 1100), (1501, 1600), (2001, 2200), (2201, 2300)], strand='+')
        for patt in t2.generate_splicing_patterns():
            t2.transcripts.append(Transcript(t2, patt))
        print(t2)
        print(self.trans_evidence.overlapping_transcripts)
        self.trans_evidence.overlapping_transcripts.add(t2)
        dist = self.trans_evidence.distance(1001, 2301)
        self.assertEqual(Interval(400, 400), dist)


class TestTransStandardize(unittest.TestCase):

    def test_shift_overaligned(self):
        # qwertyuiopas---kkkkk------dfghjklzxcvbnm
        # ..........      ................
        gene = Gene('1', 1, 1000, strand='+')
        transcript = PreTranscript(exons=[(1, 12), (20, 28)], gene=gene, strand='+')
        for spl_patt in transcript.generate_splicing_patterns():
            transcript.transcripts.append(Transcript(transcript, spl_patt))
        gene.transcripts.append(transcript)
        read = SamRead(
            reference_name='1', reference_start=0, cigar=_cigar.convert_string_to_cigar('14=7D12='),
            query_sequence='qwertyuiopasdfghjklzxcvbnm')
        evidence = TranscriptomeEvidence(
            annotations={}, reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnm')},
            bam_cache=MockObject(get_read_reference_name=lambda r: r.reference_name),
            break1=Breakpoint('1', 1, orient='L', strand='+'), break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220)
        evidence.overlapping_transcripts.add(transcript)
        new_read = evidence.standardize_read(read)
        self.assertEqual(_cigar.convert_string_to_cigar('12=7N14='), new_read.cigar)

    def test_shift_overaligned_left(self):
        # qwertyuiopasdf---kkkkkdf------ghjklzxcvbnm
        # ..........      ................
        gene = Gene('1', 1, 1000, strand='+')
        transcript = PreTranscript(exons=[(1, 14), (22, 28)], gene=gene, strand='+')
        for spl_patt in transcript.generate_splicing_patterns():
            transcript.transcripts.append(Transcript(transcript, spl_patt))
        gene.transcripts.append(transcript)
        read = SamRead(
            reference_name='1', reference_start=0, cigar=_cigar.convert_string_to_cigar('12=7D14='),
            query_sequence='qwertyuiopasdfghjklzxcvbnm')
        evidence = TranscriptomeEvidence(
            annotations={}, reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnmsbcdefhi')},
            bam_cache=MockObject(get_read_reference_name=lambda r: r.reference_name),
            break1=Breakpoint('1', 1, orient='L', strand='+'), break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220)
        evidence.overlapping_transcripts.add(transcript)
        new_read = evidence.standardize_read(read)
        self.assertEqual(_cigar.convert_string_to_cigar('14=7N12='), new_read.cigar)

    def test_shift_no_transcripts(self):
        read = SamRead(
            reference_name='1', reference_start=0, cigar=_cigar.convert_string_to_cigar('14=7D18='),
            query_sequence='qwertyuiopasdfdfghjklzxcvbnm')
        evidence = TranscriptomeEvidence(
            annotations={}, reference_genome={'1': MockObject(seq='qwertyuiopasdfkkkkkdfghjklzxcvbnm')},
            bam_cache=None, break1=Breakpoint('1', 1, orient='L', strand='+'), break2=Breakpoint('1', 10, orient='R', strand='+'),
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220)
        new_cigar = evidence.exon_boundary_shift_cigar(read)
        self.assertEqual(_cigar.convert_string_to_cigar('14=7D18='), new_cigar)


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


class TestTraverse(unittest.TestCase):

    def setUp(self):
        self.transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.POS)
        for patt in self.transcript.generate_splicing_patterns():
            self.transcript.transcripts.append(Transcript(self.transcript, patt))

        self.trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={self.transcript}
        )
        setattr(self.trans_evidence, '_select_transcripts', lambda *pos: self.trans_evidence.overlapping_transcripts)
        setattr(self.trans_evidence, 'traverse', partial(TranscriptomeEvidence.traverse, self.trans_evidence))

    def test_left_before_transcript(self):
        exp_pos = Evidence.traverse(900, 500 - 1, ORIENT.LEFT)
        self.assertEqual(exp_pos, self.trans_evidence.traverse(900, 500 - 1, ORIENT.LEFT))

    def test_left_after_transcript(self):
        exp_pos = Evidence.traverse(2200, 100, ORIENT.LEFT)
        self.assertEqual(exp_pos, self.trans_evidence.traverse(2200, 100, ORIENT.LEFT))

    def test_left_at_end(self):
        gpos = self.trans_evidence.traverse(1900, 500, ORIENT.LEFT)
        self.assertEqual(Interval(900), gpos)

    def test_left_within_transcript_exonic(self):
        gpos = self.trans_evidence.traverse(1750, 200 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(1051), gpos)

    def test_left_within_exon(self):
        gpos = self.trans_evidence.traverse(1750, 20 - 1, ORIENT.LEFT)
        self.assertEqual(1731, gpos.start)
        self.assertEqual(1731, gpos.end)

    def test_left_within_transcript_intronic(self):
        gpos = self.trans_evidence.traverse(1600, 150 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(1451), gpos)

    def test_right_before_transcript(self):
        gpos = self.trans_evidence.traverse(500, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(599), gpos)

    def test_right_before_transcript2(self):
        gpos = self.trans_evidence.traverse(901, 500 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1900), gpos)

    def test_right_after_transcript(self):
        gpos = self.trans_evidence.traverse(2201, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(2300), gpos)

    def test_right_within_transcript(self):
        gpos = self.trans_evidence.traverse(1351, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1750), gpos)

    def test_right_within_exon(self):
        gpos = self.trans_evidence.traverse(1351, 10 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1360), gpos)


class TestTraverseTransRev(unittest.TestCase):

    def setUp(self):
        self.transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.NEG)
        for patt in self.transcript.generate_splicing_patterns():
            self.transcript.transcripts.append(Transcript(self.transcript, patt))

        self.trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={self.transcript}
        )
        setattr(self.trans_evidence, '_select_transcripts', lambda *pos: self.trans_evidence.overlapping_transcripts)
        setattr(self.trans_evidence, 'traverse', partial(TranscriptomeEvidence.traverse, self.trans_evidence))

    def test_left_before_transcript(self):
        gpos = self.trans_evidence.traverse(900, 500 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(401), gpos)
        self.assertEqual(gpos, GenomeEvidence.traverse(900, 500 - 1, ORIENT.LEFT))

    def test_left_after_transcript(self):
        gpos = self.trans_evidence.traverse(2200, 100, ORIENT.LEFT)
        self.assertEqual(gpos, GenomeEvidence.traverse(2200, 100, ORIENT.LEFT))
        self.assertEqual(Interval(2100), gpos)

    def test_left_after_transcript2(self):
        gpos = self.trans_evidence.traverse(1900, 500 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(901), gpos)

    def test_left_within_transcript_exonic(self):
        gpos = self.trans_evidence.traverse(1750, 200 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(1051), gpos)

    def test_left_within_exon(self):
        gpos = self.trans_evidence.traverse(1750, 20 - 1, ORIENT.LEFT)
        self.assertEqual(1731, gpos.start)
        self.assertEqual(1731, gpos.end)

    def test_left_within_transcript_intronic(self):
        gpos = self.trans_evidence.traverse(1600, 150 - 1, ORIENT.LEFT)
        self.assertEqual(Interval(1451), gpos)

    def test_right_before_transcript(self):
        gpos = self.trans_evidence.traverse(500, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(599), gpos)

    def test_right_before_transcript2(self):
        gpos = self.trans_evidence.traverse(901, 500 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1900), gpos)

    def test_right_after_transcript(self):
        gpos = self.trans_evidence.traverse(2201, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(2300), gpos)

    def test_right_within_transcript(self):
        gpos = self.trans_evidence.traverse(1351, 100 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1750), gpos)

    def test_right_within_exon(self):
        gpos = self.trans_evidence.traverse(1351, 10 - 1, ORIENT.RIGHT)
        self.assertEqual(Interval(1360), gpos)


class TestTranscriptomeEvidenceWindow(unittest.TestCase):

    def setUp(self):
        gene = Gene('1', 1, 9999, name='KRAS', strand=STRAND.POS)
        self.pre_transcript = PreTranscript(gene=gene, exons=[(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        gene.unspliced_transcripts.append(self.pre_transcript)
        for spl in self.pre_transcript.generate_splicing_patterns():
            self.pre_transcript.transcripts.append(Transcript(self.pre_transcript, spl))
        self.annotations = {gene.chr: [gene]}
        self.genome_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11
        )
        self.trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={self.pre_transcript}
        )
        setattr(self.trans_evidence, '_select_transcripts', lambda *pos: self.trans_evidence.overlapping_transcripts)
        setattr(self.trans_evidence, 'traverse', partial(TranscriptomeEvidence.traverse, self.trans_evidence))

    def transcriptome_window(self, breakpoint, transcripts=None):
        if transcripts:
            self.trans_evidence.overlapping_transcripts.update(transcripts)
        return TranscriptomeEvidence.generate_window(self.trans_evidence, breakpoint)

    def genome_window(self, breakpoint):
        return GenomeEvidence.generate_window(self.genome_evidence, breakpoint)

    def test_before_start(self):
        b = Breakpoint(chr='1', start=100, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

        b = Breakpoint(chr='1', start=500, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_after_end(self):
        b = Breakpoint(chr='1', start=6000, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_exonic_long_exon(self):
        b = Breakpoint(chr='1', start=3200, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_intronic_long_exon(self):
        b = Breakpoint(chr='1', start=2970, orient=ORIENT.RIGHT)
        self.assertEqual(self.genome_window(b), self.transcriptome_window(b))

    def test_intronic_long_intron(self):
        b = Breakpoint(chr='1', start=1800, orient=ORIENT.RIGHT)
        print(self.genome_window(b))
        self.assertEqual(Interval(1490, 2360), self.transcriptome_window(b))

    def test_intronic_short_exon_right(self):
        b = Breakpoint(chr='1', start=1690, orient=ORIENT.RIGHT)
        print(self.genome_window(b))
        self.assertEqual(Interval(1580, 3500), self.transcriptome_window(b))

    def test_intronic_short_exon_left(self):
        b = Breakpoint(chr='1', start=2200, orient=ORIENT.LEFT)
        self.assertEqual(Interval(1440, 2310), self.transcriptome_window(b))

    def test_multiple_transcripts(self):
        #  [(1001, 1100), (1401, 1500), (1701, 1750), (3001, 4000)])
        b = Breakpoint(chr='1', start=1150, orient=ORIENT.RIGHT)
        gene = self.annotations['1'][0]
        t2 = PreTranscript(gene=gene, exons=[(1001, 1100), (1200, 1300), (2100, 2200)])
        for patt in t2.generate_splicing_patterns():
            t2.transcripts.append(Transcript(t2, patt))
        gene.transcripts.append(t2)
        # 989 - 2561
        # 989 - 3411
        self.assertEqual(Interval(1040, 3160), self.transcriptome_window(b, [self.pre_transcript, t2]))

    def test_many_small_exons(self):
        g = Gene('fake', 17271277, 17279592, strand='+')
        pre_transcript = PreTranscript(
            gene=g,
            exons=[
                (17271277, 17271984),
                (17272649, 17272709),
                (17275586, 17275681),
                (17275769, 17275930),
                (17276692, 17276817),
                (17277168, 17277388),  # 220
                (17277845, 17277888),  # 44
                (17278293, 17278378),  # 86
                (17279229, 17279592)  # 364
            ])
        g.transcripts.append(pre_transcript)
        for patt in pre_transcript.generate_splicing_patterns():
            pre_transcript.transcripts.append(Transcript(pre_transcript, patt))
        b = Breakpoint(chr='fake', start=17279591, orient=ORIENT.LEFT)
        self.assertEqual(Interval(17277321, 17279701), self.transcriptome_window(b, [pre_transcript]))


class TestNetSizeTrans(unittest.TestCase):

    def setUp(self):
        self.transcript = PreTranscript([(1001, 1100), (1301, 1400), (1701, 1800)], strand=STRAND.POS)
        for patt in self.transcript.generate_splicing_patterns():
            self.transcript.transcripts.append(Transcript(self.transcript, patt))
        self.trans_evidence = MockObject(
            annotations={},
            read_length=100,
            max_expected_fragment_size=550,
            call_error=11,
            overlapping_transcripts={self.transcript}
        )
        setattr(self.trans_evidence, '_select_transcripts', lambda *pos: self.trans_evidence.overlapping_transcripts)
        setattr(self.trans_evidence, 'distance', partial(TranscriptomeEvidence.distance, self.trans_evidence))

    def test_net_zero(self):
        bpp = BreakpointPair(
            Breakpoint('1', 1099, orient=ORIENT.LEFT),
            Breakpoint('1', 1302, orient=ORIENT.RIGHT),
            untemplated_seq='TT'
        )
        dist = partial(TranscriptomeEvidence.distance, self.trans_evidence)
        self.assertEqual(Interval(-200), bpp.net_size())
        self.assertEqual(Interval(0), bpp.net_size(dist))


class TestGenomeEvidenceWindow(unittest.TestCase):

    def test_orient_ns(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.NS)
        window = GenomeEvidence.generate_window(MockObject(read_length=100, max_expected_fragment_size=550, call_error=11), bpp)
        self.assertEqual(440, window.start)
        self.assertEqual(1560, window.end)
        self.assertEqual(1121, len(window))

    def test_orient_left(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.LEFT)
        window = GenomeEvidence.generate_window(MockObject(read_length=100, call_error=11, max_expected_fragment_size=550), bpp)
        self.assertEqual(440, window.start)
        self.assertEqual(1110, window.end)
        self.assertEqual(671, len(window))

    def test_orient_right(self):
        bpp = Breakpoint(chr='1', start=1000, end=1000, orient=ORIENT.RIGHT)
        window = GenomeEvidence.generate_window(MockObject(read_length=100, call_error=11, max_expected_fragment_size=550), bpp)
        self.assertEqual(890, window.start)
        self.assertEqual(1560, window.end)
        self.assertEqual(671, len(window))

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
