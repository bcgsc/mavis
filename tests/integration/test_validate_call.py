import unittest

from mavis.align import call_paired_read_event, select_contig_alignments
from mavis.annotate.file_io import load_reference_genome
from mavis.annotate.genomic import PreTranscript, Transcript
from mavis.bam.cache import BamCache
from mavis.bam.read import sequenced_strand, SamRead
from mavis.bam.cigar import convert_string_to_cigar
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CALL_METHOD, CIGAR, ORIENT, PYSAM_READ_FLAGS, STRAND, SVTYPE, reverse_complement
from mavis.interval import Interval
from mavis.validate import call
from mavis.validate.base import Evidence
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

from . import BAM_INPUT, FULL_BAM_INPUT, mock_read_pair, MockBamFileHandle, MockObject, MockRead, REFERENCE_GENOME_FILE, get_example_genes, MockLongString

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


class TestCallByContig(unittest.TestCase):
    def test_EGFR_small_del_transcriptome(self):
        gene = get_example_genes()['EGFR']
        reference_annotations = {gene.chr: [gene]}
        reference_genome = {gene.chr: MockObject(
            seq=MockLongString(gene.seq, offset=gene.start - 1)
        )}

        read = SamRead(
            query_sequence='CTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGGATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTACGGGGTGACCGTTTGGGAGTTGATGACCTTTGGATCCAA',
            cigar=convert_string_to_cigar('68M678D50M15D34M6472D185M10240D158M891D74M' '5875D' '6M' '1X' '29M'),
            reference_name='7',
            reference_id=6,
            reference_start=55241669
        )
        print('read.cigar', read.cigar)
        evidence = TranscriptomeEvidence(
            reference_annotations,
            Breakpoint(gene.chr, gene.start, gene.end, orient='L', strand='+'), Breakpoint(gene.chr, gene.start, gene.end, orient='R', strand='+'),
            reference_genome=reference_genome,
            read_length=75, stdev_fragment_size=75, median_fragment_size=220,
            bam_cache=MockObject(get_read_reference_name=lambda x: gene.chr, stranded=True)
        )
        evidence.contigs.append(MockObject(seq=read.query_sequence, alignments=set()))
        select_contig_alignments(evidence, {read.query_sequence: {read}})
        print('distance', evidence.distance(55219055, 55220239))
        print('selected contig alignments')
        print([c.alignments for c in evidence.contigs])
        events = call._call_by_contigs(evidence)
        self.assertEqual(1, len(events))
        self.assertEqual(Breakpoint('7', 55242465, orient='L', strand='+'), events[0].break1)
        self.assertEqual(Breakpoint('7', 55242481, orient='R', strand='+'), events[0].break2)
        print(events[0].contig_alignment.score())
        self.assertTrue(events[0].contig_alignment.score() > 0.99)


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
        self.ev = call.EventCall(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            source_evidence=self.ev1,
            event_type=SVTYPE.INV,
            call_method=CALL_METHOD.SPLIT
        )

    def test_flanking_support_empty(self):
        self.assertEqual(0, len(self.ev.flanking_pairs))

    def test_flanking_support(self):
        # 1114 ++
        # 2187 ++
        self.ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='test1',
                    reference_id=3,
                    template_length=500,
                    reference_start=1150,
                    reference_end=1200,
                    is_reverse=True),
                MockRead(
                    reference_id=3,
                    reference_start=2200,
                    reference_end=2250,
                    is_reverse=True
                )
            ))
        self.ev.flanking_pairs.add(mock_read_pair(
            MockRead(
                query_name='test2',
                reference_id=3,
                template_length=560,
                reference_start=1150,
                reference_end=1200,
                is_reverse=True
            ),
            MockRead(
                reference_id=3,
                reference_start=2200,
                reference_end=2250,
                is_reverse=True
            )
        ))
        median, stdev = self.ev.flanking_metrics()
        self.assertEqual(2, len(self.ev.flanking_pairs))
        self.assertEqual(530, median)
        self.assertEqual(30, stdev)

    def test_split_read_support_empty(self):
        self.assertEqual(0, len(self.ev.break1_split_reads) + len(self.ev.break2_split_reads))

    def test_call_by_split_delins_del_only(self):
        raise unittest.SkipTest('TODO')

    def test_call_by_split_delins_both(self):
        raise unittest.SkipTest('TODO')

    def test_call_by_split_delins_ins_only(self):
        # not implemented yet??
        raise unittest.SkipTest('TODO')


class TestPullFlankingSupport(unittest.TestCase):

    def setUp(self):
        self.bam_cache = BamCache(MockBamFileHandle({'1': 0, '2': 1}))
        self.REFERENCE_GENOME = None

    def build_genome_evidence(self, b1, b2, opposing_strands=False):
        evidence = GenomeEvidence(
            b1, b2, self.bam_cache, self.REFERENCE_GENOME,
            opposing_strands=opposing_strands,
            read_length=100, median_fragment_size=500, stdev_fragment_size=50,
            stdev_count_abnormal=3
        )
        return evidence

    def test_deletion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 1200, 1260, is_reverse=True)
            )]
        event = call.EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence, SVTYPE.DEL, CALL_METHOD.SPLIT)

        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # now test one where the read pair type is right but the positioning of the reads doesn't
        # support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 1200, 1260, is_reverse=True)
            ))
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_small_deletion_flanking_for_larger_deletion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 900, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 1500, 1260, is_reverse=True)
            )]
        event = call.EventCall(
            Breakpoint('1', 900, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence, SVTYPE.DEL, CALL_METHOD.SPLIT)

        event.add_flanking_support(flanking_pairs)
        self.assertEqual(0, len(event.flanking_pairs))

    def test_insertion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 800, orient=ORIENT.LEFT),
            Breakpoint('1', 900, orient=ORIENT.RIGHT))
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 700, 750, is_reverse=False),
                MockRead('r1', 0, 950, 1049, is_reverse=True)
            )]
        print(evidence.min_expected_fragment_size)
        event = call.EventCall(
            Breakpoint('1', 800, orient=ORIENT.LEFT),
            Breakpoint('1', 900, orient=ORIENT.RIGHT),
            evidence, SVTYPE.INS, CALL_METHOD.SPLIT)
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_inversion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.LEFT),
            opposing_strands=True
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 900, 950, is_reverse=False)
            )]
        event = call.EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.LEFT),
            evidence, SVTYPE.INV, CALL_METHOD.SPLIT)

        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # test read that is the right type but the positioning does not support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 900, 950, is_reverse=True)
            ))
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_inverted_translocation(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            opposing_strands=True
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1MockSeq, ', 0, 1100, 1150, is_reverse=True),
                MockRead('r1', 1, 1200, 1250, is_reverse=True)
            )]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            evidence, SVTYPE.ITRANS, CALL_METHOD.SPLIT)
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_translocation_rl(self):
        b1 = Breakpoint('11', 128675261, orient=ORIENT.RIGHT, strand=STRAND.POS)
        b2 = Breakpoint('22', 29683123, orient=ORIENT.LEFT, strand=STRAND.POS)
        evidence = self.build_genome_evidence(b1, b2)
        event = call.EventCall(b1, b2, evidence, SVTYPE.TRANS, CALL_METHOD.CONTIG)
        flanking_pairs = [
            mock_read_pair(
                MockRead('x', '11', 128675264, 128677087, is_reverse=False),
                MockRead('x', '22', 29683030, 29683105, is_reverse=True)
            ),
            mock_read_pair(
                MockRead('x', '11', 128675286, 128677109, is_reverse=False),
                MockRead('x', '22', 29683016, 29683091, is_reverse=True)
            ),
            mock_read_pair(
                MockRead('x', '11', 128675260, 128677083, is_reverse=False),
                MockRead('x', '22', 29683049, 29683123, is_reverse=True)
            ),
            mock_read_pair(
                MockRead('x', '11', 128675289, 128677110, is_reverse=False),
                MockRead('x', '22', 29683047, 29683122, is_reverse=True)
            ),
            mock_read_pair(
                MockRead('x', '11', 128675306, 128677129, is_reverse=False),
                MockRead('x', '22', 29683039, 29683114, is_reverse=True)
            ),
            mock_read_pair(
                MockRead('x', '11', 128675289, 128677110, is_reverse=False),
                MockRead('x', '22', 29683047, 29683122, is_reverse=True)
            )
        ]
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(len(flanking_pairs), len(event.flanking_pairs))

    def test_translocation_rl_filter_nonsupporting(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('2', 1250, orient=ORIENT.LEFT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1201, 1249, is_reverse=True),
                MockRead('r1', 1, 1201, 1249, is_reverse=False)
            )]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('2', 1250, orient=ORIENT.LEFT),
            evidence, SVTYPE.TRANS, CALL_METHOD.SPLIT)

        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # test read that is the right type but the positioning does not support the current call
        # the mate is on the wrong chromosome (not sure if this would actually be added as flanking support)
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 1200, 1249, is_reverse=True),
                MockRead('r1', 0, 1201, 1249, is_reverse=False)
            ))
        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_duplication(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            opposing_strands=False
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1205, 1250, is_reverse=True),
                MockRead('r1', 0, 1260, 1295, is_reverse=False)
            )]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            evidence, SVTYPE.DUP, CALL_METHOD.SPLIT)

        event.add_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_outside_call_range(self):
        raise unittest.SkipTest('TODO')


class TestEvidenceConsumption(unittest.TestCase):

    def setUp(self):
        self.bam_cache = BamCache(MockBamFileHandle({'1': 0, '2': 1}))
        self.REFERENCE_GENOME = None

    def build_genome_evidence(self, b1, b2, opposing_strands=False):
        evidence = GenomeEvidence(
            b1, b2, self.bam_cache, self.REFERENCE_GENOME,
            opposing_strands=opposing_strands,
            read_length=100, median_fragment_size=200, stdev_fragment_size=50,
            stdev_count_abnormal=3, min_flanking_pairs_resolution=1,
            min_splits_reads_resolution=1,
            min_spanning_reads_resolution=3,
            min_linking_split_reads=1
        )
        return evidence

    def test_call_all_methods(self):
        # DEL on 100 - 481 with contig
        # DEL from 120 - 501 with split
        # DEL from 90-299 - 591-806 with flanking
        #
        # and possible del from 30 - 501,
        # 30 - (691,806), (90,199) - 501 and ins 120 - 501
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        r1, r2 = mock_read_pair(
            MockRead(
                query_name='t1', reference_id=0, reference_name='1', reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                query_sequence='A' * 100),
            MockRead(
                query_name='t1', reference_id=0, reference_name='1', reference_start=460, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                query_sequence='A' * 100))
        contig = MockObject(
            seq='',
            alignments=[
                call_paired_read_event(r1, r2)
            ])
        contig.input_reads = {MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t4', reference_id=0, reference_start=10, reference_end=40, is_reverse=False),
                MockRead(query_name='t4', reference_id=0, reference_start=505, reference_end=540, is_reverse=True)
            ))
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=49, reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=805, reference_end=840, is_reverse=True)
            ))

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(4, len(events))
        self.assertEqual('contig', events[0].call_method)
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(481, events[0].break2.start)
        self.assertEqual('deletion', events[0].event_type)
        self.assertEqual('split reads', events[1].call_method)
        self.assertEqual(120, events[1].break1.start)
        self.assertEqual(501, events[1].break2.start)
        self.assertEqual('deletion', events[1].event_type)
        self.assertEqual('flanking reads', events[2].call_method)
        self.assertEqual(90, events[2].break1.start)
        self.assertEqual(299, events[2].break1.end)
        self.assertEqual(591, events[2].break2.start)
        self.assertEqual(806, events[2].break2.end)
        self.assertEqual('deletion', events[2].event_type)
        self.assertEqual('split reads', events[3].call_method)
        self.assertEqual(120, events[3].break1.start)
        self.assertEqual(501, events[3].break2.start)
        self.assertEqual('insertion', events[3].event_type)

    def test_call_contig_only(self):
        # event should only be 100L+, 501R+ deletion
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        r1, r2 = mock_read_pair(
            MockRead(query_name='t1', reference_id=0, reference_name='1', reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                     query_sequence='A' * 100),
            MockRead(query_name='t1', reference_id=0, reference_name='1', reference_start=480, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                     query_sequence='A' * 100))
        bpp = call_paired_read_event(r1, r2)
        contig = MockObject(
            seq='',
            alignments=[bpp])
        contig.input_reads = {MockRead(query_name='t1', reference_start=100, reference_name='1', cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_name='1', reference_start=80, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_name='1', reference_start=500, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t2', reference_name='1', reference_start=40, cigar=[(CIGAR.EQ, 50), (CIGAR.S, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t3', reference_name='1', reference_id=0, reference_start=49, reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_name='1', reference_id=0, reference_start=505, reference_end=550, is_reverse=True)
            ))

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(501, events[0].break2.start)
        self.assertEqual('contig', events[0].call_method)

    def test_call_contig_and_split(self):
        # contig breakpoint is 100L 501R, split reads is 120L 521R
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        r1, r2 = mock_read_pair(
            MockRead(query_name='t1', reference_id=0, reference_name='1', reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                     query_sequence='A' * 100),
            MockRead(query_name='t1', reference_id=0, reference_name='1', reference_start=480, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                     query_sequence='A' * 100))
        contig = MockObject(
            seq='',
            alignments=[call_paired_read_event(r1, r2)])
        contig.input_reads = {MockRead(query_name='t1', reference_name='1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, reference_name='1', cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=520, reference_name='1', cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=49, reference_name='1', reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=505, reference_name='1', reference_end=550, is_reverse=True)
            ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(3, len(events))
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(501, events[0].break2.start)
        self.assertEqual('contig', events[0].call_method)
        self.assertEqual('split reads', events[1].call_method)
        self.assertEqual(120, events[1].break1.start)
        self.assertEqual(521, events[1].break2.start)
        self.assertEqual('insertion', events[2].event_type)
        self.assertEqual('split reads', events[2].call_method)
        self.assertEqual(120, events[2].break1.start)
        self.assertEqual(521, events[2].break2.start)

    def test_call_split_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=140, cigar=[(CIGAR.EQ, 30), (CIGAR.S, 70)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=870, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=42, reference_end=140, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=885, reference_end=905, is_reverse=True)
            ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(2, len(events))
        self.assertEqual(170, events[0].break1.start)
        self.assertEqual(871, events[0].break2.start)
        self.assertEqual('split reads', events[0].call_method)
        self.assertEqual('split reads', events[1].call_method)
        self.assertEqual(170, events[1].break1.start)
        self.assertEqual(871, events[1].break2.start)
        self.assertEqual('insertion', events[1].event_type)

    def test_call_flanking_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t1', reference_id=0, reference_start=42, reference_end=140, is_reverse=False),
                MockRead(query_name='t1', reference_id=0, reference_start=885, reference_end=905, is_reverse=True)
            ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        self.assertEqual(140, events[0].break1.start)
        self.assertEqual(292, events[0].break1.end)
        self.assertEqual('flanking reads', events[0].call_method)
        self.assertEqual(656, events[0].break2.start)
        self.assertEqual(886, events[0].break2.end)


class TestCallBySupportingReads(unittest.TestCase):

    def setUp(self):
        self.ev = GenomeEvidence(
            Breakpoint('fake', 50, 150, orient=ORIENT.RIGHT),
            Breakpoint('fake', 450, 550, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=True,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1,
            min_linking_split_reads=1,
            min_spanning_reads_resolution=3
        )
        self.dup = GenomeEvidence(
            Breakpoint('fake', 50, orient=ORIENT.RIGHT),
            Breakpoint('fake', 90, orient=ORIENT.LEFT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=False,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=100,
            stdev_count_abnormal=2,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1,
            min_linking_split_reads=1,
            min_spanning_reads_resolution=3
        )

    def test_empty(self):
        with self.assertRaises(UserWarning):
            break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

    def test_call_duplication_by_split_reads_error(self):
        self.dup.split_reads[0].add(
            MockRead(query_name='t1', reference_start=30, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 20)])
        )
        self.dup.split_reads[1].add(
            MockRead(query_name='t1', reference_start=90, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )

        with self.assertRaises(UserWarning):
            call._call_by_supporting_reads(self.ev, SVTYPE.DUP)

    def test_call_both_by_split_read(self):
        self.ev.split_reads[0].add(MockRead(
            query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='A' * 40))
        self.ev.split_reads[1].add(MockRead(
            query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='G' * 40))
        self.ev.split_reads[0].add(MockRead(
            query_name='t2', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='C' * 40))
        self.ev.split_reads[1].add(MockRead(
            query_name='t2', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='A' * 40))

        events = call._call_by_supporting_reads(self.ev, SVTYPE.INV)
        self.assertEqual(1, len(events))
        event = events[0]
        self.assertEqual(4, len(event.support()))
        self.assertEqual(101, event.break1.start)
        self.assertEqual(101, event.break1.end)
        self.assertEqual(501, event.break2.start)
        self.assertEqual(501, event.break2.end)

    def test_call_by_split_read_low_resolution(self):
        self.ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='A' * 40)
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='N' * 40)
        )

        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, break1.start)
        self.assertEqual(101, break1.end)
        self.assertEqual(501, break2.start)
        self.assertEqual(501, break2.end)

    def test_call_by_split_read_resolve_untemp(self):
        self.ev.split_reads[0].add(MockRead(
            query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT'))
        self.ev.split_reads[1].add(MockRead(
            query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA', is_reverse=True))

        event = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, event.break1.start)
        self.assertEqual(101, event.break1.end)
        self.assertEqual(501, event.break2.start)
        self.assertEqual(501, event.break2.end)
        self.assertEqual('', event.untemplated_seq)

    def test_call_by_split_read_resolve_untemp_exists(self):
        self.ev.split_reads[0].add(MockRead(
            query_name='t1', reference_start=100, cigar=[(CIGAR.S, 22), (CIGAR.EQ, 18)],
            query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT'))
        self.ev.split_reads[1].add(MockRead(
            query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA', is_reverse=True))

        event = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, event.break1.start)
        self.assertEqual(101, event.break1.end)
        self.assertEqual(501, event.break2.start)
        self.assertEqual(501, event.break2.end)
        self.assertEqual('TA', event.untemplated_seq)

    def test_call_by_split_read_shift_overlap(self):
        self.ev.split_reads[0].add(MockRead(
            query_name='t1', reference_start=100, cigar=[(CIGAR.S, 18), (CIGAR.EQ, 22)],
            query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT'))
        self.ev.split_reads[1].add(MockRead(
            query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
            query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA', is_reverse=True))

        event = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(101, event.break1.start)
        self.assertEqual(101, event.break1.end)
        self.assertEqual(503, event.break2.start)
        self.assertEqual(503, event.break2.end)
        self.assertEqual('', event.untemplated_seq)

    def test_both_by_flanking_pairs(self):
        self.ev.flanking_pairs.add(mock_read_pair(
            MockRead(
                query_name='t1', reference_id=0, reference_start=150, reference_end=150
            ),
            MockRead(
                query_name='t1', reference_id=0, reference_start=500, reference_end=520
            )
        ))
        self.ev.flanking_pairs.add(mock_read_pair(
            MockRead(
                query_name='t2', reference_id=0, reference_start=120, reference_end=140
            ),
            MockRead(
                query_name='t2', reference_id=0, reference_start=520, reference_end=520
            )
        ))
        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]
        # 120-149  ..... 500-519
        # max frag = 150 - 80 = 70
        self.assertEqual(42, break1.start)
        self.assertEqual(121, break1.end)
        self.assertEqual(412, break2.start)  # 70 - 21 = 49
        self.assertEqual(501, break2.end)

    def test_call_both_by_split_reads_multiple_calls(self):
        self.ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='A' * 40)
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='T' * 40)
        )
        self.ev.split_reads[0].add(
            MockRead(query_name='t2', reference_start=110, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='T' * 40)
        )
        self.ev.split_reads[1].add(
            MockRead(query_name='t2', reference_start=520, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)], query_sequence='A' * 40)
        )

        evs = call._call_by_supporting_reads(self.ev, SVTYPE.INV)
        self.assertEqual(2, len(evs))

    def test_call_by_split_reads_consume_flanking(self):
        evidence = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            BAM_CACHE, REFERENCE_GENOME,
            opposing_strands=True,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            stdev_count_abnormal=3,
            min_flanking_pairs_resolution=1,
            min_splits_reads_resolution=1,
            min_linking_split_reads=1
        )
        evidence.split_reads[0].add(
            MockRead(
                query_name='test1', cigar=[(CIGAR.S, 110), (CIGAR.EQ, 40)],
                reference_start=1114, reference_end=1150
            ))
        evidence.split_reads[0].add(
            MockRead(
                query_name='test2', cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)],
                reference_start=1108, reference_end=1115
            ))
        evidence.split_reads[0].add(
            MockRead(
                query_name='test3', cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=1114, reference_end=1154,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)]
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name='test4', cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)], reference_start=2187
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name='test5', cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)], reference_start=2187
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name='test1', cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=2187, reference_end=2307,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)]
            ))

        evidence.flanking_pairs.add(mock_read_pair(
            MockRead(query_name='t1', reference_id=3, reference_start=1200, reference_end=1250, is_reverse=True),
            MockRead(reference_id=3, reference_start=2250, reference_end=2300, is_reverse=True)
        ))

        events = call._call_by_supporting_reads(evidence, event_type=SVTYPE.INV)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        event = events[0]
        self.assertEqual(1, len(event.flanking_pairs))
        self.assertEqual(2, len(event.break1_split_reads))
        self.assertEqual(2, len(event.break2_split_reads))
        b1 = set([read.query_name for read in event.break1_split_reads])
        b2 = set([read.query_name for read in event.break2_split_reads])
        self.assertEqual(1, len(b1 & b2))


class TestCallByFlankingReadsGenome(unittest.TestCase):

    def setUp(self):
        self.ev_LR = GenomeEvidence(
            Breakpoint('fake', 100, orient=ORIENT.LEFT),
            Breakpoint('fake', 200, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=False,
            read_length=25,
            stdev_fragment_size=25,
            median_fragment_size=100,
            stdev_count_abnormal=2,
            min_flanking_pairs_resolution=1
        )

    def test_call_coverage_too_large(self):
        with self.assertRaises(AssertionError):
            call._call_interval_by_flanking_coverage(Interval(1901459, 1902200), ORIENT.RIGHT, 725 + 150, 150, Evidence.distance, Evidence.traverse)

    def test_call_both_intrachromosomal_lr(self):
        # --LLL-100------------500-RRR-------
        # max fragment size: 100 + 2 * 25 = 150
        # max distance = 150 - read_length = 100
        # coverage ranges: 20->80 (61)   600->675 (76)
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
        self.assertEqual(80 + 64, break1.end)
        self.assertEqual(600 - 49, break2.start)
        self.assertEqual(600, break2.end)

    def test_call_both_intrachromosomal_lr_coverage_overlaps_range(self):
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
        self.assertEqual(80, break1.end)  # 119
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
            BamCache(MockBamFileHandle()), None,
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
        b2 = Breakpoint(self.ev_LR.break2.chr, 121, 150, orient=ORIENT.RIGHT)
        break1, break2 = call._call_by_flanking_pairs(
            self.ev_LR, SVTYPE.INV,
            second_breakpoint_called=b2
        )
        BreakpointPair(break1, break2, opposing_strands=False)
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

    def test_call_second_with_first_given_incompatible_error_with_overlap(self):
        evidence = GenomeEvidence(
            Breakpoint('1', 2686252, orient=ORIENT.RIGHT),
            Breakpoint('1', 2686425, 2686667, orient=ORIENT.LEFT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=102,
            median_fragment_size=431,
            min_flanking_pairs_resolution=1
        )
        evidence.flanking_pairs.add((
            MockRead(reference_start=2686251, reference_end=2686329, next_reference_start=2686290),
            MockRead(reference_start=2686290, reference_end=2686367, next_reference_start=2686251)
        ))
        evidence.flanking_pairs.add((
            MockRead(reference_start=2686218, reference_end=2686320, next_reference_start=2686240),
            MockRead(reference_start=2686240, reference_end=2686345, next_reference_start=2686218)
        ))
        with self.assertRaises(AssertionError):
            break1, break2 = call._call_by_flanking_pairs(
                evidence, SVTYPE.DUP, Breakpoint('1', 2686471, orient=ORIENT.RIGHT))

    def test_call_with_overlapping_coverage_intervals(self):
        evidence = GenomeEvidence(
            Breakpoint('1', 76185710, 76186159, orient=ORIENT.RIGHT),
            Breakpoint('1', 76186430, 76186879, orient=ORIENT.LEFT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=98,
            median_fragment_size=433,
            min_flanking_pairs_resolution=1
        )
        evidence.flanking_pairs.add((
            MockRead(reference_start=76186159, reference_end=76186309, next_reference_start=76186000),
            MockRead(reference_start=76186000, reference_end=76186150, next_reference_start=76186159)
        ))
        with self.assertRaises(AssertionError):
            break1, break2 = call._call_by_flanking_pairs(
                evidence, SVTYPE.DUP, Breakpoint('1', 76185557, orient=ORIENT.RIGHT))


class TestCallByFlankingReadsTranscriptome(unittest.TestCase):

    def build_transcriptome_evidence(self, b1, b2, opposing_strands=False):
        return TranscriptomeEvidence(
            {},  # fake the annotations
            b1, b2,
            BamCache(MockBamFileHandle(), stranded=True), None,  # bam_cache and reference_genome
            opposing_strands=opposing_strands,
            stranded=True,
            read_length=50,
            stdev_fragment_size=100,
            median_fragment_size=100,
            stdev_count_abnormal=3,
            min_splits_reads_resolution=1,
            min_flanking_pairs_resolution=1,
            strand_determining_read=2
        )

    def test_call_translocation(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        raise unittest.SkipTest('TODO')

    def test_call_inversion(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        raise unittest.SkipTest('TODO')

    def test_call_inversion_overlapping_breakpoint_calls(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        raise unittest.SkipTest('TODO')

    def test_call_deletion(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        pre_transcript = PreTranscript([(1001, 1100), (1501, 1700), (2001, 2100), (2201, 2300)], strand='+')
        for patt in pre_transcript.generate_splicing_patterns():
            pre_transcript.transcripts.append(Transcript(pre_transcript, patt))
        evidence = self.build_transcriptome_evidence(
            Breakpoint('1', 1051, 1051, 'L', '+'),
            Breakpoint('1', 1551, 1551, 'R', '+')
        )
        # now add the flanking pairs
        pair = mock_read_pair(
            MockRead('name', '1', 951, 1051, is_reverse=True),
            MockRead('name', '1', 2299, 2399, is_reverse=False)
        )
        # following help in debugging the mockup
        self.assertTrue(pair[0].is_reverse)
        self.assertTrue(pair[0].is_read1)
        self.assertFalse(pair[0].is_read2)
        self.assertFalse(pair[1].is_reverse)
        self.assertFalse(pair[1].is_read1)
        self.assertTrue(pair[1].is_read2)
        self.assertEqual(STRAND.POS, sequenced_strand(pair[0], 2))
        self.assertEqual(STRAND.POS, evidence.decide_sequenced_strand([pair[0]]))
        self.assertEqual(STRAND.POS, sequenced_strand(pair[1], 2))
        self.assertEqual(STRAND.POS, evidence.decide_sequenced_strand([pair[1]]))
        print(evidence.max_expected_fragment_size, evidence.read_length)
        evidence.flanking_pairs.add(pair)
        breakpoint1, breakpoint2 = call._call_by_flanking_pairs(evidence, SVTYPE.DEL)
        self.assertEqual(Breakpoint('1', 1051, 1301, 'L', '+'), breakpoint1)
        self.assertEqual(Breakpoint('1', 2050, 2300, 'R', '+'), breakpoint2)

        # now add the transcript and call again
        evidence.overlapping_transcripts.add(pre_transcript)
        breakpoint1, breakpoint2 = call._call_by_flanking_pairs(evidence, SVTYPE.DEL)
        self.assertEqual(Breakpoint('1', 1051, 2001, 'L', '+'), breakpoint1)
        self.assertEqual(Breakpoint('1', 1650, 2300, 'R', '+'), breakpoint2)


class TestCallBySpanningReads(unittest.TestCase):

    def test_deletion(self):
        # ATCGATCTAGATCTAGGATAGTTCTAGCAGTCATAGCTAT
        ev = GenomeEvidence(
            Breakpoint('fake', 60, orient=ORIENT.LEFT),
            Breakpoint('fake', 70, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle()), None,
            opposing_strands=False,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=180,
            min_flanking_pairs_resolution=1,
            min_spanning_reads_resolution=1
        )
        print(ev.outer_window1, ev.outer_window2)
        spanning_reads = [
            SamRead(
                query_name='name', reference_name='fake', reference_start=50,
                cigar=[(CIGAR.EQ, 15), (CIGAR.D, 5), (CIGAR.I, 2), (CIGAR.EQ, 10)],
                query_sequence='ATCGATCTAGATCTA' 'GG' 'ATAGTTCTAG'),
            SamRead(
                query_name='name', reference_name='fake', reference_start=50,
                cigar=[(CIGAR.EQ, 15), (CIGAR.I, 2), (CIGAR.D, 5), (CIGAR.EQ, 10)],
                query_sequence='ATCGATCTAGATCTA' 'GG' 'ATAGTTCTAG')
        ]
        ev.spanning_reads = set(spanning_reads)
        calls = call._call_by_spanning_reads(ev, set())
        self.assertEqual(1, len(calls))
        self.assertEqual(2, len(calls[0].support()))

    def test_insertion(self):
        pass

    def test_indel(self):
        pass

    def test_inversion(self):
        pass

    def test_duplication(self):
        pass


class TestCharacterizeRepeatRegion(unittest.TestCase):

    def test_bad_deletion_call(self):
        reference_genome = {'19': MockObject(seq=MockLongString(
            'AAATCTTTTTTCCATTATGGCTATACAAAGTGAATACATTTCCACAAGCAAATATGATAGATTAATTGGTGCATTGTATATATTTCTCAAACCATCAGCTCCTCTT'
            'TTTTTCAAAGTCTAGAATTTGTAATGGTGGATATCTCTGTTCTGTATTCTGTTGTCTAGATATCCAAGTTTAATGCAAAATTTTATGACATGGAACTTGACACTTT'
            'CTAGAAATGTTCACATATGGTTGTTTATTAAATTATCTCTCATGGAAATATTTAAATGACATGTTTATTGTCTGAAAAGGACAGATATTTAAGCTTTTTTTTTTTT'
            'TTTTTCTTTTTTTTTGAGAAAGAGTCTCGTTCTTTTGCCCAGGCTGGAGTGCAGTGGTACAATCTTGGCTCACTACAACCTTCACCTCGCAGGTTCAAGCGATTCT'
            'CCTGCCTCAGCCTCCCTAGTAGCTGGGATTACAGGTACACACCACCAGGCCCATCTAATTTTTCTATATTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTG'
            'TTCTCAAACTCCTGACCTCAGCCAATCCGCCCGCCTCAACCTCTTAAAGTGATGGGATTACAGGTGTGAGCCATTGTGCTTGGCCCCCTTTAACTATTTTATGTGA'
            'CTCTTCT', offset=23951075))}
        bpp = BreakpointPair(
            Breakpoint('19', 23951407, orient=ORIENT.LEFT),
            Breakpoint('19', 23951408, orient=ORIENT.RIGHT),
            opposing_strands=False,
            untemplated_seq='',
            event_type=SVTYPE.DEL
        )
        self.assertEqual(0, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_homopolymer_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 120, orient=ORIENT.LEFT),
            Breakpoint('1', 121, orient=ORIENT.RIGHT),
            untemplated_seq='T',
            opposing_strands=False,
            event_type=SVTYPE.INS
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATTTTGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(4, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_homopolymer_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 120, orient=ORIENT.LEFT),
            Breakpoint('1', 122, orient=ORIENT.RIGHT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DEL
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATTTTTGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(4, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_homopolymer_duplication(self):
        bpp = BreakpointPair(
            Breakpoint('1', 121, orient=ORIENT.RIGHT),
            Breakpoint('1', 121, orient=ORIENT.LEFT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DUP
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATTTTTGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(4, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_repeat_duplication(self):
        bpp = BreakpointPair(
            Breakpoint('1', 123, orient=ORIENT.RIGHT),
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DUP
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(2, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_repeat_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TAG',
            opposing_strands=False,
            event_type=SVTYPE.INS
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(3, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_repeat_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 129, orient=ORIENT.RIGHT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DEL
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATAGTAGTAGTAGGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the second breakpoint', reference_genome['1'].seq[bpp.break2.start - 10:bpp.break2.start])
        self.assertEqual(3, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_norepeat_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TTG',
            opposing_strands=False,
            event_type=SVTYPE.INS
        )
        reference_genome = {'1': MockObject(seq=MockLongString(
            'TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100
        ))}
        print('upto and including the first breakpoint', reference_genome['1'].seq[bpp.break1.start - 10:bpp.break1.start])
        self.assertEqual(0, call.EventCall.characterize_repeat_region(bpp, reference_genome))

    def test_invalid_event_type(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.RIGHT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TTG',
            event_type=SVTYPE.INV
        )
        with self.assertRaises(ValueError):
            call.EventCall.characterize_repeat_region(bpp, None)

if __name__ == '__main__':
    unittest.main()
