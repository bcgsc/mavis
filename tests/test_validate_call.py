import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from mavis.assemble import Contig #Not sure if I need this atm
from mavis.breakpoint import Breakpoint
from mavis.annotate import load_reference_genome, Gene, usTranscript, Transcript
from mavis.constants import ORIENT, STRAND, CIGAR, PYSAM_READ_FLAGS, SVTYPE, CALL_METHOD
from mavis.interval import Interval
from mavis.bam.cache import BamCache
from mavis.bam.read import sequenced_strand
from tests import MockRead, mock_read_pair, MockContig
import unittest
from tests import REFERENCE_GENOME_FILE, BAM_INPUT, FULL_BAM_INPUT, MockBamFileHandle
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence
import mavis.validate.call as call
from mavis.validate.call import EventCall

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
                query_name="test2",
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
        event = EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence, SVTYPE.DEL, CALL_METHOD.SPLIT)

        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # now test one where the read pair type is right but the positioning of the reads doesn't
        # support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 1200, 1260, is_reverse=True)
            ))
        event.pull_flanking_support(flanking_pairs)
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
        event = EventCall(
            Breakpoint('1', 900, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence, SVTYPE.DEL, CALL_METHOD.SPLIT)

        event.pull_flanking_support(flanking_pairs)
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
        event = EventCall(
            Breakpoint('1', 800, orient=ORIENT.LEFT),
            Breakpoint('1', 900, orient=ORIENT.RIGHT),
            evidence, SVTYPE.INS, CALL_METHOD.SPLIT)
        event.pull_flanking_support(flanking_pairs)
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
        event = EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.LEFT),
            evidence, SVTYPE.INV, CALL_METHOD.SPLIT)

        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # test read that is the right type but the positioning does not support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 900,950, is_reverse=True)
            ))
        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_inverted_translocation(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            opposing_strands=True
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1100, 1150, is_reverse=True),
                MockRead('r1', 1, 1200, 1250, is_reverse=True)
            )]
        event = EventCall(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            evidence, SVTYPE.ITRANS, CALL_METHOD.SPLIT)
        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_translocation_RL(self):
        b1 = Breakpoint('11', 128675261, orient=ORIENT.RIGHT, strand=STRAND.POS)
        b2 = Breakpoint('22', 29683123, orient=ORIENT.LEFT, strand=STRAND.POS)
        evidence = self.build_genome_evidence(b1, b2)
        event = EventCall(b1, b2, evidence, SVTYPE.TRANS, CALL_METHOD.CONTIG)
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
        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(len(flanking_pairs), len(event.flanking_pairs))

    def test_translocation_RL_filter_nonsupporting(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('2', 1250, orient=ORIENT.LEFT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1201, 1249, is_reverse=True),
                MockRead('r1', 1, 1201, 1249, is_reverse=False)
            )]
        event = EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('2', 1250, orient=ORIENT.LEFT),
            evidence, SVTYPE.TRANS, CALL_METHOD.SPLIT)

        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

        # test read that is the right type but the positioning does not support the current call
        # the mate is on the wrong chromosome (not sure if this would actually be added as flanking support)
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 1200, 1249, is_reverse=True),
                MockRead('r1', 0, 1201, 1249, is_reverse=False)
            ))
        event.pull_flanking_support(flanking_pairs)
        self.assertEqual(1, len(event.flanking_pairs))

    def test_duplication(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            opposing_strands = False
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1205, 1250, is_reverse=True),
                MockRead('r1', 0, 1260, 1295, is_reverse=False)
            )]
        event = EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            evidence, SVTYPE.DUP, CALL_METHOD.SPLIT)

        event.pull_flanking_support(flanking_pairs)
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
            min_splits_reads_resolution=1
        )
        return evidence

    def test_call_all_methods(self):
        #DEL on 100 - 481 with contig and possible del from 30 - 501, 30 - (691,806), (90,199) - 501 and ins 30 - 501
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        contig = MockContig('',
                            [mock_read_pair(MockRead(query_name='t1', reference_id=0, reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)], query_sequence='A'*100),
                                            MockRead(query_name='t1', reference_id=0, reference_start=460, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)], query_sequence='A'*100))])
        contig.input_reads={MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t2', reference_start=100, cigar=[(CIGAR.EQ, 50), (CIGAR.S, 50)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t5', reference_start=20, cigar=[(CIGAR.EQ, 10), (CIGAR.S, 80)])
        )
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t4', reference_id=0, reference_start=10, reference_end=40, is_reverse=False),
                MockRead(query_name='t4', reference_id=0, reference_start=505, reference_end=540, is_reverse=True)
        ))
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=49, reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=805, reference_end=840, is_reverse=True)
        ))

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(5, len(events))
        self.assertEqual(('contig','contig'), events[0].call_method)
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(481, events[0].break2.start)
        self.assertEqual(('split reads', 'split reads'), events[1].call_method)
        self.assertEqual(30, events[1].break1.start)
        self.assertEqual(501, events[1].break2.start)
        self.assertEqual(('split reads', 'flanking reads'), events[2].call_method)
        self.assertEqual(30, events[2].break1.start)
        self.assertEqual(691, events[2].break2.start)
        self.assertEqual(806, events[2].break2.end)
        self.assertEqual(('flanking reads', 'split reads'), events[3].call_method)
        self.assertEqual(90, events[3].break1.start)
        self.assertEqual(199, events[3].break1.end)
        self.assertEqual(501, events[3].break2.start)
        self.assertEqual(('split reads', 'split reads'), events[4].call_method)
        self.assertEqual(30, events[4].break1.start)
        self.assertEqual(501, events[4].break2.start)

    def test_call_contig_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        contig = MockContig('',
                            [mock_read_pair(MockRead(query_name='t1', reference_id=0, reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)], query_sequence='A'*100),
                                            MockRead(query_name='t1', reference_id=0, reference_start=460, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)], query_sequence='A'*100))])
        contig.input_reads={MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t2', reference_start=100, cigar=[(CIGAR.EQ, 50), (CIGAR.S, 50)])
        )
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=49, reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=505, reference_end=550, is_reverse=True)
        ))

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(481, events[0].break2.start)
        self.assertEqual(('contig','contig'), events[0].call_method)

    def test_call_contig_and_split(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        contig = MockContig('',
                            [mock_read_pair(MockRead(query_name='t1', reference_id=0, reference_start=40, cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)], query_sequence='A'*100),
                                            MockRead(query_name='t1', reference_id=0, reference_start=460, cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)], query_sequence='A'*100))])
        contig.input_reads={MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])}
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t2', reference_start=150, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=520, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t2', reference_start=100, cigar=[(CIGAR.EQ, 50), (CIGAR.S, 50)])
        )
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=49, reference_end=90, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=505, reference_end=550, is_reverse=True)
        ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(3, len(events))
        self.assertEqual(100, events[0].break1.start)
        self.assertEqual(481, events[0].break2.start)
        self.assertEqual(('contig','contig'), events[0].call_method)
        self.assertEqual(('split reads', 'split reads'), events[1].call_method)
        self.assertEqual(170, events[1].break1.start)
        self.assertEqual(521, events[1].break2.start)
        self.assertEqual('insertion', events[2].event_type)
        self.assertEqual(('split reads', 'split reads'), events[2].call_method)
        self.assertEqual(170, events[2].break1.start)
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
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t3', reference_id=0, reference_start=42, reference_end=140, is_reverse=False),
                MockRead(query_name='t3', reference_id=0, reference_start=885, reference_end=905, is_reverse=True)
        ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(2, len(events))
        self.assertEqual(170, events[0].break1.start)
        self.assertEqual(871, events[0].break2.start)
        self.assertEqual(('split reads', 'split reads'), events[0].call_method)
        self.assertEqual(('split reads', 'split reads'), events[1].call_method)
        self.assertEqual(170, events[1].break1.start)
        self.assertEqual(871, events[1].break2.start)
        self.assertEqual('insertion', events[1].event_type)

    def test_call_split_and_flanking(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=140, cigar=[(CIGAR.EQ, 30), (CIGAR.S, 70)])
        )
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t1', reference_id=0, reference_start=42, reference_end=140, is_reverse=False),
                MockRead(query_name='t1', reference_id=0, reference_start=885, reference_end=905, is_reverse=True)
        ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        self.assertEqual(170, events[0].break1.start)
        self.assertEqual(170, events[0].break1.end)
        self.assertEqual(('split reads', 'flanking reads'), events[0].call_method)
        self.assertEqual(756, events[0].break2.start)
        self.assertEqual(886, events[0].break2.end)
        raise unittest.SkipTest('TODO')

    def test_call_flanking_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False
        )
        evidence.flanking_pairs.add(mock_read_pair(
                MockRead(query_name='t1', reference_id=0, reference_start=42, reference_end=140, is_reverse=False),
                MockRead(query_name='t1', reference_id=0, reference_start=885, reference_end=905, is_reverse=True)
        ))
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        self.assertEqual(1, len(events))
        self.assertEqual(140, events[0].break1.start)
        self.assertEqual(192, events[0].break1.end)
        self.assertEqual(('flanking reads', 'flanking reads'), events[0].call_method)
        self.assertEqual(756, events[0].break2.start)
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

        events = call._call_by_supporting_reads(self.ev, SVTYPE.INV)
        self.assertEqual(1, len(events))
        event = events[0]
        self.assertEqual(4, len(event.supporting_reads()))
        self.assertEqual(101, event.break1.start)
        self.assertEqual(101, event.break1.end)
        self.assertEqual(501, event.break2.start)
        self.assertEqual(501, event.break2.end)

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
        self.ev.flanking_pairs.add(mock_read_pair(
            MockRead(query_name='t2', reference_id=0, reference_start=150, reference_end=150),
            MockRead(query_name='t2', reference_id=0, reference_start=505, reference_end=520)
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
        self.ev.flanking_pairs.add(mock_read_pair(
            MockRead(query_name='t2', reference_id=0, reference_start=120, reference_end=140),
            MockRead(query_name='t2', reference_id=0, reference_start=520, reference_end=520)
        ))
        break1, break2 = call._call_by_supporting_reads(self.ev, SVTYPE.INV)[0]

        self.assertEqual(71, break1.start)
        self.assertEqual(121, break1.end)
        self.assertEqual(501, break2.start)
        self.assertEqual(501, break2.end)

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
                query_name="test1", cigar=[(CIGAR.S, 110), (CIGAR.EQ, 40)],
                reference_start=1114, reference_end=1150
            ))
        evidence.split_reads[0].add(
            MockRead(
                query_name="test2", cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)],
                reference_start=1108, reference_end=1115
            ))
        evidence.split_reads[0].add(
            MockRead(
                query_name="test3", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=1114, reference_end=1154,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)]
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name="test4", cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)], reference_start=2187
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name="test5", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)], reference_start=2187
            ))
        evidence.split_reads[1].add(
            MockRead(
                query_name="test1", cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
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


class TestCallByFlankingReadsTranscriptome(unittest.TestCase):
    def build_transcriptome_evidence(self, b1, b2, opposing_strands=False):
        return TranscriptomeEvidence(
            {}, # fake the annotations
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

    def test_call_deletion_evidence_spans_exons(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        t1 = usTranscript([(1001, 1100), (1501, 1700), (2001, 2100), (2201, 2300)], strand='+')
        evidence = self.build_transcriptome_evidence(
            Breakpoint('1', 1051, 1051, 'L', '+'),
            Breakpoint('1', 1551, 1551, 'R', '+')
        )
        #evidence.overlapping_transcripts[0].add(t1)
        #evidence.overlapping_transcripts[1].add(t1)
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
        print('mock read pair', *pair)
        evidence.flanking_pairs.add(pair)
        b1, b2 = call._call_by_flanking_pairs(evidence, SVTYPE.DEL)
        self.assertEqual(Breakpoint('1', 1051, 1250, 'L', '+'), b1)
        self.assertEqual(Breakpoint('1', 2101, 2300, 'R', '+'), b2)

        evidence.flanking_pairs.update({
            mock_read_pair(
                MockRead('name', '1', 1051 - evidence.read_length + 1, 1051, is_reverse=True),
                MockRead('name', '1', 2300, 2300 + evidence.read_length - 1, is_reverse=False)
            )
        })

        #raise unittest.SkipTest('TODO')

if __name__ == "__main__":
    unittest.main()
