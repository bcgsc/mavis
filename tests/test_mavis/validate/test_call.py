from unittest import mock

import pytest
from mavis.annotate.file_io import load_reference_genome
from mavis.annotate.genomic import PreTranscript, Transcript
from mavis.bam import cigar as _cigar
from mavis.bam.cache import BamCache
from mavis.bam.cigar import convert_string_to_cigar
from mavis.bam.read import SamRead, read_pair_type, sequenced_strand
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CALL_METHOD, CIGAR, ORIENT, PYSAM_READ_FLAGS, STRAND, SVTYPE
from mavis.interval import Interval
from mavis.validate import call
from mavis.validate.align import call_paired_read_event, select_contig_alignments
from mavis.validate.base import Evidence
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence

from ...util import get_data, todo
from ..mock import MockBamFileHandle, MockLongString, MockRead, get_example_genes, mock_read_pair

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
    global FULL_BAM_CACHE
    FULL_BAM_CACHE = BamCache(get_data('mock_reads_for_events.sorted.bam'))
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
    patcher = mock.patch('mavis.breakpoint.BreakpointPair.untemplated_shift', new=lambda *x: (0, 0))
    patcher.start()


def tearDownModule():
    mock.patch.stopall()


class TestCallByContig:
    def test_EGFR_small_del_transcriptome(self):
        gene = get_example_genes()['EGFR']
        reference_annotations = {gene.chr: [gene]}
        reference_genome = {
            gene.chr: mock.Mock(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }

        read = SamRead(
            query_sequence='CTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCCGGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGACATCTCCGAAAGCCAACAAGGAAATCCTCGATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCCAGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAACTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGGATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTACGGGGTGACCGTTTGGGAGTTGATGACCTTTGGATCCAA',
            cigar=convert_string_to_cigar(
                '68M678D50M15D34M6472D185M10240D158M891D74M' '5875D' '6M' '1X' '29M'
            ),
            reference_name='7',
            reference_id=6,
            reference_start=55241669,
            alignment_rank=0,
        )
        print('read.cigar', _cigar.convert_cigar_to_string(read.cigar))
        evidence = TranscriptomeEvidence(
            reference_annotations,
            Breakpoint(gene.chr, gene.start, gene.end, orient='L', strand='+'),
            Breakpoint(gene.chr, gene.start, gene.end, orient='R', strand='+'),
            reference_genome=reference_genome,
            read_length=75,
            stdev_fragment_size=75,
            median_fragment_size=220,
            bam_cache=mock.Mock(get_read_reference_name=lambda x: gene.chr, stranded=True),
        )
        evidence.contigs.append(mock.Mock(seq=read.query_sequence, alignments=set()))
        select_contig_alignments(evidence, {read.query_sequence: {read}})
        print('distance', evidence.distance(55219055, 55220239))
        print('selected contig alignments')
        for contig in evidence.contigs:
            print(contig)
            for aln in contig.alignments:
                print(aln.alignment_id())
        events = call._call_by_contigs(evidence)
        for ev in events:
            print(ev)
            print(evidence.distance(ev.break1.start, ev.break2.start))
        assert len(events) == 1
        assert events[0].break1 == Breakpoint('7', 55242465, orient='L', strand='+')
        assert events[0].break2 == Breakpoint('7', 55242481, orient='R', strand='+')
        print(events[0].contig_alignment.score())
        assert events[0].contig_alignment.score() > 0.99


class TestEventCall:
    def test_bad_deletion(self):
        evidence = GenomeEvidence(
            Breakpoint('reference3', 16, orient='L'),
            Breakpoint('reference3', 90, orient='R'),
            BAM_CACHE,
            REFERENCE_GENOME,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
        )
        with pytest.raises(ValueError):
            call.EventCall(
                Breakpoint('reference3', 43, orient='L'),
                Breakpoint('reference3', 44, orient='R'),
                evidence,
                event_type=SVTYPE.DEL,
                call_method=CALL_METHOD.SPLIT,
            )

    def test_flanking_support_empty(self):

        ev = call.EventCall(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            source_evidence=GenomeEvidence(
                Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
                Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
                BAM_CACHE,
                REFERENCE_GENOME,
                opposing_strands=True,
                read_length=125,
                stdev_fragment_size=100,
                median_fragment_size=380,
                stdev_count_abnormal=3,
                min_flanking_pairs_resolution=3,
            ),
            event_type=SVTYPE.INV,
            call_method=CALL_METHOD.SPLIT,
        )
        assert len(ev.flanking_pairs) == 0

    def test_flanking_support(self):
        # 1114 ++
        # 2187 ++
        ev = call.EventCall(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            source_evidence=GenomeEvidence(
                Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
                Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
                BAM_CACHE,
                REFERENCE_GENOME,
                opposing_strands=True,
                read_length=125,
                stdev_fragment_size=100,
                median_fragment_size=380,
                stdev_count_abnormal=3,
                min_flanking_pairs_resolution=3,
            ),
            event_type=SVTYPE.INV,
            call_method=CALL_METHOD.SPLIT,
        )
        ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='test1',
                    reference_id=3,
                    template_length=500,
                    reference_start=1150,
                    reference_end=1200,
                    is_reverse=True,
                ),
                MockRead(reference_id=3, reference_start=2200, reference_end=2250, is_reverse=True),
            )
        )
        ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='test2',
                    reference_id=3,
                    template_length=560,
                    reference_start=1150,
                    reference_end=1200,
                    is_reverse=True,
                ),
                MockRead(reference_id=3, reference_start=2200, reference_end=2250, is_reverse=True),
            )
        )
        median, stdev = ev.flanking_metrics()
        assert len(ev.flanking_pairs) == 2
        assert median == 530
        assert stdev == 30

    def test_split_read_support_empty(self):
        ev = call.EventCall(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            source_evidence=GenomeEvidence(
                Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
                Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
                BAM_CACHE,
                REFERENCE_GENOME,
                opposing_strands=True,
                read_length=125,
                stdev_fragment_size=100,
                median_fragment_size=380,
                stdev_count_abnormal=3,
                min_flanking_pairs_resolution=3,
            ),
            event_type=SVTYPE.INV,
            call_method=CALL_METHOD.SPLIT,
        )
        assert len(ev.break1_split_reads) + len(ev.break2_split_reads) == 0

    @todo
    def test_call_by_split_delins_del_only(self):
        pass

    @todo
    def test_call_by_split_delins_both(self):
        pass

    @todo
    def test_call_by_split_delins_ins_only(self):
        # not implemented yet??
        pass


class TestPullFlankingSupport:
    def build_genome_evidence(self, b1, b2, opposing_strands=False):
        evidence = GenomeEvidence(
            b1,
            b2,
            BamCache(MockBamFileHandle({'1': 0, '2': 1})),
            None,
            opposing_strands=opposing_strands,
            read_length=100,
            median_fragment_size=500,
            stdev_fragment_size=50,
            stdev_count_abnormal=3,
        )
        return evidence

    def test_deletion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 500, orient=ORIENT.LEFT), Breakpoint('1', 1000, orient=ORIENT.RIGHT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 1200, 1260, is_reverse=True),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence,
            SVTYPE.DEL,
            CALL_METHOD.SPLIT,
        )

        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

        # now test one where the read pair type is right but the positioning of the reads doesn't
        # support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 1200, 1260, is_reverse=True),
            )
        )
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    def test_small_deletion_flanking_for_larger_deletion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 900, orient=ORIENT.LEFT), Breakpoint('1', 1000, orient=ORIENT.RIGHT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 1500, 1260, is_reverse=True),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 900, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.RIGHT),
            evidence,
            SVTYPE.DEL,
            CALL_METHOD.SPLIT,
        )

        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 0

    def test_insertion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 800, orient=ORIENT.LEFT), Breakpoint('1', 900, orient=ORIENT.RIGHT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 700, 750, is_reverse=False),
                MockRead('r1', 0, 950, 1049, is_reverse=True),
            )
        ]
        print(evidence.min_expected_fragment_size)
        event = call.EventCall(
            Breakpoint('1', 800, orient=ORIENT.LEFT),
            Breakpoint('1', 900, orient=ORIENT.RIGHT),
            evidence,
            SVTYPE.INS,
            CALL_METHOD.SPLIT,
        )
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    def test_inversion(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.LEFT),
            opposing_strands=True,
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 400, 450, is_reverse=False),
                MockRead('r1', 0, 900, 950, is_reverse=False),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 500, orient=ORIENT.LEFT),
            Breakpoint('1', 1000, orient=ORIENT.LEFT),
            evidence,
            SVTYPE.INV,
            CALL_METHOD.SPLIT,
        )

        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

        # test read that is the right type but the positioning does not support the current call
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 501, 600, is_reverse=False),
                MockRead('r1', 0, 900, 950, is_reverse=True),
            )
        )
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    def test_inverted_translocation(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            opposing_strands=True,
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1MockSeq, ', 0, 1100, 1150, is_reverse=True),
                MockRead('r1', 1, 1200, 1250, is_reverse=True),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.LEFT),
            Breakpoint('2', 1300, orient=ORIENT.LEFT),
            evidence,
            SVTYPE.ITRANS,
            CALL_METHOD.SPLIT,
        )
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    def test_translocation_rl(self):
        b1 = Breakpoint('11', 128675261, orient=ORIENT.RIGHT, strand=STRAND.POS)
        b2 = Breakpoint('22', 29683123, orient=ORIENT.LEFT, strand=STRAND.POS)
        evidence = self.build_genome_evidence(b1, b2)
        event = call.EventCall(b1, b2, evidence, SVTYPE.TRANS, CALL_METHOD.CONTIG)
        flanking_pairs = [
            mock_read_pair(
                MockRead('x', '11', 128675264, 128677087, is_reverse=False),
                MockRead('x', '22', 29683030, 29683105, is_reverse=True),
            ),
            mock_read_pair(
                MockRead('x', '11', 128675286, 128677109, is_reverse=False),
                MockRead('x', '22', 29683016, 29683091, is_reverse=True),
            ),
            mock_read_pair(
                MockRead('x', '11', 128675260, 128677083, is_reverse=False),
                MockRead('x', '22', 29683049, 29683123, is_reverse=True),
            ),
            mock_read_pair(
                MockRead('x', '11', 128675289, 128677110, is_reverse=False),
                MockRead('x', '22', 29683047, 29683122, is_reverse=True),
            ),
            mock_read_pair(
                MockRead('x', '11', 128675306, 128677129, is_reverse=False),
                MockRead('x', '22', 29683039, 29683114, is_reverse=True),
            ),
            mock_read_pair(
                MockRead('x', '11', 128675289, 128677110, is_reverse=False),
                MockRead('x', '22', 29683047, 29683122, is_reverse=True),
            ),
        ]
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == len(flanking_pairs)

    def test_translocation_rl_filter_nonsupporting(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT), Breakpoint('2', 1250, orient=ORIENT.LEFT)
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1201, 1249, is_reverse=True),
                MockRead('r1', 1, 1201, 1249, is_reverse=False),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('2', 1250, orient=ORIENT.LEFT),
            evidence,
            SVTYPE.TRANS,
            CALL_METHOD.SPLIT,
        )

        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

        # test read that is the right type but the positioning does not support the current call
        # the mate is on the wrong chromosome (not sure if this would actually be added as flanking support)
        flanking_pairs.append(
            mock_read_pair(
                MockRead('r1', 0, 1200, 1249, is_reverse=True),
                MockRead('r1', 0, 1201, 1249, is_reverse=False),
            )
        )
        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    def test_duplication(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            opposing_strands=False,
        )
        flanking_pairs = [
            mock_read_pair(
                MockRead('r1', 0, 1205, 1250, is_reverse=True),
                MockRead('r1', 0, 1260, 1295, is_reverse=False),
            )
        ]
        event = call.EventCall(
            Breakpoint('1', 1200, orient=ORIENT.RIGHT),
            Breakpoint('1', 1300, orient=ORIENT.LEFT),
            evidence,
            SVTYPE.DUP,
            CALL_METHOD.SPLIT,
        )

        event.add_flanking_support(flanking_pairs)
        assert len(event.flanking_pairs) == 1

    @todo
    def test_outside_call_range(self):
        pass


class TestEvidenceConsumption:
    def build_genome_evidence(self, b1, b2, opposing_strands=False):
        evidence = GenomeEvidence(
            b1,
            b2,
            BamCache(MockBamFileHandle({'1': 0, '2': 1})),
            None,
            opposing_strands=opposing_strands,
            read_length=100,
            median_fragment_size=200,
            stdev_fragment_size=50,
            config={
                'validate.stdev_count_abnormal': 3,
                'validate.min_flanking_pairs_resolution': 1,
                'validate.min_splits_reads_resolution': 1,
                'validate.min_spanning_reads_resolution': 3,
                'validate.min_linking_split_reads': 1,
                'validate.min_call_complexity': 0,
            },
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
            opposing_strands=False,
        )
        r1, r2 = mock_read_pair(
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=40,
                cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                query_sequence='A' * 100,
            ),
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=460,
                cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                query_sequence='A' * 100,
            ),
        )
        contig = mock.Mock(
            **{
                'seq': '',
                'complexity.return_value': 1,
                'alignments': [call_paired_read_event(r1, r2)],
            }
        )
        contig.input_reads = {
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        }
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=100, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=500, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t4',
                    reference_id=0,
                    reference_start=10,
                    reference_end=40,
                    is_reverse=False,
                    query_alignment_length=30,
                ),
                MockRead(
                    query_name='t4',
                    reference_id=0,
                    reference_start=505,
                    reference_end=540,
                    is_reverse=True,
                    query_alignment_length=35,
                ),
            )
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=49,
                    reference_end=90,
                    is_reverse=False,
                    query_alignment_length=41,
                ),
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=805,
                    reference_end=840,
                    is_reverse=True,
                    query_alignment_length=35,
                ),
            )
        )

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 4
        assert events[0].call_method == 'contig'
        assert events[0].break1.start == 100
        assert events[0].break2.start == 481
        assert events[0].event_type == 'deletion'
        assert events[1].call_method == 'split reads'
        assert events[1].break1.start == 120
        assert events[1].break2.start == 501
        assert events[1].event_type == 'deletion'
        assert events[2].call_method == 'flanking reads'
        assert events[2].break1.start == 90
        assert events[2].break1.end == 299
        assert events[2].break2.start == 591
        assert events[2].break2.end == 806
        assert events[2].event_type == 'deletion'
        assert events[3].call_method == 'split reads'
        assert events[3].break1.start == 120
        assert events[3].break2.start == 501
        assert events[3].event_type == 'insertion'

    def test_call_contig_only(self):
        # event should only be 100L+, 501R+ deletion
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        r1, r2 = mock_read_pair(
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=40,
                cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                query_sequence='A' * 100,
                query_alignment_length=100,
            ),
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=480,
                cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                query_sequence='A' * 100,
                query_alignment_length=100,
            ),
        )
        bpp = call_paired_read_event(r1, r2)
        contig = mock.Mock(**{'seq': '', 'complexity.return_value': 1, 'alignments': [bpp]})
        contig.input_reads = {
            MockRead(
                query_name='t1',
                reference_start=100,
                reference_name='1',
                cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)],
            )
        }
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_name='1',
                reference_start=80,
                cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)],
            )
        )
        evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_name='1',
                reference_start=500,
                cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)],
            )
        )
        evidence.split_reads[0].add(
            MockRead(
                query_name='t2',
                reference_name='1',
                reference_start=40,
                cigar=[(CIGAR.EQ, 50), (CIGAR.S, 50)],
            )
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t3',
                    reference_name='1',
                    reference_id=0,
                    reference_start=49,
                    reference_end=90,
                    is_reverse=False,
                    query_alignment_length=100,
                ),
                MockRead(
                    query_name='t3',
                    reference_name='1',
                    reference_id=0,
                    reference_start=505,
                    reference_end=550,
                    is_reverse=True,
                    query_alignment_length=100,
                ),
            )
        )

        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 1
        assert events[0].break1.start == 100
        assert events[0].break2.start == 501
        assert events[0].call_method == 'contig'

    def test_call_contig_and_split(self):
        # contig breakpoint is 100L 501R, split reads is 120L 521R
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 450, 500, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        r1, r2 = mock_read_pair(
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=40,
                cigar=[(CIGAR.EQ, 60), (CIGAR.S, 40)],
                query_sequence='A' * 100,
            ),
            MockRead(
                query_name='t1',
                reference_id=0,
                reference_name='1',
                reference_start=480,
                cigar=[(CIGAR.S, 40), (CIGAR.EQ, 60)],
                query_sequence='A' * 100,
            ),
        )
        contig = mock.Mock(
            **{
                'seq': '',
                'complexity.return_value': 1,
                'alignments': [call_paired_read_event(r1, r2)],
            }
        )
        contig.input_reads = {
            MockRead(
                query_name='t1',
                reference_name='1',
                reference_start=100,
                cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)],
            )
        }
        evidence.contigs.append(contig)

        evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                reference_name='1',
                cigar=[(CIGAR.EQ, 20), (CIGAR.S, 80)],
            )
        )
        evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=520,
                reference_name='1',
                cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)],
            )
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=49,
                    reference_name='1',
                    reference_end=90,
                    is_reverse=False,
                    query_alignment_length=100,
                ),
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=505,
                    reference_name='1',
                    reference_end=550,
                    is_reverse=True,
                    query_alignment_length=100,
                ),
            )
        )
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 3
        assert events[0].break1.start == 100
        assert events[0].break2.start == 501
        assert events[0].call_method == 'contig'
        assert events[1].call_method == 'split reads'
        assert events[1].break1.start == 120
        assert events[1].break2.start == 521
        assert events[2].event_type == 'insertion'
        assert events[2].call_method == 'split reads'
        assert events[2].break1.start == 120
        assert events[2].break2.start == 521

    def test_call_split_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        evidence.split_reads[0].add(
            MockRead(query_name='t1', reference_start=140, cigar=[(CIGAR.EQ, 30), (CIGAR.S, 70)])
        )
        evidence.split_reads[1].add(
            MockRead(query_name='t1', reference_start=870, cigar=[(CIGAR.S, 50), (CIGAR.EQ, 50)])
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=42,
                    reference_end=140,
                    is_reverse=False,
                    query_alignment_length=100,
                ),
                MockRead(
                    query_name='t3',
                    reference_id=0,
                    reference_start=885,
                    reference_end=905,
                    is_reverse=True,
                    query_alignment_length=100,
                ),
            )
        )
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 2
        assert events[0].break1.start == 170
        assert events[0].break2.start == 871
        assert events[0].call_method == 'split reads'
        assert events[1].call_method == 'split reads'
        assert events[1].break1.start == 170
        assert events[1].break2.start == 871
        assert events[1].event_type == 'insertion'

    def test_call_flanking_only(self):
        evidence = self.build_genome_evidence(
            Breakpoint('1', 50, 150, orient=ORIENT.LEFT),
            Breakpoint('1', 850, 900, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t1',
                    reference_id=0,
                    reference_start=42,
                    reference_end=140,
                    is_reverse=False,
                    query_alignment_length=98,
                ),
                MockRead(
                    query_name='t1',
                    reference_id=0,
                    reference_start=885,
                    reference_end=905,
                    is_reverse=True,
                    query_alignment_length=20,
                ),
            )
        )
        events = call.call_events(evidence)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 1
        assert events[0].break1.start == 140
        assert events[0].break1.end == 292
        assert events[0].call_method == 'flanking reads'
        assert events[0].break2.start == 656
        assert events[0].break2.end == 886


@pytest.fixture
def duplication_ev():
    return GenomeEvidence(
        Breakpoint('fake', 50, orient=ORIENT.RIGHT),
        Breakpoint('fake', 90, orient=ORIENT.LEFT),
        BamCache(MockBamFileHandle()),
        None,
        opposing_strands=False,
        read_length=40,
        stdev_fragment_size=25,
        median_fragment_size=100,
        config={
            'validate.stdev_count_abnormal': 2,
            'validate.min_splits_reads_resolution': 1,
            'validate.min_flanking_pairs_resolution': 1,
            'validate.min_linking_split_reads': 1,
            'validate.min_spanning_reads_resolution': 3,
            'validate. min_call_complexity': 0,
        },
    )


@pytest.fixture
def inversion_evidence():
    return GenomeEvidence(
        Breakpoint('fake', 50, 150, orient=ORIENT.RIGHT),
        Breakpoint('fake', 450, 550, orient=ORIENT.RIGHT),
        BamCache(MockBamFileHandle()),
        None,
        opposing_strands=True,
        read_length=40,
        stdev_fragment_size=25,
        median_fragment_size=100,
        config={
            'validate.stdev_count_abnormal': 2,
            'validate.min_splits_reads_resolution': 1,
            'validate.min_flanking_pairs_resolution': 1,
            'validate.min_linking_split_reads': 1,
            'validate.min_spanning_reads_resolution': 3,
            'validate. min_call_complexity': 0,
        },
    )


class TestCallBySupportingReads:
    def test_empty(self, inversion_evidence):
        with pytest.raises(AssertionError):
            call._call_by_flanking_pairs(inversion_evidence, SVTYPE.INV)[0]

    def test_call_no_duplication_by_split_reads(self, duplication_ev, inversion_evidence):
        duplication_ev.split_reads[0].add(
            MockRead(query_name='t1', reference_start=30, cigar=[(CIGAR.EQ, 20), (CIGAR.S, 20)])
        )
        duplication_ev.split_reads[1].add(
            MockRead(query_name='t1', reference_start=90, cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)])
        )

        bpps = call._call_by_split_reads(inversion_evidence, SVTYPE.DUP)
        assert len(bpps) == 0

    def test_by_split_read(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='A' * 40,
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='G' * 40,
            )
        )
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t2',
                reference_start=100,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='C' * 40,
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t2',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='A' * 40,
            )
        )

        events = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)
        assert len(events) == 1
        event = events[0]
        assert len(event.support()) == 4
        assert event.break1.start == 101
        assert event.break1.end == 101
        assert event.break2.start == 501
        assert event.break2.end == 501

    def test_call_by_split_read_low_resolution(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='A' * 40,
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='N' * 40,
            )
        )

        bpp = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)
        assert len(bpp) == 1
        bpp = bpp[0]

        assert bpp.break1.start == 101
        assert bpp.break1.end == 101
        assert bpp.break2.start == 501
        assert bpp.break2.end == 501

    def test_call_by_split_read_resolve_untemp(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT',
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA',
                is_reverse=True,
            )
        )

        event = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)[0]

        assert event.break1.start == 101
        assert event.break1.end == 101
        assert event.break2.start == 501
        assert event.break2.end == 501
        assert event.untemplated_seq == ''

    def test_call_by_split_read_resolve_untemp_exists(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 22), (CIGAR.EQ, 18)],
                query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT',
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA',
                is_reverse=True,
            )
        )

        event = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)[0]

        assert event.break1.start == 101
        assert event.break1.end == 101
        assert event.break2.start == 501
        assert event.break2.end == 501
        assert event.untemplated_seq == 'TA'

    def test_call_by_split_read_shift_overlap(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 18), (CIGAR.EQ, 22)],
                query_sequence='TCGGCTCCCGTACTTGTGTATAAGGGGCTTCTGATGTTAT',
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='ATAACATCAGAAGCCCCTTATACACAAGTACGGGAGCCGA',
                is_reverse=True,
            )
        )

        event = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)[0]

        assert event.break1.start == 101
        assert event.break1.end == 101
        assert event.break2.start == 503
        assert event.break2.end == 503
        assert event.untemplated_seq == ''

    def test_both_by_flanking_pairs(self, inversion_evidence):
        inversion_evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t1', reference_id=0, reference_start=150, reference_end=150),
                MockRead(query_name='t1', reference_id=0, reference_start=500, reference_end=520),
            )
        )
        inversion_evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(query_name='t2', reference_id=0, reference_start=120, reference_end=140),
                MockRead(query_name='t2', reference_id=0, reference_start=520, reference_end=520),
            )
        )
        bpp = call._call_by_flanking_pairs(inversion_evidence, SVTYPE.INV)
        # 120-149  ..... 500-519
        # max frag = 150 - 80 = 70
        assert bpp.break1.start == 42
        assert bpp.break1.end == 120
        assert bpp.break2.start == 412  # 70 - 21 = 49
        assert bpp.break2.end == 500

    def test_by_split_reads_multiple_calls(self, inversion_evidence):
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t1',
                reference_start=100,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='A' * 40,
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t1',
                reference_start=500,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='T' * 40,
            )
        )
        inversion_evidence.split_reads[0].add(
            MockRead(
                query_name='t2',
                reference_start=110,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='T' * 40,
            )
        )
        inversion_evidence.split_reads[1].add(
            MockRead(
                query_name='t2',
                reference_start=520,
                cigar=[(CIGAR.S, 20), (CIGAR.EQ, 20)],
                query_sequence='A' * 40,
            )
        )

        evs = call._call_by_split_reads(inversion_evidence, SVTYPE.INV)
        assert len(evs) == 2

    def test_call_by_split_reads_consume_flanking(self):
        evidence = GenomeEvidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            BAM_CACHE,
            REFERENCE_GENOME,
            opposing_strands=True,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            config={
                'validate.stdev_count_abnormal': 3,
                'validate.min_flanking_pairs_resolution': 1,
                'validate.min_splits_reads_resolution': 1,
                'validate.min_linking_split_reads': 1,
            },
        )
        evidence.split_reads[0].add(
            MockRead(
                query_name='test1',
                cigar=[(CIGAR.S, 110), (CIGAR.EQ, 40)],
                reference_start=1114,
                reference_end=1150,
            )
        )
        evidence.split_reads[0].add(
            MockRead(
                query_name='test2',
                cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)],
                reference_start=1108,
                reference_end=1115,
            )
        )
        evidence.split_reads[0].add(
            MockRead(
                query_name='test3',
                cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=1114,
                reference_end=1154,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)],
            )
        )
        evidence.split_reads[1].add(
            MockRead(
                query_name='test4', cigar=[(CIGAR.EQ, 30), (CIGAR.S, 120)], reference_start=2187
            )
        )
        evidence.split_reads[1].add(
            MockRead(
                query_name='test5', cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)], reference_start=2187
            )
        )
        evidence.split_reads[1].add(
            MockRead(
                query_name='test1',
                cigar=[(CIGAR.S, 30), (CIGAR.EQ, 120)],
                reference_start=2187,
                reference_end=2307,
                tags=[(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1)],
            )
        )

        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    query_name='t1',
                    reference_id=3,
                    reference_start=1200,
                    reference_end=1250,
                    is_reverse=True,
                ),
                MockRead(reference_id=3, reference_start=2250, reference_end=2300, is_reverse=True),
            )
        )

        events = call._call_by_split_reads(evidence, event_type=SVTYPE.INV)
        for ev in events:
            print(ev, ev.event_type, ev.call_method)
        assert len(events) == 1
        event = events[0]
        assert len(event.flanking_pairs) == 1
        assert len(event.break1_split_reads) == 2
        assert len(event.break2_split_reads) == 2
        b1 = set([read.query_name for read in event.break1_split_reads])
        b2 = set([read.query_name for read in event.break2_split_reads])
        assert len(b1 & b2) == 1


@pytest.fixture
def left_right_ev():
    return GenomeEvidence(
        Breakpoint('fake', 100, orient=ORIENT.LEFT),
        Breakpoint('fake', 200, orient=ORIENT.RIGHT),
        BamCache(MockBamFileHandle()),
        None,
        opposing_strands=False,
        read_length=25,
        stdev_fragment_size=25,
        median_fragment_size=100,
        config={
            'validate.stdev_count_abnormal': 2,
            'validate.min_flanking_pairs_resolution': 1,
            'validate.min_call_complexity': 0,
        },
    )


class TestCallByFlankingReadsGenome:
    def test_call_coverage_too_large(self):
        with pytest.raises(AssertionError):
            call._call_interval_by_flanking_coverage(
                Interval(1901459, 1902200),
                ORIENT.RIGHT,
                725 + 150,
                150,
                Evidence.distance,
                Evidence.traverse,
            )

    def test_intrachromosomal_lr(self, left_right_ev):
        # --LLL-100------------500-RRR-------
        # max fragment size: 100 + 2 * 25 = 150
        # max distance = 150 - read_length = 125
        # coverage ranges: 20->80 (61)   600->675 (76)
        assert left_right_ev.max_expected_fragment_size == 150
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=19,
                    reference_end=60,
                    next_reference_start=599,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=599,
                    reference_end=650,
                    next_reference_start=19,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=39,
                    reference_end=80,
                    next_reference_start=649,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=649,
                    reference_end=675,
                    next_reference_start=39,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        # add a pair that will be ignored
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=39,
                    reference_end=50,
                    next_reference_start=91,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=91,
                    reference_end=110,
                    next_reference_start=39,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        bpp = call._call_by_flanking_pairs(left_right_ev, SVTYPE.DEL)
        print(bpp, bpp.flanking_pairs)
        assert bpp.break1.start == 80
        assert bpp.break1.end == 80 + 125 - 45
        assert bpp.break2.start == 600 - 125 + 75
        assert bpp.break2.end == 600

    def test_intrachromosomal_lr_coverage_overlaps_range(self, left_right_ev):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=21,
                    reference_end=60,
                    next_reference_start=80,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=80,
                    reference_end=120,
                    next_reference_start=21,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=41,
                    reference_end=80,
                    next_reference_start=110,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=110,
                    reference_end=140,
                    next_reference_start=41,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        # pair to skip
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=39,
                    reference_end=80,
                    next_reference_start=649,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=649,
                    reference_end=675,
                    next_reference_start=39,
                    query_alignment_length=25,
                    is_reverse=True,
                ),
            )
        )
        break1, break2 = call._call_by_flanking_pairs(left_right_ev, SVTYPE.INS)
        assert break1.start == 80
        assert break1.end == 80  # 119
        assert break2.start == 81
        assert break2.end == 81

    def test_intrachromosomal_flanking_coverage_overlap_error(self, left_right_ev):
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=19,
                    reference_end=60,
                    next_reference_start=599,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=599,
                    reference_end=650,
                    next_reference_start=19,
                    query_alignment_length=25,
                ),
            )
        )
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=620,
                    reference_end=80,
                    next_reference_start=780,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=780,
                    reference_end=820,
                    next_reference_start=620,
                    query_alignment_length=25,
                ),
            )
        )
        with pytest.raises(AssertionError):
            call._call_by_flanking_pairs(left_right_ev, SVTYPE.DEL)

    def test_coverage_larger_than_max_expected_variance_error(self, left_right_ev):
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=19,
                    reference_end=60,
                    next_reference_start=599,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=599,
                    reference_end=650,
                    next_reference_start=19,
                    query_alignment_length=25,
                ),
            )
        )
        left_right_ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=301,
                    reference_end=350,
                    next_reference_start=780,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=780,
                    reference_end=820,
                    next_reference_start=301,
                    query_alignment_length=25,
                ),
            )
        )
        with pytest.raises(AssertionError):
            call._call_by_flanking_pairs(left_right_ev, SVTYPE.DEL)

    def test_close_to_zero(self, left_right_ev):
        # this test is for ensuring that if a theoretical window calculated for the
        # first breakpoint overlaps the actual coverage for the second breakpoint (or the reverse)
        # that we adjust the theoretical window accordingly
        ev = GenomeEvidence(
            Breakpoint('fake', 100, orient=ORIENT.RIGHT),
            Breakpoint('fake', 500, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle()),
            None,
            opposing_strands=True,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=180,
            config={
                'validate.stdev_count_abnormal': 2,
                'validate.min_flanking_pairs_resolution': 1,
            },
        )
        ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=19,
                    reference_end=60,
                    next_reference_start=149,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=149,
                    reference_end=150,
                    next_reference_start=19,
                    query_alignment_length=25,
                ),
            )
        )
        ev.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=39,
                    reference_end=80,
                    next_reference_start=199,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=199,
                    reference_end=200,
                    next_reference_start=39,
                    query_alignment_length=25,
                ),
            )
        )
        break1, break2 = call._call_by_flanking_pairs(ev, SVTYPE.INV)

        assert break1.start == 1
        assert break1.end == 20
        assert break2.start == 65
        assert break2.end == 150

    def test_call_with_overlapping_coverage_intervals(self, left_right_ev):
        evidence = GenomeEvidence(
            Breakpoint('1', 76185710, 76186159, orient=ORIENT.RIGHT),
            Breakpoint('1', 76186430, 76186879, orient=ORIENT.LEFT),
            BamCache(MockBamFileHandle()),
            None,
            opposing_strands=False,
            read_length=150,
            stdev_fragment_size=98,
            median_fragment_size=433,
            config={'validate.min_flanking_pairs_resolution': 1},
        )
        evidence.flanking_pairs.add(
            mock_read_pair(
                MockRead(
                    reference_start=76186159,
                    reference_end=76186309,
                    next_reference_start=76186000,
                    query_alignment_length=25,
                ),
                MockRead(
                    reference_start=76186000,
                    reference_end=76186150,
                    next_reference_start=76186159,
                    query_alignment_length=25,
                ),
            )
        )
        with pytest.raises(AssertionError):
            call._call_by_flanking_pairs(evidence, SVTYPE.DUP)


class TestCallByFlankingReadsTranscriptome:
    def build_transcriptome_evidence(self, b1, b2, opposing_strands=False):
        return TranscriptomeEvidence(
            {},  # fake the annotations
            b1,
            b2,
            BamCache(MockBamFileHandle(), stranded=True),
            None,  # bam_cache and reference_genome
            opposing_strands=opposing_strands,
            stranded=True,
            read_length=50,
            stdev_fragment_size=100,
            median_fragment_size=100,
            config={
                'validate.stdev_count_abnormal': 3,
                'validate.min_splits_reads_resolution': 1,
                'validate.min_flanking_pairs_resolution': 1,
                'validate.strand_determining_read': 2,
                'validate.min_call_complexity': 0,
            },
        )

    @todo
    def test_call_translocation(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        pass

    @todo
    def test_call_inversion(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        pass

    @todo
    def test_call_inversion_overlapping_breakpoint_calls(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        pass

    def test_call_deletion(self):
        # transcriptome test will use exonic coordinates for the associated transcripts
        pre_transcript = PreTranscript(
            [(1001, 1100), (1501, 1700), (2001, 2100), (2201, 2300)], strand='+'
        )
        for patt in pre_transcript.generate_splicing_patterns():
            pre_transcript.transcripts.append(Transcript(pre_transcript, patt))
        evidence = self.build_transcriptome_evidence(
            Breakpoint('1', 1051, 1051, 'L', '+'), Breakpoint('1', 1551, 1551, 'R', '+')
        )
        # now add the flanking pairs
        pair = mock_read_pair(
            MockRead(
                'name', '1', 951, 1051, is_reverse=False, query_alignment_length=50, is_read1=False
            ),
            MockRead('name', '1', 2299, 2399, is_reverse=True, query_alignment_length=50),
        )
        print(read_pair_type(pair[0]))
        # following help in debugging the mockup
        assert not pair[0].is_reverse
        assert not pair[0].is_read1
        assert pair[0].is_read2
        assert pair[1].is_reverse
        assert pair[1].is_read1
        assert not pair[1].is_read2
        assert sequenced_strand(pair[0], 2) == STRAND.POS
        assert evidence.decide_sequenced_strand([pair[0]]) == STRAND.POS
        assert sequenced_strand(pair[1], 2) == STRAND.POS
        assert evidence.decide_sequenced_strand([pair[1]]) == STRAND.POS
        print(evidence.max_expected_fragment_size, evidence.read_length)
        evidence.flanking_pairs.add(pair)
        breakpoint1, breakpoint2 = call._call_by_flanking_pairs(evidence, SVTYPE.DEL)
        print(breakpoint1, breakpoint2)
        assert breakpoint1 == Breakpoint('1', 1051, 1351, 'L', '+')
        assert breakpoint2 == Breakpoint('1', 2000, 2300, 'R', '+')

        # now add the transcript and call again
        evidence.overlapping_transcripts.add(pre_transcript)
        breakpoint1, breakpoint2 = call._call_by_flanking_pairs(evidence, SVTYPE.DEL)
        print(breakpoint1, breakpoint2)
        assert breakpoint1 == Breakpoint('1', 1051, 2051, 'L', '+')
        assert breakpoint2 == Breakpoint('1', 1600, 2300, 'R', '+')


class TestCallBySpanningReads:
    def test_deletion(self):
        # ATCGATCTAGATCTAGGATAGTTCTAGCAGTCATAGCTAT
        ev = GenomeEvidence(
            Breakpoint('fake', 60, orient=ORIENT.LEFT),
            Breakpoint('fake', 70, orient=ORIENT.RIGHT),
            BamCache(MockBamFileHandle()),
            None,
            opposing_strands=False,
            read_length=40,
            stdev_fragment_size=25,
            median_fragment_size=180,
            config={
                'validate.min_flanking_pairs_resolution': 1,
                'validate.min_spanning_reads_resolution': 1,
            },
        )
        print(ev.outer_window1, ev.outer_window2)
        spanning_reads = [
            SamRead(
                query_name='name',
                reference_name='fake',
                reference_start=50,
                cigar=[(CIGAR.EQ, 15), (CIGAR.D, 5), (CIGAR.I, 2), (CIGAR.EQ, 10)],
                query_sequence='ATCGATCTAGATCTA' 'GG' 'ATAGTTCTAG',
            ),
            SamRead(
                query_name='name2',
                reference_name='fake',
                reference_start=50,
                cigar=[(CIGAR.EQ, 15), (CIGAR.I, 2), (CIGAR.D, 5), (CIGAR.EQ, 10)],
                query_sequence='ATCGATCTAGATCTA' 'GG' 'ATAGTTCTAG',
            ),
        ]
        ev.spanning_reads = set(spanning_reads)
        calls = call._call_by_spanning_reads(ev, set())
        assert len(calls) == 1
        assert len(calls[0].support()) == 2

    def test_insertion(self):
        pass

    def test_indel(self):
        pass

    def test_inversion(self):
        pass

    def test_duplication(self):
        pass


class TestCharacterizeRepeatRegion:
    def test_bad_deletion_call(self):
        reference_genome = {
            '19': mock.Mock(
                seq=MockLongString(
                    'AAATCTTTTTTCCATTATGGCTATACAAAGTGAATACATTTCCACAAGCAAATATGATAGATTAATTGGTGCATTGTATATATTTCTCAAACCATCAGCTCCTCTT'
                    'TTTTTCAAAGTCTAGAATTTGTAATGGTGGATATCTCTGTTCTGTATTCTGTTGTCTAGATATCCAAGTTTAATGCAAAATTTTATGACATGGAACTTGACACTTT'
                    'CTAGAAATGTTCACATATGGTTGTTTATTAAATTATCTCTCATGGAAATATTTAAATGACATGTTTATTGTCTGAAAAGGACAGATATTTAAGCTTTTTTTTTTTT'
                    'TTTTTCTTTTTTTTTGAGAAAGAGTCTCGTTCTTTTGCCCAGGCTGGAGTGCAGTGGTACAATCTTGGCTCACTACAACCTTCACCTCGCAGGTTCAAGCGATTCT'
                    'CCTGCCTCAGCCTCCCTAGTAGCTGGGATTACAGGTACACACCACCAGGCCCATCTAATTTTTCTATATTTAGTAGAGACAGGGTTTCACCATGTTGGCCAGGCTG'
                    'TTCTCAAACTCCTGACCTCAGCCAATCCGCCCGCCTCAACCTCTTAAAGTGATGGGATTACAGGTGTGAGCCATTGTGCTTGGCCCCCTTTAACTATTTTATGTGA'
                    'CTCTTCT',
                    offset=23951075,
                )
            )
        }
        bpp = BreakpointPair(
            Breakpoint('19', 23951407, orient=ORIENT.LEFT),
            Breakpoint('19', 23951408, orient=ORIENT.RIGHT),
            opposing_strands=False,
            untemplated_seq='',
            event_type=SVTYPE.DEL,
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (0, '')

    def test_homopolymer_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 120, orient=ORIENT.LEFT),
            Breakpoint('1', 121, orient=ORIENT.RIGHT),
            untemplated_seq='T',
            opposing_strands=False,
            event_type=SVTYPE.INS,
        )
        reference_genome = {
            '1': mock.Mock(seq=MockLongString('TCGATTCAGGATCAGATTTTGAACAAGTACATACG', offset=100))
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (4, 'T')

    def test_homopolymer_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 120, orient=ORIENT.LEFT),
            Breakpoint('1', 122, orient=ORIENT.RIGHT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DEL,
        )
        reference_genome = {
            '1': mock.Mock(seq=MockLongString('TCGATTCAGGATCAGATTTTTGAACAAGTACATACG', offset=100))
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (4, 'T')

    def test_homopolymer_duplication(self):
        bpp = BreakpointPair(
            Breakpoint('1', 121, orient=ORIENT.RIGHT),
            Breakpoint('1', 121, orient=ORIENT.LEFT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DUP,
        )
        reference_genome = {
            '1': mock.Mock(seq=MockLongString('TCGATTCAGGATCAGATTTTTGAACAAGTACATACG', offset=100))
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (4, 'T')

    def test_repeat_duplication(self):
        bpp = BreakpointPair(
            Breakpoint('1', 123, orient=ORIENT.RIGHT),
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DUP,
        )
        reference_genome = {
            '1': mock.Mock(
                seq=MockLongString('TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100)
            )
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (2, 'TAG')

    def test_repeat_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TAG',
            opposing_strands=False,
            event_type=SVTYPE.INS,
        )
        reference_genome = {
            '1': mock.Mock(
                seq=MockLongString('TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100)
            )
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (3, 'TAG')

    def test_repeat_deletion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 129, orient=ORIENT.RIGHT),
            untemplated_seq='',
            opposing_strands=False,
            event_type=SVTYPE.DEL,
        )
        reference_genome = {
            '1': mock.Mock(
                seq=MockLongString('TCGATTCAGGATCAGATAGTAGTAGTAGGAACAAGTACATACG', offset=100)
            )
        }
        print(
            'upto and including the second breakpoint',
            reference_genome['1'].seq[bpp.break2.start - 10 : bpp.break2.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (3, 'TAG')

    def test_norepeat_insertion(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.LEFT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TTG',
            opposing_strands=False,
            event_type=SVTYPE.INS,
        )
        reference_genome = {
            '1': mock.Mock(
                seq=MockLongString('TCGATTCAGGATCAGATAGTAGTAGGAACAAGTACATACG', offset=100)
            )
        }
        print(
            'upto and including the first breakpoint',
            reference_genome['1'].seq[bpp.break1.start - 10 : bpp.break1.start],
        )
        assert call.EventCall.characterize_repeat_region(bpp, reference_genome) == (0, 'TTG')

    def test_invalid_event_type(self):
        bpp = BreakpointPair(
            Breakpoint('1', 125, orient=ORIENT.RIGHT),
            Breakpoint('1', 126, orient=ORIENT.RIGHT),
            untemplated_seq='TTG',
            event_type=SVTYPE.INV,
        )
        with pytest.raises(ValueError):
            call.EventCall.characterize_repeat_region(bpp, None)
