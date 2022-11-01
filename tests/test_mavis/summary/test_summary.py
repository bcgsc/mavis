import pytest
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import CALL_METHOD, COLUMNS, PROTOCOL, STRAND, SVTYPE
from mavis.summary.summary import filter_by_annotations, annotate_dgv
from mavis.annotate.file_io import load_masking_regions

from ...util import todo, get_data


@pytest.fixture
def genomic_event1():
    return BreakpointPair(
        Breakpoint('1', 1),
        Breakpoint('1', 10),
        opposing_strands=True,
        **{
            COLUMNS.event_type: SVTYPE.DEL,
            COLUMNS.call_method: CALL_METHOD.CONTIG,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.protocol: PROTOCOL.GENOME,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.fusion_cdna_coding_start: None,
        }
    )


@pytest.fixture
def genomic_event2():
    return BreakpointPair(
        Breakpoint('1', 1),
        Breakpoint('1', 100),
        opposing_strands=True,
        **{
            COLUMNS.event_type: SVTYPE.DEL,
            COLUMNS.call_method: CALL_METHOD.CONTIG,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.protocol: PROTOCOL.GENOME,
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
        }
    )


@pytest.fixture
def genomic_event3():
    return BreakpointPair(
        Breakpoint('1', 10001),
        Breakpoint('1', 22118),
        opposing_strands=True,
        **{
            COLUMNS.event_type: SVTYPE.DEL,
            COLUMNS.call_method: CALL_METHOD.CONTIG,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.protocol: PROTOCOL.GENOME,
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.tracking_id: "genomic_event3",
        }
    )


@pytest.fixture
def best_transcripts():
    return {'ABCA': True, 'ABCD': True}


@pytest.fixture
def dgv_event2():
    return load_masking_regions(get_data("mock_dgv_annotation.txt"))


class TestFilterByAnnotations:
    def test_filter_by_annotations_two_best_transcripts(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = 'ABC'
        genomic_event1.data[COLUMNS.gene2] = 'ABC'
        genomic_event1.data[COLUMNS.transcript1] = 'ABCA'
        genomic_event1.data[COLUMNS.transcript2] = 'ABCA'
        genomic_event2.data[COLUMNS.gene1] = 'ABC'
        genomic_event2.data[COLUMNS.gene2] = 'ABC'
        genomic_event2.data[COLUMNS.transcript1] = 'ABCD'
        genomic_event2.data[COLUMNS.transcript2] = 'ABCD'
        result, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        bpp = result[0]
        print(bpp.data)
        assert bpp == genomic_event1
        assert bpp.data[COLUMNS.transcript1] == 'ABCA'

    def test_filter_by_annotations_two_transcripts(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = 'XYZ'
        genomic_event1.data[COLUMNS.gene2] = 'XYS'
        genomic_event1.data[COLUMNS.transcript1] = 'XYZB'
        genomic_event1.data[COLUMNS.transcript2] = 'XYSZ'
        genomic_event2.data[COLUMNS.gene1] = 'XYZ'
        genomic_event2.data[COLUMNS.gene2] = 'XYS'
        genomic_event2.data[COLUMNS.transcript1] = 'XYZA'
        genomic_event2.data[COLUMNS.transcript2] = 'XYSB'
        bpps, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        print(bpps)
        bpp = bpps[0]
        print(bpp, bpp.data)
        assert bpp == genomic_event2
        assert bpp.data[COLUMNS.transcript1] == 'XYZA'

    def test_filter_by_annotations_two_fusion_cdna(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = 'XYZ'
        genomic_event1.data[COLUMNS.gene2] = 'XYS'
        genomic_event1.data[COLUMNS.transcript1] = 'XYZB'
        genomic_event1.data[COLUMNS.transcript2] = 'XYSZ'
        genomic_event2.data[COLUMNS.gene1] = 'XYZ'
        genomic_event2.data[COLUMNS.gene2] = 'XYS'
        genomic_event2.data[COLUMNS.transcript1] = 'XYZB'
        genomic_event2.data[COLUMNS.transcript2] = 'XYSZ'
        genomic_event1.data[COLUMNS.fusion_cdna_coding_start] = 1
        genomic_event1.data[COLUMNS.fusion_cdna_coding_end] = 20
        genomic_event2.data[COLUMNS.fusion_cdna_coding_start] = 1
        genomic_event2.data[COLUMNS.fusion_cdna_coding_end] = 40
        result, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        bpp = result[0]
        assert bpp == genomic_event2

    def test_filter_by_annotations_one_transcript(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = None
        genomic_event1.data[COLUMNS.gene2] = 'XYS'
        genomic_event1.data[COLUMNS.transcript1] = None
        genomic_event1.data[COLUMNS.transcript2] = 'XYSZ'
        genomic_event2.data[COLUMNS.gene1] = 'XYZ'
        genomic_event2.data[COLUMNS.gene2] = 'XYS'
        genomic_event2.data[COLUMNS.transcript1] = 'XYZA'
        genomic_event2.data[COLUMNS.transcript2] = 'XYSB'
        result, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        bpp = result[0]
        assert bpp == genomic_event2

    def test_filter_by_annotations_one_best_transcripts(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = 'XYZ'
        genomic_event1.data[COLUMNS.gene2] = 'ABC'
        genomic_event1.data[COLUMNS.transcript1] = 'XYZB'
        genomic_event1.data[COLUMNS.transcript2] = 'ABCA'
        genomic_event2.data[COLUMNS.gene1] = 'XYZ'
        genomic_event2.data[COLUMNS.gene2] = 'ABC'
        genomic_event2.data[COLUMNS.transcript1] = 'XYZA'
        genomic_event2.data[COLUMNS.transcript2] = 'ABCB'
        result, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        bpp = result[0]
        assert bpp == genomic_event1
        assert bpp.data[COLUMNS.transcript1] == 'XYZB'

    def test_filter_by_annotations_no_transcripts(
        self, genomic_event1, genomic_event2, best_transcripts
    ):
        genomic_event1.data[COLUMNS.gene1] = None
        genomic_event1.data[COLUMNS.gene2] = None
        genomic_event1.data[COLUMNS.transcript1] = None
        genomic_event1.data[COLUMNS.transcript2] = None
        genomic_event2.data[COLUMNS.gene1] = None
        genomic_event2.data[COLUMNS.gene2] = None
        genomic_event2.data[COLUMNS.transcript1] = None
        genomic_event2.data[COLUMNS.transcript2] = None
        genomic_event1.break1.strand = STRAND.POS
        result, removed = filter_by_annotations([genomic_event1, genomic_event2], best_transcripts)
        bpp = result[0]
        assert bpp.data[COLUMNS.transcript1] is None

    @todo
    def test_combine_events(self):
        pass

    @todo
    def test_filtering_events_contigs(self):
        pass

    @todo
    def test_filtering_events_none(self):
        pass

    @todo
    def test_filtering_events_flanking(self):
        pass

    @todo
    def test_filtering_events_spanning(self):
        pass

    @todo
    def test_filtering_events_split(self):
        pass

    @todo
    def test_get_pairing_state(self):
        pass


class TestFilterByCallMethod:
    def test_annotate_dgv_distance_bed(self, genomic_event3, dgv_event2):
        bpps = [genomic_event3]
        annotate_dgv(bpps, dgv_event2, 103)
        print(bpps[0].data)
        assert len(bpps[0].data['dgv']) == 3
        assert (
            bpps[0].data
            == "['dgv1n82(1:10001-22118)', 'rgv2n98(1:10001-22120)', 'rgv2n99(1:10001-22221)']"
        )
