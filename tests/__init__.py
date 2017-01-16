from structural_variant.constants import CIGAR
from structural_variant.blat import BlatAlignedSegment
import os


filedir = os.path.join(os.path.dirname(__file__), 'files')
REFERENCE_GENOME_FILE = os.path.join(filedir, 'mock_reference_genome.fa')
REFERENCE_ANNOTATIONS_FILE = os.path.join(filedir, 'mock_reference_annotations.tsv')
BAM_INPUT = os.path.join(filedir, 'mock_reads_for_events.sorted.bam')
BASE_EVENTS = os.path.join(filedir, 'mock_sv_events.svmerge.tsv')
BLAT_INPUT = os.path.join(filedir, 'blat_input.fa')
BLAT_OUTPUT = os.path.join(filedir, 'blat_output.pslx')

class MockRead:
    def __init__(
        self,
        query_name=None,
        reference_id=None,
        reference_start=None,
        reference_end=None,
        cigar=None,
        is_reverse=False,
        mate_is_reverse=True,
        next_reference_start=None,
        next_reference_id=None,
        reference_name=None,
        query_sequence=None,
        template_length=None,
        query_alignment_sequence=None
    ):
        self.query_name = query_name
        self.reference_id = reference_id
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.cigar = cigar
        if self.reference_end is None and cigar and reference_start is not None:
            self.reference_end = reference_start + sum([f for v, f in cigar if v not in [CIGAR.S, CIGAR.I]])
        self.is_reverse = is_reverse
        self.mate_is_reverse = mate_is_reverse
        self.next_reference_start = next_reference_start
        self.next_reference_id = next_reference_id
        self.reference_name = reference_name
        self.query_sequence = query_sequence
        self.query_alignment_sequence = query_alignment_sequence
        if query_alignment_sequence is None and cigar and query_sequence:
            s = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]
            t = len(query_sequence)
            if cigar[-1][0] == CIGAR.S:
                t -= cigar[-1][1]
            self.query_alignment_sequence = query_sequence[s:t]
        if cigar and query_sequence:
            assert(len(query_sequence) == sum([f for v, f in cigar if v not in [CIGAR.H, CIGAR.N, CIGAR.D]]))
        if template_length is None and reference_end and next_reference_start:
            self.template_length = next_reference_start - reference_end
        else:
            self.template_length = template_length

    def query_coverage_interval(self):
        return BlatAlignedSegment.query_coverage_interval(self)


class MockBamFileHandle:
    def __init__(self, chrom_to_tid={}):
        self.chrom_to_tid = chrom_to_tid

    def fetch(self, *pos):
        return []

    def get_tid(self, chrom):
        if chrom in self.chrom_to_tid:
            return self.chrom_to_tid[chrom]
        else:
            return -1

    def get_reference_name(self, input_tid):
        for chrom, tid in self.chrom_to_tid.items():
            if input_tid == tid:
                return chrom
        raise KeyError('invalid id')


class MockSeq:
    def __init__(self, seq=None):
        self.seq = seq


class MockString:
    def __init__(self, char=' '):
        self.char = char

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.char * (index.stop - index.start)
        else:
            return self.char
