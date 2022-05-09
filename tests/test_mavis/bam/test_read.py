from typing import List, Union

import pysam
import pytest
from black import read_cache
from build.lib.mavis.bam.cigar import convert_cigar_to_string
from build.lib.mavis.bam.cigar import join as join_cigar
from mavis.bam.cigar import QUERY_ALIGNED_STATES, REFERENCE_ALIGNED_STATES
from mavis.bam.read import SamRead, simplify_long_read
from mavis.constants import CIGAR


def read_repr(read):
    return f'{read.query_name}:{read.query_alignment_start}-{read.query_alignment_end}|{read.reference_name}:{read.reference_start}-{read.reference_end}'


class DiscontinuousAlignment:
    segments: List[Union[pysam.AlignedSegment, str, None]]  # none indicates unknown untemplated seq

    def __init__(self, segments: List[Union[pysam.AlignedSegment, str, None]]) -> None:
        self.segments = segments  # connected in order

    def seq(self):
        seq = []
        for segment in self.segments:
            if hasattr(segment, 'seq'):
                seq.append(segment.seq)
            else:
                seq.append(segment)
        return ''.join(seq)

    def __repr__(self):
        rep = []
        for segment in self.segments:
            if hasattr(segment, 'query_alignment_start'):
                rep.append(read_repr(segment))
            else:
                rep.append(segment)
        return 'DA[' + ', '.join(rep) + ']'

    def query_coverage(self):
        pass

    def query_name(self) -> Union[str, None]:
        for segment in self.segments:
            if hasattr(segment, 'query_name'):
                return self.query_name
        return None


def simplify_read(read, reference_seq, min_event_size):
    """
    Given some read where we expect low quality indel artifacts to be present (ex long read) remove all
    events below a given threshold to simplify the read for downstream processing
    """
    pass


def is_hardclipped(read):
    return bool(read.cigartuples[0][0] == CIGAR.H or read.cigartuples[-1][0] == CIGAR.H)


def soften_clipping(read_group):
    """
    For a group of alignments on the same query_name (read) that are supplementary alignments of the same sequence. Change the hard clipped reads to be soft clipped versions
    """
    full_seq_reads = []
    hard_clipped = []

    assert len({r.query_name for r in read_group}) == 1

    for read in read_group:
        if is_hardclipped(read):
            hard_clipped.append(read)
        else:
            full_seq_reads.append(read)

    if not hard_clipped:
        return read_group

    if hard_clipped and not full_seq_reads:
        raise ValueError(
            'cannot resolve hard clipping as there are no reads in this group that are not hard-clipped'
        )

    full_seq = {r.seq for r in full_seq_reads}

    if len(full_seq) > 1:
        raise AssertionError(
            'reads that are not hard clipped with the same names have different sequences'
        )

    query_sequence = list(full_seq)[0]

    result = full_seq_reads

    for read in hard_clipped:
        if read.seq not in query_sequence:
            raise AssertionError(
                'hard clipped sequence is not a subset of the unclipped query sequence'
            )
        sc_read = SamRead.copy(read)
        sc_read.seq = query_sequence
        sc_read.cigartuples = [
            (state if state != CIGAR.H else CIGAR.S, count) for state, count in read.cigartuples
        ]
        result.append(sc_read)
    return result


def shift_query_alignment_start(read: pysam.AlignedSegment, query_start_boundary: int):
    """
    adjust the read CIGAR tuples to ensure that everything before the new query_alignment_start is alterted to soft-clipping
    """
    if read.query_alignment_end <= query_start_boundary:
        raise ValueError(
            f'read ({read.query_name}:{read.query_alignment_start}-{read.query_alignment_end}) does not cover the expected query_start_boundary ({query_start_boundary}) region'
        )

    if read.query_alignment_start >= query_start_boundary:
        return read
    new_cigar = []
    current_query_start = read.query_alignment_start
    current_ref_start = read.reference_start

    for cigar_state, size in read.cigartuples:
        if current_query_start >= query_start_boundary:
            new_cigar.append((cigar_state, size))
            continue

        if cigar_state in {CIGAR.H, CIGAR.S}:
            new_cigar.append((cigar_state, size))
        elif cigar_state in {CIGAR.D, CIGAR.N}:
            # shift the reference pos but do not keep the CIGAR tuple
            current_ref_start += size
        elif cigar_state in QUERY_ALIGNED_STATES:
            # need to unalign part of all of this segment
            shift = query_start_boundary - current_query_start
            if shift < size:
                # only need part of the section

                if cigar_state in REFERENCE_ALIGNED_STATES:
                    current_ref_start += shift
                current_query_start += size
                new_cigar.append((CIGAR.S, shift))
                new_cigar.append((cigar_state, size - shift))
            else:
                # use the entire segment
                current_query_start += size
                if cigar_state in REFERENCE_ALIGNED_STATES:
                    current_ref_start += size
                new_cigar.append((CIGAR.S, size))
        else:
            raise NotImplementedError(f'unexpected query cigar state {cigar_state}')

    new_read = SamRead.copy(read)
    new_read.reference_start = current_ref_start
    new_read.cigar = join_cigar(new_cigar)
    return new_read


def merge_alignments(reads: List[pysam.AlignedSegment]) -> DiscontinuousAlignment:
    print([read_repr(r) for r in reads])
    for read in reads:
        print(
            f'{read.query_name}:{read.query_alignment_start}-{read.query_alignment_end}',
            f'{read.reference_name}:{read.reference_start}-{read.reference_end}',
        )

    consumed_seq = -1
    segments = []
    for read in sorted(reads, key=lambda read: read.query_alignment_start):
        # greedy consume the sequence. If overlaps a previously consumed segment of the query sequence, that segment should be skipped
        shifted_read = shift_query_alignment_start(read, consumed_seq)
        consumed_seq = shifted_read.query_alignment_end
        if segments:
            last_segment = segments[-1]
            if last_segment.query_alignment_end < shifted_read.query_alignment_start:
                # add seq
                seq = read.seq[
                    last_segment.query_alignment_end : shifted_read.query_alignment_start
                ]
                segments.append(seq)
        segments.append(shifted_read)
    print([read_repr(r) for r in segments])
    return DiscontinuousAlignment(segments)


def convert_supplementary_to_discontinuous(
    reads: List[pysam.AlignedSegment],
) -> List[DiscontinuousAlignment]:
    reads_by_name = {}
    for read in reads:
        reads_by_name.setdefault(read.query_name, []).append(read)

    result = []
    for read_group in reads_by_name.values():
        if len(read_group) < 2:
            result.append(DiscontinuousAlignment(soften_clipping(read_group)))
            continue
        print('merging read group', len(read_group))
        result.append(merge_alignments(soften_clipping(read_group)))
    return result


@pytest.fixture(scope='module')
def sample_long_reads():
    with pysam.AlignmentFile(
        '/projects/jfan_prj/jfan_prj/Nanopore_Testing/2021_nanopore_sv_testing/scratch/depth_testing/POG/COLO829/minimap2_bam/assembly_testing/april_2022/redbeans_assembly/xd_wtdbg2_sorted.bam',
        'rb',
    ) as fh:
        all_reads = list(fh.fetch())
    return all_reads


def test_simplify_long_read_pair(sample_long_reads):
    assert len(sample_long_reads) == 4
    print([read_repr(r) for r in sample_long_reads])
    alignments = convert_supplementary_to_discontinuous(sample_long_reads)
    assert len(alignments) == 3


@pytest.mark.parametrize(
    'query_start_boundary,read_query_alignment_start', [[19851, 19851], [19000, 19850]]
)
def test_shift_query_alignment_start(
    query_start_boundary, read_query_alignment_start, sample_long_reads
):
    ctg1_reads = soften_clipping([r for r in sample_long_reads if r.query_name == 'ctg1'])
    target_reads = [r for r in ctg1_reads if r.reference_name == 'chr10']
    assert len(target_reads) == 1
    target_read = target_reads[0]

    assert target_read.query_alignment_start == 19850
    target_read = shift_query_alignment_start(target_read, query_start_boundary)
    assert target_read.query_alignment_start == read_query_alignment_start
