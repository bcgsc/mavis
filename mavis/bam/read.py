from ..constants import ORIENT, CIGAR, DNA_ALPHABET, STRAND, READ_PAIR_TYPE, SVTYPE
from . import cigar as cigar_tools
from .cigar import EVENT_STATES, REFERENCE_ALIGNED_STATES, ALIGNED_STATES, QUERY_ALIGNED_STATES
import pysam
import subprocess
import re
from copy import copy


def get_samtools_version():
    proc = subprocess.getoutput(['samtools'])
    for line in proc.split('\n'):
        m = re.search('Version: (?P<major>\d+)(\.(?P<mid>\d+)(\.(?P<minor>\d+))?)?', line)
        if m:
            major = int(m.group('major'))
            mid = int(m.group('mid')) if m.group('mid') else 0
            minor = int(m.group('minor')) if m.group('minor') else 0
            return major, mid, minor
    raise ValueError('unable to parse samtools version number')


def samtools_v0_sort(input_bam, output_bam):
    prefix = re.sub('\.bam$', '', output_bam)
    return 'samtools sort {} {}'.format(input_bam, prefix)


def samtools_v1_sort(input_bam, output_bam):
    return 'samtools sort {} -o {}'.format(input_bam, output_bam)


def breakpoint_pos(read, orient=ORIENT.NS):
    """
    assumes the breakpoint is the position following softclipping on the side with more
    softclipping (unless and orientation has been specified)

    Args:
        read (:class:`~pysam.AlignedSegment`): the read object
        orient (ORIENT): the orientation

    Returns:
        int: the position of the breakpoint in the input read
    """
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping', read.cigar)

    if orient == ORIENT.NS:
        if (typ == CIGAR.S and end_typ == CIGAR.S and freq > end_freq) \
                or typ == CIGAR.S and end_typ != CIGAR.S:
            orient = ORIENT.RIGHT
            # soft clipped to the left
        else:
            # soft clipped to the right
            orient = ORIENT.LEFT

    if orient == ORIENT.RIGHT:
        if typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint', repr(orient), read.cigar, read.get_tags())
        return read.reference_start
    else:
        if end_typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint', orient, read.cigar, read.get_tags())
        return read.reference_end - 1


def calculate_alignment_score(read):
    """
    calculates a score for comparing alignments

    Args:
        read (pysam.AlignedSegment): the input read

    Returns:
        float: the score
    """
    score = 0
    qlen = read.reference_end - read.reference_start
    max_score = 2 * qlen - 1
    for c, v in read.cigar:
        if c == CIGAR.M:
            raise ValueError('cannot calculate the alignment score if mismatch v match has not been specified')
        elif c == CIGAR.EQ:
            score += 2 * v - 1
    return score / max_score


def nsb_align(
        ref, seq, 
        weight_of_score=0.5, 
        min_overlap_percent=1, 
        min_match=0, 
        min_consecutive_match=1, 
        scoring_function=calculate_alignment_score):
    """
    given some reference string and a smaller sequence string computes the best non-space-breaking alignment
    i.e. an alignment that does not allow for indels (straight-match). Positions in the aligned segments are
    given relative to the length of the reference sequence (1-based)

    Args:
        ref (str): the reference sequence
        seq (str): the sequence being aligned
        weight_of_score (float): when scoring alignments this determines the amount
            of weight to place on the cigar match. Should be a number between 0 and 1
        min_overlap_percent (float): the minimum amount of overlap of the input sequence to the reference
            should be a number between 0 and 1
        min_match (float): the minimum number of matches compared to total
        scoring_function (callable): any function that will take a read as input and return a float
          used in comparing alignments to choose the best alignment

    Returns:
        :class:`list` of :class:`~pysam.AlignedSegment`: list of aligned segments

    Note:
        using a higher min_match may improve performance as low quality alignments are rejected more quickly. However
        this may also result in no match being returned when there is no high quality match to be found.
    """
    ref = str(ref)
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence')
    if min_match < 0 or min_match > 1:
        raise AttributeError('min_match must be between 0 and 1')

    if min_overlap_percent <= 0 or min_overlap_percent > 1:
        raise AttributeError('percent must be greater than 0 and up to 1', min_overlap_percent)

    min_overlap = int(round(min_overlap_percent * len(seq), 0))

    # store to improve speed and space (don't need to store all alignments)
    best_score = (0, 0)
    results = []

    putative_start_positions = range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap)
    if min_consecutive_match > 1:
        putative_start_positions = set()
        kmers_checked = {}
        for i in range(0, len(seq) - min_consecutive_match):
            current_kmer = seq[i:i + min_consecutive_match]
            if current_kmer in kmers_checked:
                putative_start_positions.update([p - i for p in kmers_checked[current_kmer]])
                continue
            rp = [m.start() for m in re.finditer(current_kmer, ref)]
            kmers_checked[current_kmer] = rp
            putative_start_positions.update([p - i for p in rp])
    for ref_start in putative_start_positions:
        score = 0
        cigar = []
        mismatches = 0
        length = len(seq)
        for i in range(0, len(seq)):
            if length == 0:
                break
            r = ref_start + i
            if r < 0 or r >= len(ref):  # outside the length of the reference seq
                cigar.append((CIGAR.S, 1))
                length -= 1
                continue
            if DNA_ALPHABET.match(ref[r], seq[i]):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))
                mismatches += 1
                if mismatches / length > 1 - min_match:
                    break
        if length == 0 or mismatches / length > 1 - min_match:
            continue

        cigar = cigar_tools.join(cigar)
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])

        qstart = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]

        a = pysam.AlignedSegment()
        a.query_sequence = str(seq)
        a.reference_start = ref_start + qstart
        a.cigar = cigar
        qlen = a.reference_end - a.reference_start
        score = (scoring_function(a), qlen)  # this way for equal identity matches we take the longer alignment
        if qlen < min_overlap:
            continue

        if score >= best_score:
            best_score = score
            results.append((a, score))

    filtered = [x for x, y in results if y == best_score]

    return filtered


def sequenced_strand(read, strand_determining_read=2):
    """
    determines the strand that was sequenced

    Args:
        read (:class:`~pysam.AlignedSegment`): the read being used to determine the strand
        strand_determining_read (int): which read in the read pair is the same as the sequenced strand

    Returns:
        STRAND: the strand that was sequenced

    Raises:
        ValueError: if strand_determining_read is not 1 or 2

    Warning:
        if the input pair is unstranded the information will not be representative of the
        strand sequenced since the assumed convention is not followed
    """
    if read.is_unmapped or not read.is_paired:
        raise ValueError('cannot determine strand if the read is unmapped or unpaired')
    strand = None
    if strand_determining_read == 1:
        if read.is_read1:
            strand = STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            strand = STRAND.NEG if not read.is_reverse else STRAND.POS
    elif strand_determining_read == 2:
        if read.is_read2:
            strand = STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            strand = STRAND.NEG if not read.is_reverse else STRAND.POS
    else:
        raise ValueError('unexpected value. Expected 1 or 2, found:', strand_determining_read)
    return strand


def read_pair_type(read):
    # check if the read pair is in the expected orientation
    """
    assumptions based on illumina pairs: only 4 possible combinations

    Args:
        read (:class:`~pysam.AlignedSegment`): the input read

    Returns:
        READ_PAIR_TYPE: the type of input read pair

    Raises:
        NotImplementedError: for any read that does not fall into the four expected configurations (see below)

    ::

        ++++> <---- is LR same-strand
        ++++> ++++> is LL opposite
        <---- <---- is RR opposite
        <---- ++++> is RL same-strand
    """
    reverse = False
    if read.reference_id > read.next_reference_id or \
            (read.reference_id == read.next_reference_id and read.reference_start > read.next_reference_start):
        reverse = True

    if not read.is_reverse and read.mate_is_reverse:  # LR
        return READ_PAIR_TYPE.RL if reverse else READ_PAIR_TYPE.LR
    elif not read.is_reverse and not read.mate_is_reverse:  # LL opp
        return READ_PAIR_TYPE.LL
    elif read.is_reverse and read.mate_is_reverse:  # RR opp
        return READ_PAIR_TYPE.RR
    elif read.is_reverse and not read.mate_is_reverse:  # RL
        return READ_PAIR_TYPE.LR if reverse else READ_PAIR_TYPE.RL
    else:
        raise NotImplementedError('unexpected orientation for pair')


def orientation_supports_type(read, event_type):
    """
    checks if the orientation is compatible with the type of event

    Args:
        read (:class:`~pysam.AlignedSegment`): a read from the pair
        event_type (SVTYPE): the type of event to check

    Returns:
        bool:
            - ``True`` - the read pair is in the correct orientation for this event type
            - ``False`` - the read is not in the correct orientation
    """
    if event_type == SVTYPE.DEL or event_type == SVTYPE.INS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR:
            return False
    elif event_type == SVTYPE.TRANS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR and \
                read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    elif event_type == SVTYPE.ITRANS or event_type == SVTYPE.INV:
        if read_pair_type(read) != READ_PAIR_TYPE.LL and \
                read_pair_type(read) != READ_PAIR_TYPE.RR:
            return False
    elif event_type == SVTYPE.DUP:
        if read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    else:
        raise ValueError('unexpected event type', event_type)
    return True


def convert_events_to_softclipping(read, orientation, max_event_size, min_anchor_size=None):
    """
    given an alignment, simplifies the alignment by grouping everything past the first anchor and including the
    first event considered too large and unaligning them turning them into softclipping
    
    """
    if min_anchor_size is None:
        min_anchor_size = max_event_size
    
    if orientation == ORIENT.LEFT:
        event_size = 0
        adjusted_cigar = []
        anchor = 0
        for state, count in read.cigar:
            if state == CIGAR.M:
                raise NotImplementedError('match v mismatch must be specified')
            elif anchor < min_anchor_size:
                if state == CIGAR.EQ:
                    anchor += count
            elif state in EVENT_STATES:
                event_size += count
                if event_size > max_event_size:
                    break
            else:
                event_size = 0
            adjusted_cigar.append((state, count))
        if event_size > max_event_size:
            while adjusted_cigar[-1][0] in EVENT_STATES:
                del adjusted_cigar[-1]
            aligned = sum([y for x, y in adjusted_cigar if x in QUERY_ALIGNED_STATES] + [0])
            sc = len(read.query_sequence) - aligned
            adjusted_cigar.append((CIGAR.S, sc))
            read = copy(read)
            read.cigar = adjusted_cigar
    elif orientation == ORIENT.RIGHT:
        # more complicated than left b/c need to also adjust the start position
        event_size = 0
        anchor = 0
        adjusted_cigar = []
        for state, count in read.cigar[::-1]:  # first event from the right
            if state == CIGAR.M:
                raise NotImplementedError('match v mismatch must be specified')
            elif anchor < min_anchor_size:
                if state == CIGAR.EQ:
                    anchor += count
            elif state in EVENT_STATES:
                event_size += count
                if event_size > max_event_size:
                    break
            else:
                event_size = 0
            adjusted_cigar.append((state, count))
        if event_size > max_event_size:
            while adjusted_cigar[-1][0] in EVENT_STATES:
                del adjusted_cigar[-1]
            originally_refaligned = sum([y for x, y in read.cigar if x in REFERENCE_ALIGNED_STATES] + [0])
            refaligned = sum([y for x, y in adjusted_cigar if x in REFERENCE_ALIGNED_STATES] + [0])
            aligned = sum([y for x, y in adjusted_cigar if x in QUERY_ALIGNED_STATES] + [0])
            sc = len(read.query_sequence) - aligned
            adjusted_cigar = [(CIGAR.S, sc)] + adjusted_cigar[::-1]
            read = copy(read)
            read.cigar = adjusted_cigar
            read.reference_start += originally_refaligned - refaligned
    else:
        raise ValueError('orientation must be specified', orientation)
    return read


