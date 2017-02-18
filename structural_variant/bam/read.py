from ..constants import ORIENT, CIGAR, DNA_ALPHABET, STRAND, READ_PAIR_TYPE, SVTYPE
from . import cigar as cigar_tools
import pysam


def breakpoint_pos(read, orient=ORIENT.NS):
    """
    assumes the breakpoint is the position following softclipping on the side with more
    softclipping (unless and orientation has been specified)

    Args:
        read (psyam.AlignedSegment): the read object
        orient (ORIENT): the orientation

    Returns:
        int: the position of the breakpoint in the input read
    """
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping')

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
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint')
        return read.reference_start
    else:
        if end_typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint')
        return read.reference_end - 1


def nsb_align(ref, seq, weight_of_score=0.5, min_overlap_percent=100, min_match=0):
    """
    given some reference string and a smaller sequence string computes the best non-space-breaking alignment
    i.e. an alignment that does not allow for indels (straight-match)
    """
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence')
    if min_match < 0 or min_match > 1:
        raise AttributeError('min_match must be between 0 and 1')

    if min_overlap_percent <= 0 or min_overlap_percent > 100:
        raise AttributeError('percent must be greater than 0 and up to 100', min_overlap_percent)

    min_overlap = int(round(min_overlap_percent * len(seq) / 100, 0))

    # store to improve speed and space (don't need to store all alignments)
    best_score = 0
    results = []

    for ref_start in range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap):
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

        score = cigar_tools.score(cigar) * weight_of_score + \
            cigar_tools.longest_exact_match(cigar) * (1 - weight_of_score)
        score = cigar_tools.score(cigar)
        a = pysam.AlignedSegment()
        a.query_sequence = str(seq)
        a.reference_start = ref_start + qstart
        a.cigar = cigar

        if score >= best_score:
            best_score = score
            results.append((a, score))

    filtered = [x for x, y in results if y == best_score]

    return filtered


def read_pair_strand(read, strand_determining_read=2):
    if not read.is_paired:
        return STRAND.NEG if read.is_reverse else STRAND.POS
    elif strand_determining_read == 1:
        if read.is_read1:
            return STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            return STRAND.NEG if read.mate_is_reverse else STRAND.POS
    elif strand_determining_read == 2:
        if read.is_read2:
            return STRAND.NEG if read.is_reverse else STRAND.POS
        else:
            return STRAND.NEG if read.mate_is_reverse else STRAND.POS
    else:
        raise ValueError('unexpected value. Expected 1 or 2, found:', strand_determining_read)


def read_pair_type(read):
    # check if the read pair is in the expected orientation
    """
    assumptions based on illumina pairs: only 4 possible combinations

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
    checks if the orientation is compatible with the type and also
    if the insert size makes sense
    """
    if event_type == SVTYPE.DEL:
        if read_pair_type(read) != READ_PAIR_TYPE.LR:
            return False
    elif event_type == SVTYPE.INS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR:
            return False
    elif event_type == SVTYPE.TRANS:
        if read_pair_type(read) != READ_PAIR_TYPE.LR and \
                read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    elif event_type == SVTYPE.ITRANS:
        if read_pair_type(read) != READ_PAIR_TYPE.LL and \
                read_pair_type(read) != READ_PAIR_TYPE.RR:
            return False
    elif event_type == SVTYPE.INV:
        if read_pair_type(read) != READ_PAIR_TYPE.LL and \
                read_pair_type(read) != READ_PAIR_TYPE.RR:
            return False
    elif event_type == SVTYPE.DUP:
        if read_pair_type(read) != READ_PAIR_TYPE.RL:
            return False
    else:
        raise ValueError('unexpected event type', event_type)
    return True
