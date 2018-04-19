"""
holds methods related to processing cigar tuples. Cigar tuples are generally
an iterable list of tuples where the first element in each tuple is the
CIGAR value (i.e. 1 for an insertion), and the second value is the frequency
"""
import re
from ..constants import CIGAR, DNA_ALPHABET, GAP

EVENT_STATES = {CIGAR.D, CIGAR.I, CIGAR.X}
ALIGNED_STATES = {CIGAR.M, CIGAR.X, CIGAR.EQ}
REFERENCE_ALIGNED_STATES = ALIGNED_STATES | {CIGAR.D, CIGAR.N}
QUERY_ALIGNED_STATES = ALIGNED_STATES | {CIGAR.I, CIGAR.S}
CLIPPING_STATE = {CIGAR.S, CIGAR.H}


def recompute_cigar_mismatch(read, ref):
    """
    for cigar tuples where M is used, recompute to replace with X/= for increased
    utility and specificity

    Args:
        read (pysam.AlignedSegment): the input read
        ref (str): the reference sequence

    Returns:
        :class:`list` of :class:`tuple` of :class:`int` and :class:`int`: the cigar tuple
    """
    result = []
    offset = 0

    ref_pos = read.reference_start
    seq_pos = 0

    for cigar_value, freq in read.cigar:
        if cigar_value in ALIGNED_STATES:
            for offset in range(0, freq):
                if DNA_ALPHABET.match(ref[ref_pos], read.query_sequence[seq_pos]):
                    if len(result) == 0 or result[-1][0] != CIGAR.EQ:
                        result.append((CIGAR.EQ, 1))
                    else:
                        result[-1] = (CIGAR.EQ, result[-1][1] + 1)
                else:
                    if len(result) == 0 or result[-1][0] != CIGAR.X:
                        result.append((CIGAR.X, 1))
                    else:
                        result[-1] = (CIGAR.X, result[-1][1] + 1)
                ref_pos += 1
                seq_pos += 1
            continue
        if cigar_value in QUERY_ALIGNED_STATES:
            seq_pos += freq
        if cigar_value in REFERENCE_ALIGNED_STATES:
            ref_pos += freq
        result.append((cigar_value, freq))
    assert sum([x[1] for x in result]) == sum(x[1] for x in read.cigar)
    return result


def longest_fuzzy_match(cigar, max_fuzzy_interupt=1):
    """
    computes the longest sequence of exact matches allowing for 'x' event interrupts

    Args:
        cigar: cigar tuples
        max_fuzzy_interupt (int): number of mismatches allowed

    """
    temp = join(cigar)
    longest_fuzzy_match = 0
    for pos, c in enumerate(temp):
        if c[0] != CIGAR.EQ:
            continue
        current_fuzzy_match = c[1]
        pos += 1
        fuzzy_count = 0
        while pos < len(temp) and fuzzy_count <= max_fuzzy_interupt:
            v, f = temp[pos]
            if v != CIGAR.EQ:
                fuzzy_count += 1
            else:
                current_fuzzy_match += f
            pos += 1
        if current_fuzzy_match > longest_fuzzy_match:
            longest_fuzzy_match = current_fuzzy_match
    return longest_fuzzy_match


def longest_exact_match(cigar):
    """
    returns the longest consecutive exact match

    Args:
        cigar (:class:`list` of :class:`tuple` of :class:`int` and :class:`int`): the cigar tuples
    """
    return longest_fuzzy_match(cigar, 0)


def score(cigar, **kwargs):
    """scoring based on sw alignment properties with gap extension penalties

    Args:
        cigar (:class:`list` of :class:`~mavis.constants.CIGAR` and :class:`int`):
          list of cigar tuple values
        MISMATCH (int): mismatch penalty
        MATCH (int): match penalty
        GAP (int): initial gap penalty
        GAP_EXTEND (int): gap extension penalty

    Returns:
        int: the score value
    """

    mismatch = kwargs.pop('MISMATCH', -1)
    match = kwargs.pop('MATCH', 2)
    gap = kwargs.pop('GAP', -4)
    gap_extend = kwargs.pop('GAP_EXTEND', -1)

    score = 0
    for v, freq in cigar:
        if v == CIGAR.EQ:
            score += match * freq
        elif v == CIGAR.X:
            score += mismatch * freq
        elif v in [CIGAR.I, CIGAR.D]:
            score += gap + gap_extend * (freq - 1)
        elif v in [CIGAR.S, CIGAR.N]:
            pass
        else:
            raise AssertionError('unexpected cigar value', v)
    return score


def match_percent(cigar):
    """
    calculates the percent of aligned bases (matches or mismatches) that are matches
    """
    matches = 0
    mismatches = 0
    for v, f in cigar:
        if v == CIGAR.EQ:
            matches += f
        elif v == CIGAR.X:
            mismatches += f
        elif v == CIGAR.M:
            raise AttributeError('cannot calculate match percent with non-specific alignments', cigar)
    if matches + mismatches == 0:
        raise AttributeError('input cigar str does not have any aligned sections (X or =)', cigar)
    else:
        return matches / (matches + mismatches)


def join(*pos):
    """
    given a number of cigar lists, joins them and merges any consecutive tuples
    with the same cigar value

    Example:
        >>> join([(1, 1), (4, 7)], [(4, 3), (2, 4)])
        [(1, 1), (4, 10), (2, 4)]
    """
    result = []
    for cigar in pos:
        for v, f in cigar:
            if len(result) > 0 and result[-1][0] == v:
                result[-1] = (v, f + result[-1][1])
            else:
                result.append((v, f))
    return result


def extend_softclipping(cigar, min_exact_to_stop_softclipping):
    """
    given some input cigar, extends softclipping if there are mismatches/insertions/deletions
    close to the end of the aligned portion. The stopping point is defined by the
    min_exact_to_stop_softclipping parameter. this function will throw an error if there is no
    exact match aligned portion to signal stop

    Args:
        original_cigar (:class:`list` of :class:`~mavis.constants.CIGAR` and :class:`int`): the input cigar
        min_exact_to_stop_softclipping (int): number of exact matches to terminate extension

    Returns:
        tuple:
            - :class:`list` of :class:`~mavis.constants.CIGAR` and :class:`int` - new cigar list
            - :class:`int` - shift from the original start position
    """
    new_cigar = []
    anchors = [i for i, (v, f) in enumerate(cigar) if v in {CIGAR.EQ, CIGAR.M} and f >= min_exact_to_stop_softclipping]
    if not anchors:
        raise AttributeError('cannot compute on this cigar as there is no stop point')
    start_anchor = min(anchors)
    end_anchor = max(anchors)
    if cigar[0][0] == CIGAR.H:
        start_anchor = 0
    if cigar[-1][0] == CIGAR.H:
        end_anchor = len(cigar)
    start_query_aligned = sum([f for v, f in cigar[:start_anchor] if v in QUERY_ALIGNED_STATES] + [0])
    start_ref_aligned = sum([f for v, f in cigar[:start_anchor] if v in REFERENCE_ALIGNED_STATES] + [0])
    end_query_aligned = sum([f for v, f in cigar[end_anchor + 1:] if v in QUERY_ALIGNED_STATES] + [0])
    new_cigar = []
    if start_query_aligned:
        new_cigar.append((CIGAR.S, start_query_aligned))
    new_cigar.extend(cigar[start_anchor:end_anchor + 1])
    if end_query_aligned:
        new_cigar.append((CIGAR.S, end_query_aligned))
    return new_cigar, start_ref_aligned


def compute(ref, alt, force_softclipping=True, min_exact_to_stop_softclipping=6):
    """
    given a ref and alt sequence compute the cigar string representing the alt

    returns the cigar tuples along with the start position of the alt relative to the ref
    """
    if not force_softclipping:
        min_exact_to_stop_softclipping = 1

    if len(ref) != len(alt):
        raise AttributeError('ref and alt must be the same length')
    cigar = []
    for r, a in zip(ref, alt):
        if r == a and r == GAP:
            pass
        elif r == GAP:
            cigar.append((CIGAR.I, 1))
        elif a == GAP:
            cigar.append((CIGAR.D, 1))
        elif DNA_ALPHABET.match(r, a):
            cigar.append((CIGAR.EQ, 1))
        else:
            cigar.append((CIGAR.X, 1))
    cigar = join(cigar)

    try:
        c, rs = extend_softclipping(cigar, min_exact_to_stop_softclipping)
        return c, rs
    except AttributeError:
        return cigar, 0


def convert_for_igv(cigar):
    """
    igv does not support the extended CIGAR values for match v mismatch

    Example:
        >>> convert_for_igv([(7, 4), (8, 1), (7, 5)])
        [(0, 10)]
    """
    result = []
    for v, f in cigar:
        if v in ALIGNED_STATES:
            v = CIGAR.M
        result.append((v, f))
    return join(result)


def alignment_matches(cigar):
    """
    counts the number of aligned bases irrespective of match/mismatch
    this is equivalent to counting all CIGAR.M
    """
    result = 0
    for v, f in cigar:
        if v in ALIGNED_STATES:
            result += f
    return result


def merge_indels(cigar):
    """
    For a given cigar tuple, merges adjacent insertions/deletions

    Example:
        >>> merge_indels([(CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.D, 4), (CIGAR.I, 2), (CIGAR.D, 2), (CIGAR.EQ, 10)])
        [(CIGAR.EQ, 10), (CIGAR.I, 5), (CIGAR.D, 6), (CIGAR.EQ, 10)]
    """
    new_cigar = cigar[:]
    for i in range(0, len(new_cigar)):
        t = i - 1  # for bubbling
        if new_cigar[i][0] == CIGAR.I:
            while t >= 0:
                if new_cigar[t][0] != CIGAR.D:
                    break
                t -= 1
            if t < i - 1:
                t += 1
                new_cigar[i], new_cigar[t] = new_cigar[t], new_cigar[i]
    new_cigar = join(new_cigar)
    return new_cigar


def hgvs_standardize_cigar(read, reference_seq):
    """
    extend alignments as long as matches are possible.
    call insertions before deletions
    """
    cigar = join(read.cigar)
    new_cigar = []
    # ensure that any del ins become ins del
    for i in range(0, len(cigar)):
        convert = False
        if cigar[i][0] == CIGAR.X:
            if i > 0:
                if cigar[i - 1][0] in [CIGAR.I, CIGAR.D, CIGAR.N]:
                    convert = True
            if i < len(cigar) - 1:
                if cigar[i + 1][0] in [CIGAR.I, CIGAR.D, CIGAR.N]:
                    convert = True
        if cigar[i][0] == CIGAR.N:
            new_cigar.append((CIGAR.D, cigar[i][1]))
        elif convert:
            new_cigar.append((CIGAR.I, cigar[i][1]))
            new_cigar.append((CIGAR.D, cigar[i][1]))
        else:
            new_cigar.append(cigar[i])
    new_cigar = merge_indels(new_cigar)
    # now we need to extend any insertions
    rpos = read.reference_start
    qpos = 0
    cigar = [new_cigar[0]]
    if cigar[0][0] in REFERENCE_ALIGNED_STATES:
        rpos += cigar[0][1]
    if cigar[0][0] in QUERY_ALIGNED_STATES:
        qpos += cigar[0][1]
    i = 1
    while i < len(new_cigar):
        if i < len(new_cigar) - 1:
            c, v = new_cigar[i]
            next_c, next_v = new_cigar[i + 1]
            prev_c, prev_v = new_cigar[i - 1]

            if c == CIGAR.I:
                qseq = read.query_sequence[qpos:qpos + v]
                if next_c == CIGAR.EQ and prev_c == CIGAR.EQ:
                    rseq = reference_seq[rpos:rpos + next_v]
                    t = 0
                    while t < next_v and rseq[t] == read.query_sequence[qpos + t]:
                        t += 1
                    if t > 0:
                        cigar.append((CIGAR.EQ, t))
                        rpos += t
                        qpos += t
                        if t == next_v:
                            del new_cigar[i + 1]
                        else:
                            new_cigar[i + 1] = next_c, next_v - t
                        continue
                elif next_c == CIGAR.D and prev_c == CIGAR.EQ:
                    # reduce the insertion and deletion by extending the alignment if possible
                    delseq = reference_seq[rpos: rpos + next_v]
                    start = 0
                    end = 0
                    shift = True
                    while max(start, end) < min(len(delseq), len(qseq)) and shift:
                        shift = False
                        if qseq[start] == delseq[start]:
                            start += 1
                            shift = True
                        if qseq[-1 - end] == delseq[-1 - end]:
                            end += 1
                            shift = True
                    if start:
                        cigar.append((CIGAR.EQ, start))
                        if start < next_v:
                            new_cigar[i + 1] = (next_c, next_v - start)
                        else:
                            del new_cigar[i + 1]
                        if start < v:
                            new_cigar[i] = (c, v - start)
                        else:
                            del new_cigar[i]
                        qpos += start
                        rpos += start
                        continue
                    elif end:
                        if end < v:
                            cigar.append((c, v - end))
                        if end < next_v:
                            cigar.append((next_c, next_v - end))
                        cigar.append((CIGAR.EQ, end))
                        qpos += v
                        rpos += next_v
                        i += 2
                        continue
                qpos += v

            elif c == CIGAR.D:
                rseq = reference_seq[rpos:rpos + v]
                if next_c == CIGAR.EQ and prev_c == CIGAR.EQ:
                    qseq = read.query_sequence[qpos:qpos + next_v]
                    t = 0
                    while t < next_v and qseq[t] == reference_seq[rpos + t]:
                        t += 1
                    if t > 0:
                        cigar.append((CIGAR.EQ, t))
                        qpos += t
                        rpos += t
                        if t == next_v:
                            del new_cigar[i + 1]
                        else:
                            new_cigar[i + 1] = next_c, next_v - t
                        continue
                rpos += v
            else:
                if c in QUERY_ALIGNED_STATES:
                    qpos += v
                if c in REFERENCE_ALIGNED_STATES:
                    rpos += v
        cigar.append(new_cigar[i])
        i += 1
    return join(cigar)


def convert_string_to_cigar(string):
    """
    Given a cigar string, converts it to the appropriate cigar tuple

    Example:
        >>> convert_string_to_cigar('8M2I1D9X')
        [(CIGAR.M, 8), (CIGAR.I, 2), (CIGAR.D, 1), (CIGAR.X, 9)]
    """
    patt = r'(\d+(\D))'
    cigar = [m[0] for m in re.findall(patt, string)]
    cigar = [(CIGAR[match[-1]] if match[-1] != '=' else CIGAR.EQ, int(match[:-1])) for match in cigar]
    return cigar


def convert_cigar_to_string(cigar):
    return ''.join(['{}{}'.format(f, CIGAR.reverse(s) if s != CIGAR.EQ else '=') for s, f in cigar])


def merge_internal_events(cigar, inner_anchor=10, outer_anchor=10):
    """
    merges events (insertions, deletions, mismatches) within a cigar if they are
    between exact matches on either side (anchors) and separated by less exact
    matches than the given parameter

    does not merge two mismatches, must contain a deletion/insertion

    Args:
        cigar (list): a list of tuples of cigar states and counts
        inner_anchor (int): minimum number of consecutive exact matches separating events
        outer_anchor (int): minimum consecutively aligned exact matches to anchor an end for merging

    Returns:
        list: new list of cigar tuples with merged events

    Example:
        >>> merge_internal_events([(CIGAR.EQ, 10), (CIGAR.X, 1), (CIGAR.EQ, 2), (CIGAR.D, 1), (CIGAR.EQ, 10)])
        [(CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.D, 4), (CIGAR.EQ, 10)]
    """
    read_cigar = join(cigar)
    prefix = []
    # get the initial anchors
    exact_match_pos = []
    for i, tup in enumerate(read_cigar):
        if tup[0] in [CIGAR.EQ, CIGAR.M] and tup[1] >= outer_anchor:
            exact_match_pos.append((i, tup[1]))

    if len(exact_match_pos) < 2:
        return read_cigar

    ppos = exact_match_pos[0][0]
    prefix = read_cigar[0: ppos + 1]
    spos = exact_match_pos[-1][0]
    suffix = read_cigar[spos:None]
    read_cigar = read_cigar[len(prefix):len(suffix) * -1]
    new_cigar = prefix[:]
    for state, count in read_cigar:
        last_state, last_value = new_cigar[-1]

        if state == CIGAR.X:
            if last_state in EVENT_STATES - {CIGAR.X}:  # any event that is not a mismatch
                new_cigar.extend([(CIGAR.I, count), (CIGAR.D, count)])
            else:
                new_cigar.append((state, count))
        elif state in EVENT_STATES:
            if last_state == CIGAR.X:  # last event was a mismatch, convert the mismatch to an indel
                new_cigar[-1] = CIGAR.I, last_value
                new_cigar.append((CIGAR.D, last_value))
            # count the inner block anchor size
            last_event = len(new_cigar)
            for i, (new_state, new_count) in enumerate(new_cigar[::-1]):
                if new_state in EVENT_STATES:
                    last_event = len(new_cigar) - i - 1
                elif last_event != len(new_cigar):
                    break
            anchor = sum([v for c, v in new_cigar[last_event + 1:]] + [0])
            # if the anchor is too small, merge it
            if last_event >= len(prefix) and anchor < inner_anchor:
                temp = new_cigar[last_event:]
                new_cigar = new_cigar[:last_event]
                for new_state, new_count in temp:
                    if new_state in {CIGAR.EQ, CIGAR.M, CIGAR.X}:
                        new_cigar.extend([(CIGAR.D, new_count), (CIGAR.I, new_count)])
                    else:
                        new_cigar.append((new_state, new_count))
            new_cigar.append((state, count))
        elif state in [CIGAR.EQ, CIGAR.M]:
            if count >= inner_anchor or last_state not in EVENT_STATES - {CIGAR.X}:
                new_cigar.append((state, count))
            else:
                if last_state == CIGAR.X:
                    new_cigar[-1] = CIGAR.I, last_value
                    new_cigar.append((CIGAR.D, last_value))
                new_cigar.append((CIGAR.I, count))
                new_cigar.append((CIGAR.D, count))
        else:
            new_cigar.append((state, count))
    new_cigar.extend(suffix)
    return merge_indels(new_cigar)
