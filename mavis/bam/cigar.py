"""
holds methods related to processing cigar tuples. Cigar tuples are generally
an iterable list of tuples where the first element in each tuple is the
CIGAR value (i.e. 1 for an insertion), and the second value is the frequency
"""
from ..constants import CIGAR, DNA_ALPHABET, GAP


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
    temp = []
    offset = 0

    ref_pos = read.reference_start
    seq_pos = 0

    for cigar_value, freq in read.cigar:
        if cigar_value in [CIGAR.S, CIGAR.I]:
            temp.append((cigar_value, freq))
            seq_pos += freq
        elif cigar_value in [CIGAR.D, CIGAR.N]:
            temp.append((cigar_value, freq))
            ref_pos += freq
        elif cigar_value in [CIGAR.M, CIGAR.X, CIGAR.EQ]:
            for offset in range(0, freq):
                if DNA_ALPHABET.match(ref[ref_pos], read.query_sequence[seq_pos]):
                    if len(temp) == 0 or temp[-1][0] != CIGAR.EQ:
                        temp.append((CIGAR.EQ, 1))
                    else:
                        temp[-1] = (CIGAR.EQ, temp[-1][1] + 1)
                else:
                    if len(temp) == 0 or temp[-1][0] != CIGAR.X:
                        temp.append((CIGAR.X, 1))
                    else:
                        temp[-1] = (CIGAR.X, temp[-1][1] + 1)
                ref_pos += 1
                seq_pos += 1
        else:
            raise NotImplementedError('unexpected CIGAR value {0} is not supported currently'.format(cigar_value))
    assert(sum([x[1] for x in temp]) == sum(x[1] for x in read.cigar))
    return temp


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
        cigar (:class:`list` of :class:`~structural_variant.constants.CIGAR` and :class:`int`):
          list of cigar tuple values
        MISMATCH (int): mismatch penalty
        MATCH (int): match penalty
        GAP (int): initial gap penalty
        GAP_EXTEND (int): gap extension penalty

    Returns:
        int: the score value
    """

    MISMATCH = kwargs.pop('MISMATCH', -1)
    MATCH = kwargs.pop('MATCH', 2)
    GAP = kwargs.pop('GAP', -4)
    GAP_EXTEND = kwargs.pop('GAP_EXTEND', -1)

    score = 0
    for v, freq in cigar:
        if v == CIGAR.EQ:
            score += MATCH * freq
        elif v == CIGAR.X:
            score += MISMATCH * freq
        elif v in [CIGAR.I, CIGAR.D]:
            score += GAP + GAP_EXTEND * (freq - 1)
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
        original_cigar (:class:`list` of :class:`~structural_variant.constants.CIGAR` and :class:`int`): the input cigar
        min_exact_to_stop_softclipping (int): number of exact matches to terminate extension

    Returns:
        (list of CIGAR and int, int): the new cigar string and a number representing the shift from the original
            start position
    """
    ref_start_shift = 0
    # determine how far to scoop for the front softclipping
    new_cigar = []
    temp = [(v, f) for v, f in cigar if (v in [CIGAR.EQ, CIGAR.M] and f >= min_exact_to_stop_softclipping)]
    if len(temp) == 0:
        raise AttributeError('cannot compute on this cigar as there is no stop point')

    match_satisfied = False
    for v, f in cigar:
        if v in [CIGAR.M, CIGAR.EQ] and f >= min_exact_to_stop_softclipping:
            match_satisfied = True
        if match_satisfied and (len(new_cigar) == 0 or new_cigar[-1][0] != CIGAR.S):
            new_cigar.append((v, f))
        elif match_satisfied:  # first after SC
            if v in [CIGAR.D, CIGAR.X, CIGAR.N]:
                ref_start_shift += f
                new_cigar.append((CIGAR.S, f))
            elif v == CIGAR.I:
                new_cigar.append((CIGAR.S, f))
            else:
                new_cigar.append((v, f))
        else:
            if v in [CIGAR.D, CIGAR.N]:
                ref_start_shift += f
                pass
            elif v in [CIGAR.I, CIGAR.S]:
                new_cigar.append((CIGAR.S, f))
            else:
                new_cigar.append((CIGAR.S, f))
                ref_start_shift += f
    cigar = new_cigar[::-1]

    new_cigar = []

    match_satisfied = False
    for v, f in cigar:
        if v in [CIGAR.M, CIGAR.EQ] and f >= min_exact_to_stop_softclipping:
            match_satisfied = True
        if match_satisfied and (len(new_cigar) == 0 or new_cigar[-1][0] != CIGAR.S):
            new_cigar.append((v, f))
        elif match_satisfied:  # first after SC
            if v in [CIGAR.D, CIGAR.X, CIGAR.N]:
                new_cigar.append((CIGAR.S, f))
            elif v == CIGAR.I:
                new_cigar.append((CIGAR.S, f))
            else:
                new_cigar.append((v, f))
        else:
            if v in [CIGAR.D, CIGAR.N]:
                pass
            elif v in [CIGAR.I, CIGAR.S]:
                new_cigar.append((CIGAR.S, f))
            else:
                new_cigar.append((CIGAR.S, f))
    new_cigar.reverse()
    return new_cigar, ref_start_shift


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
        if v in [CIGAR.X, CIGAR.EQ]:
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
        if v in [CIGAR.X, CIGAR.EQ, CIGAR.M]:
            result += f
    return result


def smallest_nonoverlapping_repeat(s):
    for repsize in range(1, len(s) + 1):
        if len(s) % repsize == 0:
            substrings = [s[i:i + repsize] for i in range(0, len(s), repsize)]
            if len(set(substrings)) == 1:
                return substrings[0]
    return s


def hgvs_standardize_cigar(read, reference_seq):
    """
    extend alignments as long as matches are possible.
    call insertions before deletions
    """
    ci = 0
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
    # now sort any consecutive ins/dels so that that insertions come first
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
    # now we need to extend any insertions
    rpos = read.reference_start
    qpos = 0
    cigar = []
    i = 0
    while i < len(new_cigar):
        if i < len(new_cigar) - 1:
            c, v = new_cigar[i]
            next_c, next_v = new_cigar[i + 1]
            
            if c == CIGAR.I:
                qpos += v
                qseq = read.query_sequence[qpos - v:qpos]
                qrep = smallest_nonoverlapping_repeat(qseq)
                if next_c == CIGAR.EQ and next_v >= len(qrep):
                    rseq = reference_seq[rpos:rpos + next_v]
                    t = 0
                    while t + len(qrep) <= next_v and rseq[t:t + len(qrep)] == qrep:
                        t += len(qrep)
                    if t > 0:
                        cigar.append((CIGAR.EQ, t))
                        rpos += t
                        if t == next_v:
                            del new_cigar[i + 1]
                        else:
                            new_cigar[i + 1] = next_c, next_v - t
                        continue
            elif c == CIGAR.D:
                rpos += v
                rseq = reference_seq[rpos - v:rpos]
                rrep = smallest_nonoverlapping_repeat(rseq)
                if next_c == CIGAR.EQ and next_v >= len(rrep):
                    qseq = read.query_sequence[qpos:qpos + next_v]
                    t = 0
                    while t + len(rrep) <= next_v and qseq[t:t + len(rrep)] == rrep:
                        t += len(rrep)
                    if t > 0:
                        cigar.append((CIGAR.EQ, t))
                        qpos += t
                        if t == next_v:
                            del new_cigar[i + 1]
                        else:
                            new_cigar[i + 1] = next_c, next_v - t
                        continue
            elif c == CIGAR.S:
                qpos += v
            elif c != CIGAR.H:
                qpos += v
                rpos += v
        cigar.append(new_cigar[i])
        i += 1
    return join(cigar)
