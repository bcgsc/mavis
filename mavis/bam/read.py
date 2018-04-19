from copy import copy
import itertools
import re
import subprocess

import pysam
from Bio.Data import IUPACData as iupac

from . import cigar as _cigar
from .cigar import EVENT_STATES, QUERY_ALIGNED_STATES, REFERENCE_ALIGNED_STATES, convert_cigar_to_string
from ..constants import CIGAR, DNA_ALPHABET, ORIENT, READ_PAIR_TYPE, STRAND, SVTYPE, NA_MAPPING_QUALITY
from ..interval import Interval


class SamRead(pysam.AlignedSegment):
    """
    Subclass to extend the pysam.AlignedSegment class adding some utility methods and convenient representations

    Allows next_reference_name and reference_name to be set directly so that is does not depend on a bam header
    """

    def __init__(self, reference_name=None, next_reference_name=None, alignment_score=None, **kwargs):
        pysam.AlignedSegment.__init__(self)
        self._reference_name = reference_name
        self._next_reference_name = next_reference_name
        self.alignment_score = alignment_score
        self.mapping_quality = NA_MAPPING_QUALITY
        self.alignment_rank = None
        self._key = None
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def set_key(self):
        """
        Warning:
            Using this method sets the _key attribute which is used for comparison and hashing. If you alter
            this attribute while items are in a hashed state it may lead to unexpected results such as duplicates
            of a single object within a set
        """
        self._key = (self.query_name, self.query_sequence, self.reference_id, self.reference_start, self.is_supplementary)

    def key(self):
        """
        uses a stored _key attribute, if available. This is to avoid the hash changing if the reference start (for example)
        is changed but also allow this attribute to be used and calculated for non SamRead objects

        This way to change the hash behaviour the user must be explicit and use the set_key method
        """
        if hasattr(self, '_key') and self._key is not None:
            return self._key
        return (self.query_name, self.query_sequence, self.reference_id, self.reference_start, self.is_supplementary)

    def alignment_id(self):
        return '{}:{}[{}]{}'.format(self.reference_name, self.reference_start, self.query_name, convert_cigar_to_string(self.cigar))

    @classmethod
    def copy(cls, pysamread):
        return cls.copy_onto(pysamread)

    @classmethod
    def copy_onto(cls, pysamread, copyread=None):
        cp = cls() if copyread is None else copyread
        cp._reference_name = pysamread.reference_name
        cp.query_sequence = pysamread.query_sequence
        cp.reference_start = pysamread.reference_start
        cp.reference_id = pysamread.reference_id
        cp.cigar = pysamread.cigar[:]
        cp.query_name = pysamread.query_name
        cp.mapping_quality = pysamread.mapping_quality
        cp.query_qualities = pysamread.query_qualities
        cp.template_length = pysamread.template_length
        try:
            cp.alignment_rank = pysamread.alignment_rank
        except AttributeError:
            pass
        cp.set_tags(pysamread.get_tags())
        cp.flag = pysamread.flag

        if pysamread.is_paired:
            cp.next_reference_id = pysamread.next_reference_id
            cp.next_reference_start = pysamread.next_reference_start
            cp._next_reference_name = pysamread.next_reference_name
        try:
            cp.alignment_score = pysamread.alignment_score
        except AttributeError:
            pass
        cp.set_key()
        return cp

    def __copy__(self):
        return self.__class__.copy(self)

    @property
    def reference_name(self):
        return self._reference_name

    @property
    def next_reference_name(self):
        return self._next_reference_name

    def deletion_sequences(self, reference_genome):
        """returns the reference sequences for all deletions"""
        rpos = self.reference_start
        result = []
        for state, freq in self.cigar:
            if state in REFERENCE_ALIGNED_STATES:
                if state not in QUERY_ALIGNED_STATES:
                    result.append(reference_genome[self.reference_name].seq[rpos:rpos + freq])
                rpos += freq
        return result

    def insertion_sequences(self):
        """returns the inserted sequence for all insertions"""
        qpos = 0
        result = []
        for state, freq in self.cigar:
            if state in QUERY_ALIGNED_STATES:
                if state not in REFERENCE_ALIGNED_STATES:
                    result.append(self.query_sequence[qpos:qpos + freq])
                qpos += freq
        return result

    def __eq__(self, other):
        return self.key() == SamRead.key(other)

    def __hash__(self):
        return hash(self.key())


def pileup(reads, filter_func=None):
    """
    For a given set of reads generate a pileup of all reads (excluding those for which the filter_func returns True)

    Args:
        reads (iterable of pysam.AlignedSegment): reads to pileup
        filter_func (callable): function which takes in a  read and returns True if it should be ignored and False otherwise

    Returns:
        iterable of tuple of int and int: tuples of genomic position and read count at that position

    Note:
        returns positions using 1-based indexing
    """
    hist = {}  # genome position => frequency count
    for read in reads:
        if filter_func and filter_func(read):
            continue
        for pos in read.get_reference_positions():
            hist[pos + 1] = hist.get(pos + 1, 0) + 1
    return sorted(hist.items())


def map_ref_range_to_query_range(read, ref_range):
    """
    Args:
        ref_range (Interval): 1-based inclusive
        read (pysam.AlignedSegment): read used for the mapping
    Returns:
        Interval: 1-based inclusive range
    """
    rpos = read.reference_start
    qpos = 0
    qstart = None
    qend = None
    for state, value in read.cigar:
        for i in range(0, value):
            if state in QUERY_ALIGNED_STATES:
                qpos += 1
            if state in REFERENCE_ALIGNED_STATES:
                rpos += 1
            if qstart is None and ref_range.start <= rpos:
                qstart = qpos
            if ref_range.end >= rpos:
                qend = qpos
    if qstart is None or qend is None:
        raise ValueError('reference range is not mapped by input read', ref_range, read.reference_start, read.cigar)
    return Interval(qstart, qend)


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


def calculate_alignment_score(read, consec_bonus=1):
    """
    calculates a score for comparing alignments

    Args:
        read (pysam.AlignedSegment): the input read

    Returns:
        float: the score
    """
    score = 0
    qlen = read.reference_end - read.reference_start
    max_score = qlen + (qlen - 1) * consec_bonus
    for c, v in read.cigar:
        if c == CIGAR.M:
            raise ValueError('cannot calculate the alignment score if mismatch v match has not been specified')
        elif c == CIGAR.EQ:
            score += v + (v - 1) * consec_bonus
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
        raise AttributeError('cannot overlap on an empty sequence: len(ref)={}, len(seq)={}'.format(len(ref), len(seq)))
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

        cigar = _cigar.join(cigar)
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])

        qstart = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]

        a = SamRead(
            query_sequence=str(seq),
            reference_start=ref_start + qstart,
            cigar=cigar
        )
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


def sequence_complexity(seq):
    """
    basic measure of sequence complexity
    """
    if not seq:
        return 0
    hist = {c: 0 for c in iupac.unambiguous_dna_letters}
    for ambig_base in seq.upper():
        values = iupac.ambiguous_dna_values[ambig_base]
        for base in values:  # ignore N's etc
            hist[base] += 1 / len(values)
    total = sum(hist.values())
    scores = [(hist[base1] + hist[base2]) / total for base1, base2 in itertools.combinations(iupac.unambiguous_dna_letters, 2)]
    return min(scores)
