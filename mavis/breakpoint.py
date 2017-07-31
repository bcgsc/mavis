from __future__ import division

from copy import copy as sys_copy
from .constants import ORIENT, STRAND, COLUMNS, CIGAR, SVTYPE, reverse_complement, DNA_ALPHABET
from .error import *
from .interval import Interval
import TSV
import re
import itertools


class Breakpoint(Interval):
    """
    class for storing information about a SV breakpoint
    coordinates are given as 1-indexed
    """
    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand)

    def __init__(self, chr, start, end=None, orient=ORIENT.NS, strand=STRAND.NS, seq=None):
        """
        Args:
            chr (str): the chromosome
            start (int): the genomic position of the breakpoint
            end (int): if the breakpoint is uncertain (a range) then specify the end of the range here
            orient (ORIENT): the orientation (which side is retained at the break)
            strand (STRAND): the strand
            seq (str): the seq

        Examples:
            >>> Breakpoint('1', 1, 2)
            >>> Breakpoint('1', 1)
            >>> Breakpoint('1', 1, 2, 'R', )
            >>> Breakpoint('1', 1, orient='R')
        """
        from .annotate.base import ReferenceName

        Interval.__init__(self, start, end)
        self.orient = ORIENT.enforce(orient)
        self.chr = ReferenceName(chr)
        self.strand = STRAND.enforce(strand)
        self.seq = seq

    def __repr__(self):
        strand = '' if self.strand == STRAND.NS else self.strand
        orient = '' if self.orient == ORIENT.NS else self.orient

        return 'Breakpoint({0}:{1}{2}{3}{4})'.format(
            self.chr,
            self.start,
            '-' + str(self.end) if self.end != self.start else '',
            orient,
            strand
        )

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def to_dict(self):
        d = {
            'chr': self.chr,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'seq': self.seq,
            'orientation': self.orient,
            'type': self.__class__.__name__
        }
        return d


class BreakpointPair:
    """
    """

    def __getattr__(self, attr):
        data = object.__getattribute__(self, 'data')
        try:
            return data[COLUMNS[attr]]
        except KeyError:
            try:
                return data[attr]
            except KeyError:
                pass
        raise AttributeError(attr)

    def __getitem__(self, index):
        try:
            index = int(index)
        except ValueError:
            raise IndexError('index input accessor must be an integer', index)
        if index == 0:
            return self.break1
        elif index == 1:
            return self.break2
        raise IndexError('index input accessor is out of bounds: 1 or 2 only', index)

    def __hash__(self):
        return hash((self.break1.key, self.break2.key, self.opposing_strands, self.stranded, self.untemplated_seq))

    def __eq__(self, other):
        for attr in ['break1', 'break2', 'opposing_strands', 'stranded', 'untemplated_seq']:
            if not hasattr(other, attr):
                return False
            elif getattr(self, attr) != getattr(other, attr):
                return False
        return True

    @property
    def interchromosomal(self):
        """:class:`bool`: True if the breakpoints are on different chromosomes, False otherwise"""
        if self.break1.chr == self.break2.chr:
            return False
        return True

    def copy(self):
        temp = sys_copy(self)
        temp.break1 = sys_copy(self.break1)
        temp.break2 = sys_copy(self.break2)
        return temp

    def __init__(self, b1, b2, stranded=False, opposing_strands=None, untemplated_seq=None, data=None, **kwargs):
        """
        Args:
            b1 (Breakpoint): the first breakpoint
            b2 (Breakpoint): the second breakpoint
            stranded (bool): if not stranded then +/- is equivalent to -/+
            opposing_strands (bool): are the strands at the breakpoint opposite? i.e. +/- instead of +/+
            untemplated_seq (str): seq between the breakpoints that is not part of either breakpoint
            data (dict): optional dictionary of attributes associated with this pair

        Note:
            untemplated_seq should always be given wrt to the positive/forward reference strand

        Example:
            >>> BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
            >>> BreakpointPair(Breakpoint('1', 1, strand='+'), Breakpoint('1', 9999, strand='-'))
        """

        if b1.key > b2.key:
            self.break1 = b2
            self.break2 = b1
        else:
            self.break1 = b1
            self.break2 = b2
        self.stranded = stranded
        self.opposing_strands = opposing_strands
        # between break1 and break2 not in either
        self.untemplated_seq = untemplated_seq
        self.data = {}
        if data is not None:
            self.data.update(data)
            conflicts = set(data.keys()) & set(kwargs.keys())
            if len(conflicts) > 0:
                raise TypeError('data got multiple values for data elements:', conflicts)
        self.data.update(kwargs)

        if self.break1.strand != STRAND.NS and self.break2.strand != STRAND.NS:
            opposing = self.break1.strand != self.break2.strand
            if self.opposing_strands is None:
                self.opposing_strands = opposing
            elif self.opposing_strands != opposing:
                raise AssertionError(
                    'conflict in input arguments, opposing_strands must agree with input breakpoints'
                    ' when the strand has been specified'
                )
        if self.break1.orient != ORIENT.NS and self.break2.orient != ORIENT.NS:
            if self.opposing_strands is not None:
                if (self.break1.orient == self.break2.orient and not self.opposing_strands) \
                        or (self.break1.orient != self.break2.orient and self.opposing_strands):
                    raise InvalidRearrangement(
                        'invalid breakpoint pair cannot form a valid combination', b1, b2, self.opposing_strands)

        if self.opposing_strands is None:
            raise NotSpecifiedError('must specify if opposing_strands')
        if self.stranded and STRAND.NS in [self.break1.strand, self.break2.strand]:
            raise NotSpecifiedError('if stranded is specified, breakpoint strands cannot be \'not specified\'')

        # try classifying to make sure it's a valid combination
        BreakpointPair.classify(self)

    def __str__(self):
        return 'BPP({}, {}{}{})'.format(
            str(self.break1),
            str(self.break2),
            ', opposing={}'.format(self.opposing_strands) if not self.stranded else '',
            ', seq=' + repr(self.untemplated_seq) if self.untemplated_seq is not None else ''
        )

    def flatten(self):
        """
        returns the key-value self for the breakpoint self information as
        can be written directly as a TSV row
        """
        row = {}
        row.update(self.data)
        temp = {
            COLUMNS.break1_chromosome: self.break1.chr,
            COLUMNS.break1_position_start: self.break1.start,
            COLUMNS.break1_position_end: self.break1.end,
            COLUMNS.break1_orientation: self.break1.orient,
            COLUMNS.break1_strand: self.break1.strand,
            COLUMNS.break2_chromosome: self.break2.chr,
            COLUMNS.break2_position_start: self.break2.start,
            COLUMNS.break2_position_end: self.break2.end,
            COLUMNS.break2_orientation: self.break2.orient,
            COLUMNS.break2_strand: self.break2.strand,
            COLUMNS.opposing_strands: self.opposing_strands,
            COLUMNS.stranded: self.stranded,
            COLUMNS.untemplated_seq: self.untemplated_seq,
            COLUMNS.break1_seq: self.break1.seq,
            COLUMNS.break2_seq: self.break2.seq
        }
        for c in temp:
            temp[c] = str(temp[c])
        row.update(temp)
        return row

    @classmethod
    def classify(cls, pair):
        """
        uses the chr, orientations and strands to determine the
        possible structural_variant types that this pair could support

        Args:
            pair (BreakpointPair): the pair to classify
        Returns:
            :class:`list` of :any:`SVTYPE`: a list of possible SVTYPE

        Example:
            >>> bpp = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
            >>> BreakpointPair.classify(bpp)
            ['inversion']
            >>> bpp = BreakpointPair(Breakpoint('1', 1, orient='L'), Breakpoint('1', 9999, orient='R'), opposing_strands=False)
            >>> BreakpointPair.classify(bpp)
            ['deletion', 'insertion']

        see :ref:`related theory documentation <theory-classifying-events>`
        """
        if pair.break1.chr == pair.break2.chr:  # intrachromosomal
            if pair.opposing_strands:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.RIGHT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.LEFT):
                    raise InvalidRearrangement(pair)
                return [SVTYPE.INV]
            else:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.LEFT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.RIGHT):
                    raise InvalidRearrangement(pair)
                elif pair.break1.orient == ORIENT.LEFT or pair.break2.orient == ORIENT.RIGHT:
                    return [SVTYPE.DEL, SVTYPE.INS]
                elif pair.break1.orient == ORIENT.RIGHT or pair.break2.orient == ORIENT.LEFT:
                    return [SVTYPE.DUP]
        else:  # interchromosomal
            if pair.opposing_strands:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.RIGHT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.LEFT):
                    raise InvalidRearrangement(pair)
                return [SVTYPE.ITRANS]
            else:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.LEFT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.RIGHT):
                    raise InvalidRearrangement(pair)
                return [SVTYPE.TRANS]

    @classmethod
    def call_breakpoint_pair(cls, read1, read2=None, REFERENCE_GENOME=None):
        """
        calls a set of breakpoints from a single or a pair of pysam style read(s)

        Args:
            read1 (pysam.AlignedSegment): the first read
            read2 (pysam.AlignedSegment): the second read

        Returns:
            BreakpointPair: the newly called breakpoint pair from the contig

        .. todo::

            return multiple events not just the major event
        """
        if read2 is None:
            return cls._call_from_single_contig(read1, REFERENCE_GENOME=REFERENCE_GENOME)
        return cls._call_from_paired_contig(read1, read2, REFERENCE_GENOME=REFERENCE_GENOME)

    @classmethod
    def _call_from_single_contig(cls, read, REFERENCE_GENOME=None):
        """
        calls a set of breakpoints from a pysam style read using the cigar values

        .. todo::

            return multiple events not just the major event

        Returns:
            BreakpointPair: the newly called breakpoint pair from the contig

        Raises:
            UserWarning: if the contig does not contain insertions/deletions

        if the duplicated sequence is longer than the untemplated sequence then this should
        be called as a duplication and not an insertion
        """
        read_events = []
        for i, t in enumerate(read.cigar):
            v, f = t
            if v not in [CIGAR.I, CIGAR.D, CIGAR.N]:
                continue
            if i > 0 and read.cigar[i - 1][0] in [CIGAR.I, CIGAR.D, CIGAR.N]:
                start, last, size = read_events[-1]
                read_events[-1] = (start, i, size + f)
            else:
                read_events.append((i, i, f))

        if len(read_events) == 0:
            raise UserWarning(
                'Cannot call event breakpoints from single contig. Contig does not contain any ins/del events')
        biggest = max(read_events, key=lambda x: x[2])

        first_breakpoint = read.reference_start - 1
        second_breakpoint = first_breakpoint
        untemplated_seq = ''

        seq_pos = 0
        seq_first_stop = -1
        seq_second_start = -1

        for i, t in enumerate(read.cigar):
            v, f = t
            if i < biggest[0]:  # before the event
                if v in [CIGAR.I, CIGAR.S]:
                    seq_pos += f
                elif v in [CIGAR.D, CIGAR.N]:
                    first_breakpoint += f
                    second_breakpoint += f
                else:
                    first_breakpoint += f
                    second_breakpoint += f
                    seq_pos += f
                seq_first_stop = seq_pos
            elif i >= biggest[0] and i <= biggest[1]:  # inside the event
                if v in [CIGAR.I, CIGAR.S]:
                    if v == CIGAR.I:
                        untemplated_seq += read.query_sequence[seq_pos:seq_pos + f]
                    seq_pos += f
                elif v in [CIGAR.D, CIGAR.N]:
                    second_breakpoint += f
                else:
                    second_breakpoint += f
                    seq_pos += f
            else:   # after the event
                seq_second_start = seq_pos
                break

        break1 = Breakpoint(
            read.reference_name,
            first_breakpoint + 1,  # 1-based coordinates
            orient=ORIENT.LEFT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS,
            seq=read.query_sequence[0:seq_first_stop]
        )
        break2 = Breakpoint(
            read.reference_name,
            second_breakpoint + 2,  # need to add to be past the event
            orient=ORIENT.RIGHT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS,
            seq=read.query_sequence[seq_second_start:]
        )
        # determine if the inserted sequence is a repeat of the aligned portion
        # assume that events with deletions cannot be duplications
        if first_breakpoint == second_breakpoint and REFERENCE_GENOME:  # insertion or duplication
            refseq = REFERENCE_GENOME[read.reference_name].seq[
                first_breakpoint - len(untemplated_seq) + 1:first_breakpoint + 1]
            refseq = str(refseq)
            midpoint = len(untemplated_seq) // 2 + 1  # more than the untemplated
            for dup_len in range(len(untemplated_seq), midpoint - 1, -1):
                subseq = untemplated_seq[0:dup_len]
                if subseq == refseq[-1 * dup_len:]:
                    pos = first_breakpoint - len(subseq) + 1
                    break1 = Breakpoint(
                        read.reference_name,
                        pos + 1,
                        orient=ORIENT.RIGHT,
                        strand=STRAND.NEG if read.is_reverse else STRAND.POS,
                        seq=read.query_sequence[pos:]
                    )
                    break2 = Breakpoint(
                        read.reference_name,
                        first_breakpoint + 1,
                        orient=ORIENT.LEFT,
                        strand=STRAND.NEG if read.is_reverse else STRAND.POS,
                        seq=read.query_sequence[:first_breakpoint + 1]
                    )
                    untemplated_seq = untemplated_seq[dup_len:]
                    break

        return BreakpointPair(break1, break2, opposing_strands=False, untemplated_seq=untemplated_seq)

    @classmethod
    def _call_from_paired_contig(cls, read1, read2, REFERENCE_GENOME=None):
        """
        calls a set of breakpoints from a pair of pysam style reads using their softclipping to
        find the breakpoints and orientations
        if the first read and the second read have overlapping range on their mutual query sequence
        then the breakpoint is called as is on the first and shifted on the second

        .. todo::
            return multiple events not just the major event
        """
        from .align import query_coverage_interval
        if read1.reference_id > read2.reference_id:
            read1, read2 = (read2, read1)
        elif read1.reference_id == read2.reference_id and read1.reference_start > read2.reference_start:
            read1, read2 = (read2, read1)

        r1_qci = query_coverage_interval(read1)
        r2_qci = query_coverage_interval(read2)

        if read1.is_reverse == read2.is_reverse:
            assert(read1.query_sequence == read2.query_sequence)
        else:
            assert(read1.query_sequence == reverse_complement(read2.query_sequence))
            l = len(read1.query_sequence) - 1
            r2_qci = Interval(l - r2_qci.end, l - r2_qci.start)
        b1 = None
        b2 = None

        o1 = ORIENT.NS
        o2 = ORIENT.NS
        s1 = STRAND.NEG if read1.is_reverse else STRAND.POS
        s2 = STRAND.NEG if read2.is_reverse else STRAND.POS

        # <==============
        r1_st = read1.cigar[0][1] if read1.cigar[0][0] == CIGAR.S else 0
        r1_end = read1.cigar[-1][1] if read1.cigar[-1][0] == CIGAR.S else 0

        if r1_st > r1_end:
            o1 = ORIENT.RIGHT
        elif r1_end > r1_st:
            o1 = ORIENT.LEFT

        r2_st = read2.cigar[0][1] if read2.cigar[0][0] == CIGAR.S else 0
        r2_end = read2.cigar[-1][1] if read2.cigar[-1][0] == CIGAR.S else 0

        if r2_st > r2_end:
            o2 = ORIENT.RIGHT
        elif r2_end > r2_st:
            o2 = ORIENT.LEFT

        if o1 == ORIENT.NS or o2 == ORIENT.NS:
            raise AssertionError(
                'read does not have softclipping on either end and cannot therefore determine orientation',
                read1.cigar, read2.cigar)

        if o1 == ORIENT.LEFT:  # ========++++
            b1 = Breakpoint(
                read1.reference_name,
                read1.reference_end,  # 1-based from 0-based exclusive end so don't need to -1
                strand=s1,
                orient=o1,
                seq=read1.query_alignment_sequence
            )
        else:
            b1 = Breakpoint(
                read1.reference_name,
                read1.reference_start + 1,  # 1-based from 0-based
                strand=s1,
                orient=o1,
                seq=read1.query_alignment_sequence
            )

        overlap = 0
        if Interval.overlaps(r1_qci, r2_qci):
            # adjust the second read to remove the overlapping query region
            overlap = len(r1_qci & r2_qci)

        if o2 == ORIENT.RIGHT:
            b2 = Breakpoint(
                read2.reference_name,
                read2.reference_start + overlap + 1,  # 1-based from 0-based
                strand=s2,
                orient=o2,
                seq=read2.query_alignment_sequence[overlap:]
            )
        else:
            s = read2.query_alignment_sequence
            if overlap > 0:
                s = read2.query_alignment_sequence[:-1 * overlap]
            b2 = Breakpoint(
                read2.reference_name,
                read2.reference_end - overlap,  # 1-based from 0-based exclusive end so don't need to -1
                strand=s2,
                orient=o2,
                seq=s
            )
        # now check for untemplated sequence
        untemplated_seq = ''
        d = Interval.dist(r1_qci, r2_qci)
        if d > 0:  # query coverage for read1 is after query coverage for read2
            untemplated_seq = read1.query_sequence[r2_qci[1] + 1:r1_qci[0]]
        elif d < 0:  # query coverage for read2 is after query coverage for read1
            untemplated_seq = read1.query_sequence[r1_qci[1] + 1:r2_qci[0]]
        else:  # query coverage overlaps
            pass

        return BreakpointPair(b1, b2, untemplated_seq=untemplated_seq)

    def breakpoint_sequence_homology(self, REFERENCE_GENOME):
        """
        for a given set of breakpoints matches the sequence opposite the partner breakpoint
        this sequence comparison is done with reference to a reference genome and does not
        use novel or untemplated sequence in the comparison. For this reason, insertions
        will never return any homologous sequence

        ::

            small duplication event CTT => CTTCTT

            GATACATTTCTTCTTGAAAA reference
            ---------<========== first breakpoint
            ===========>-------- second breakpoint
            ---------CT-CT------ first break homology
            -------TT-TT-------- second break homology

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by template/chr name

        Returns:
            tuple:
                - :class:`str` - homologous sequence at the first breakpoint
                - :class:`str` - homologous sequence at the second breakpoint

        Raises:
            AttributeError: for non specific breakpoints
        """

        b1_refseq = REFERENCE_GENOME[self.break1.chr].seq
        b2_refseq = REFERENCE_GENOME[self.break2.chr].seq
        if len(self.break1) > 1 or len(self.break2) > 1:
            raise AttributeError('cannot call shared sequence for non-specific breakpoints')

        # find for the first breakpoint
        p1 = self.break1.start - 1
        p2 = self.break2.start - 1
        shift1 = -1 if self.break1.orient == ORIENT.LEFT else 1
        shift2 = -1 if self.break2.orient == ORIENT.RIGHT else 1
        p2 += shift2

        first_seq = []
        while all([
            p1 >= 0,
            p1 < len(b1_refseq),
            p2 >= 0,
            p2 < len(b2_refseq),
            self.interchromosomal or p2 > self.break1.start - 1,
            self.interchromosomal or p1 < self.break2.start - 1
        ]):
            b1 = b1_refseq[p1] if self.break1.strand == STRAND.POS else reverse_complement(b1_refseq[p1])
            b2 = b2_refseq[p2] if self.break2.strand == STRAND.POS else reverse_complement(b2_refseq[p2])
            if DNA_ALPHABET.match(b1, b2):
                first_seq.append(b1_refseq[p1])
            else:
                break
            p1 += shift1
            p2 += shift2

        if shift1 < 0:
            first_seq.reverse()
        # now go over the second breakpoint
        p1 = self.break1.start - 1  # reset the start points
        p2 = self.break2.start - 1
        shift1 *= -1  # flip the directions
        shift2 *= -1
        p1 += shift1

        second_seq = []
        while all([
            p1 >= 0,
            p1 < len(b1_refseq),
            p2 >= 0,
            p2 < len(b2_refseq),
            self.interchromosomal or p2 > self.break1.start - 1,
            self.interchromosomal or p1 < self.break2.start - 1
        ]):
            b1 = b1_refseq[p1] if self.break1.strand == STRAND.POS else reverse_complement(b1_refseq[p1])
            b2 = b2_refseq[p2] if self.break2.strand == STRAND.POS else reverse_complement(b2_refseq[p2])
            if DNA_ALPHABET.match(b1, b2):
                second_seq.append(b2_refseq[p2])
            else:
                break
            p1 += shift1
            p2 += shift2

        if shift1 < 0:
            second_seq.reverse()

        return ''.join(first_seq).upper(), ''.join(second_seq).upper()

    def get_bed_repesentation(self):
        bed = []
        if self.interchromosomal:
            bed.append((self.break1.chr, self.break1.start, self.break1.end, self.data.get(COLUMNS.cluster_id, None)))
            bed.append((self.break2.chr, self.break2.start, self.break2.end, self.data.get(COLUMNS.cluster_id, None)))
        else:
            bed.append((self.break1.chr, self.break1.start, self.break2.end, self.data.get(COLUMNS.cluster_id, None)))
        return bed


def read_bpp_from_input_file(filename, expand_ns=True, explicit_strand=False, **kwargs):
    """
    reads a file using the TSV module. Each row is converted to a breakpoint pair and
    other column data is stored in the data attribute

    Args:
        filename (str): path to the input file
        expand_ns (bool): expand not specified orient/strand settings to all specific version
            (for strand this is only applied if the bam itself is stranded)
        explicit_strand (bool): used to stop unstranded breakpoint pairs from losing input strand information
    Returns:
        :class:`list` of :any:`BreakpointPair`: a list of pairs

    Example:
        >>> read_bpp_from_input_file('filename')
        [BreakpointPair(), BreakpointPair(), ...]

    One can also validate other expected columns that will go in the data attribute using the usual arguments
    to the TSV.read_file function

    .. code-block:: python

        >>> read_bpp_from_input_file('filename', cast={'index': int})
        [BreakpointPair(), BreakpointPair(), ...]
    """

    def soft_null_cast(value):
        try:
            value = TSV.null(value)
        except TypeError:
            pass
        return value

    def soft_boolean_cast(value):
        try:
            return TSV.tsv_boolean(value)
        except TypeError:
            if value == '?':
                value = 'null'
        return TSV.null(value)

    kwargs.setdefault('cast', {}).update(
        {
            COLUMNS.break1_position_start: int,
            COLUMNS.break1_position_end: int,
            COLUMNS.break2_position_start: int,
            COLUMNS.break2_position_end: int,
            COLUMNS.opposing_strands: soft_boolean_cast,
            COLUMNS.stranded: TSV.tsv_boolean,
            COLUMNS.untemplated_seq: soft_null_cast,
            COLUMNS.break1_chromosome: lambda x: re.sub('^chr', '', x),
            COLUMNS.break2_chromosome: lambda x: re.sub('^chr', '', x)
        })
    kwargs.setdefault('require', [])
    kwargs.setdefault('add', {}).update({
        COLUMNS.untemplated_seq: None,
        COLUMNS.break1_orientation: ORIENT.NS,
        COLUMNS.break1_strand: STRAND.NS,
        COLUMNS.break2_orientation: ORIENT.NS,
        COLUMNS.break2_strand: STRAND.NS,
        COLUMNS.opposing_strands: None
    })
    kwargs.setdefault('in_', {}).update(
        {
            COLUMNS.break1_orientation: ORIENT,
            COLUMNS.break1_strand: STRAND,
            COLUMNS.break2_orientation: ORIENT,
            COLUMNS.break2_strand: STRAND
        })
    header, rows = TSV.read_file(
        filename, suppress_index=True,
        **kwargs
    )
    restricted = [
        COLUMNS.break1_chromosome,
        COLUMNS.break1_position_start,
        COLUMNS.break1_position_end,
        COLUMNS.break1_strand,
        COLUMNS.break1_orientation,
        COLUMNS.break2_chromosome,
        COLUMNS.break2_position_start,
        COLUMNS.break2_position_end,
        COLUMNS.break2_strand,
        COLUMNS.break2_orientation,
        COLUMNS.stranded,
        COLUMNS.opposing_strands,
        COLUMNS.untemplated_seq
    ]
    pairs = []
    for row in rows:
        for attr, val in row.items():
            row[attr] = soft_null_cast(val)
        for attr in row:
            if attr in [COLUMNS.cluster_id, COLUMNS.annotation_id, COLUMNS.validation_id]:
                if not re.match('^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$', row[attr]):
                    raise AssertionError(
                        'error in column', attr, 'All mavis pipeline step ids must satisfy the regex:',
                        '^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$', row[attr])
        stranded = row[COLUMNS.stranded]
        opp = row[COLUMNS.opposing_strands]

        strand1 = row[COLUMNS.break1_strand] if (stranded or explicit_strand) else STRAND.NS
        strand2 = row[COLUMNS.break2_strand] if (stranded or explicit_strand) else STRAND.NS

        if explicit_strand and not expand_ns and {strand1, strand2} & {STRAND.NS}:
            raise AssertionError('cannot use explicit strand and not expand unknowns unless the strand is given')

        temp = []
        expand_strand = (stranded or explicit_strand) and expand_ns
        event_type = [None]
        if row.get(COLUMNS.event_type, None) not in [None, 'None']:
            try:
                event_type = row[COLUMNS.event_type].split(';')
                for et in event_type:
                    SVTYPE.enforce(event_type)
            except KeyError:
                pass

        for o1, o2, opp, s1, s2, et in itertools.product(
            ORIENT.expand(row[COLUMNS.break1_orientation]) if expand_ns else [row[COLUMNS.break1_orientation]],
            ORIENT.expand(row[COLUMNS.break2_orientation]) if expand_ns else [row[COLUMNS.break2_orientation]],
            [True, False] if opp is None and expand_ns else [opp],
            STRAND.expand(strand1) if expand_strand else [strand1],
            STRAND.expand(strand2) if expand_strand else [strand2],
            event_type
        ):
            try:
                b1 = Breakpoint(
                    row[COLUMNS.break1_chromosome],
                    row[COLUMNS.break1_position_start],
                    row[COLUMNS.break1_position_end],
                    strand=s1,
                    orient=o1
                )
                b2 = Breakpoint(
                    row[COLUMNS.break2_chromosome],
                    row[COLUMNS.break2_position_start],
                    row[COLUMNS.break2_position_end],
                    strand=s2,
                    orient=o2
                )

                data = {k: v for k, v in row.items() if k not in restricted}
                bpp = BreakpointPair(
                    b1,
                    b2,
                    opposing_strands=opp,
                    untemplated_seq=row[COLUMNS.untemplated_seq],
                    stranded=row[COLUMNS.stranded],
                )
                bpp.data.update(data)
                if et is not None:
                    bpp.data[COLUMNS.event_type] = et
                    if et not in BreakpointPair.classify(bpp):
                        raise InvalidRearrangement(
                            'error: expected one of', BreakpointPair.classify(bpp), 'but found', et, str(bpp), row)
                temp.append(bpp)
            except InvalidRearrangement as err:
                if not expand_ns:
                    raise err
            except AssertionError as err:
                if not expand_ns and not explicit_strand:
                    raise err
        if len(temp) == 0:
            raise InvalidRearrangement('could not produce a valid rearrangement', row)
        else:
            pairs.extend(temp)
    return pairs
