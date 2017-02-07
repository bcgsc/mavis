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
            strand (STRAND): the strand
            orient (ORIENT): the orientation (which side is retained at the break)

        Examples:
            >>> Breakpoint('1', 1, 2)
            >>> Breakpoint('1', 1)
            >>> Breakpoint('1', 1, 2, '+', 'R')
            >>> Breakpoint('1', 1, orient='R')
        """
        Interval.__init__(self, start, end)
        self.orient = ORIENT.enforce(orient)
        self.chr = str(chr)
        self.strand = STRAND.enforce(strand)
        self.seq = seq

    def __repr__(self):
        temp = '{0}:{1}{2}{3}{4}'.format(
            self.chr, self.start, '-' + str(self.end) if self.end != self.start else '', self.orient, self.strand)
        return 'Breakpoint(' + temp + ')'

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
            'sequence': self.seq,
            'orientation': self.orient,
            'type': self.__class__.__name__
        }
        return d


class BreakpointPair:
    """
    """
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
        return hash((self.break1.key, self.break2.key, self.opposing_strands, self.stranded, self.untemplated_sequence))

    def __eq__(self, other):
        for attr in ['break1', 'break2', 'opposing_strands', 'stranded', 'untemplated_sequence']:
            if not hasattr(other, attr):
                return False
            elif getattr(self, attr) != getattr(other, attr):
                return False
        return True

    @property
    def interchromosomal(self):
        """(:class:`bool`): True if the breakpoints are on different chromosomes, False otherwise"""
        if self.break1.chr == self.break2.chr:
            return False
        return True

    def copy(self):
        temp = sys_copy(self)
        temp.break1 = sys_copy(self.break1)
        temp.break2 = sys_copy(self.break2)
        return temp

    def __init__(self, b1, b2, stranded=False, opposing_strands=None, untemplated_sequence=None, data={}):
        """
        Args:
            b1 (Breakpoint): the first breakpoint
            b2 (Breakpoint): the second breakpoint
            stranded (bool): if not stranded then +/- is equivalent to -/+
            opposing_strands (bool): are the strands at the breakpoint opposite? i.e. +/- instead of +/+
            untemplated_sequence (str): sequence between the breakpoints that is not part of either breakpoint
            data (dict): optional dictionary of attributes associated with this pair

        Note:
            untemplated_sequence should always be given wrt to the positive/forward reference strand

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
        self.untemplated_sequence = untemplated_sequence
        self.data = {}
        self.data.update(data)

        if self.break1.strand != STRAND.NS and self.break2.strand != STRAND.NS:
            opposing = self.break1.strand != self.break2.strand
            if self.opposing_strands is None:
                self.opposing_strands = opposing
            elif self.opposing_strands != opposing:
                raise AttributeError('conflict in input arguments, opposing_strands must agree with input breakpoints'
                                     'when the strand has been specified')
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
        return 'BPP<{}==>{}; opposing={} seq={}>'.format(
            str(self.break1),
            str(self.break2),
            self.opposing_strands,
            repr(self.untemplated_sequence)
        )

    def flatten(self):
        """
        returns the key-value self for the breakpoint self information as
        can be written directly as a TSV row
        """
        row = {}
        row.update(self.data)
        temp = {
            COLUMNS.break1_chromosome.name: self.break1.chr,
            COLUMNS.break1_position_start.name: self.break1.start,
            COLUMNS.break1_position_end.name: self.break1.end,
            COLUMNS.break1_orientation.name: self.break1.orient,
            COLUMNS.break1_strand.name: self.break1.strand,
            COLUMNS.break2_chromosome.name: self.break2.chr,
            COLUMNS.break2_position_start.name: self.break2.start,
            COLUMNS.break2_position_end.name: self.break2.end,
            COLUMNS.break2_orientation.name: self.break2.orient,
            COLUMNS.break2_strand.name: self.break2.strand,
            COLUMNS.opposing_strands.name: self.opposing_strands,
            COLUMNS.stranded.name: self.stranded,
            COLUMNS.untemplated_sequence.name: self.untemplated_sequence
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
            >>> BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
            [SVTYPE.INV]
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
    def call_breakpoint_pair(cls, read1, read2=None):
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
            return cls._call_from_single_contig(read1)
        return cls._call_from_paired_contig(read1, read2)

    @classmethod
    def _call_from_single_contig(cls, read):
        """
        calls a set of breakpoints from a pysam style read using the cigar values

        .. todo::

            return multiple events not just the major event

        Returns:
            BreakpointPair: the newly called breakpoint pair from the contig

        Raises:
            UserWarning: if the contig does not contain insertions/deletions
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
        return BreakpointPair(break1, break2, opposing_strands=False, untemplated_sequence=untemplated_seq)

    @classmethod
    def _call_from_paired_contig(cls, read1, read2):
        """
        calls a set of breakpoints from a pair of pysam style reads using their softclipping to
        find the breakpoints and orientations
        if the first read and the second read have overlapping range on their mutual query sequence
        then the breakpoint is called as is on the first and shifted on the second

        .. todo::
            return multiple events not just the major event
        """

        if read1.reference_id > read2.reference_id:
            read1, read2 = (read2, read1)
        elif read1.reference_id == read2.reference_id and read1.reference_start > read2.reference_start:
            read1, read2 = (read2, read1)

        r1_qci = read1.query_coverage_interval()
        r2_qci = read2.query_coverage_interval()

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
        if read1.cigar[0][0] == CIGAR.S and read1.cigar[-1][0] == CIGAR.S:
            if read1.cigar[0][1] > read1.cigar[-1][1]:
                o1 = ORIENT.RIGHT
            elif read1.cigar[0][1] < read1.cigar[-1][1]:
                o1 = ORIENT.LEFT
        elif read1.cigar[0][0] == CIGAR.S:
            o1 = ORIENT.RIGHT
        elif read1.cigar[-1][0] == CIGAR.S:
            o1 = ORIENT.LEFT

        if read2.cigar[0][0] == CIGAR.S and read2.cigar[-1][0] == CIGAR.S:
            if read2.cigar[0][1] > read2.cigar[-1][1]:
                o2 = ORIENT.RIGHT
            elif read2.cigar[0][1] < read2.cigar[-1][1]:
                o2 = ORIENT.LEFT
        elif read2.cigar[0][0] == CIGAR.S:
            o2 = ORIENT.RIGHT
        elif read2.cigar[-1][0] == CIGAR.S:
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
        untemplated_sequence = ''
        d = Interval.dist(r1_qci, r2_qci)
        if d > 0:  # query coverage for read1 is after query coverage for read2
            untemplated_sequence = read1.query_sequence[r2_qci[1] + 1:r1_qci[0]]
        elif d < 0:  # query coverage for read2 is after query coverage for read1
            untemplated_sequence = read1.query_sequence[r1_qci[1] + 1:r2_qci[0]]
        else:  # query coverage overlaps
            pass

        return BreakpointPair(b1, b2, untemplated_sequence=untemplated_sequence)

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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            :class:`tuple` of :class:`str` and :class:`str`: tuple contains
                - str: homologous sequence at the first breakpoint
                - str: homologous sequence at the second breakpoint

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


def read_bpp_from_input_file(filename, expand_ns=True, **kwargs):
    """
    reads a file using the TSV module. Each row is converted to a breakpoint pair and
    other column data is stored in the data attribute

    Args:
        filename (str): path to the input file
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
            COLUMNS.untemplated_sequence: soft_null_cast,
            COLUMNS.break1_chromosome: lambda x: re.sub('^chr', '', x),
            COLUMNS.break2_chromosome: lambda x: re.sub('^chr', '', x)
        })
    kwargs.setdefault('require', [])
    kwargs.setdefault('add', {}).update({
        COLUMNS.untemplated_sequence: None,
        COLUMNS.break1_orientation: ORIENT.NS,
        COLUMNS.break1_strand: STRAND.NS,
        COLUMNS.break2_orientation: ORIENT.NS,
        COLUMNS.break2_strand: STRAND.NS,
        COLUMNS.opposing_strands: None
    })
    kwargs.setdefault('_in').update(
        {
            COLUMNS.break1_orientation: ORIENT,
            COLUMNS.break1_strand: STRAND,
            COLUMNS.break2_orientation: ORIENT,
            COLUMNS.break2_strand: STRAND
        })
    header, rows = TSV.read_file(
        filename,
        **kwargs
    )
    pairs = []
    for row in rows:
        stranded = row[COLUMNS.stranded]
        opp = row[COLUMNS.opposing_strands]
        for o1, o2, opp, s1, s2 in itertools.product(
            ORIENT.expand(row[COLUMNS.break1_orientation]) if expand_ns else [row[COLUMNS.break1_orientation]],
            ORIENT.expand(row[COLUMNS.break2_orientation]) if expand_ns else [row[COLUMNS.break2_orientation]],
            [True, False] if opp is None and expand_ns else [opp],
            STRAND.expand(row[COLUMNS.break1_strand]) if stranded and expand_ns else [row[COLUMNS.break1_strand]],
            STRAND.expand(row[COLUMNS.break2_strand]) if stranded and expand_ns else [row[COLUMNS.break2_strand]]
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
                bpp = BreakpointPair(
                    b1,
                    b2,
                    opposing_strands=opp,
                    untemplated_sequence=row[COLUMNS.untemplated_sequence],
                    stranded=row[COLUMNS.stranded],
                    data=row
                )
                pairs.append(bpp)
            except (AttributeError, InvalidRearrangement) as err:
                if not expand_ns:
                    raise err
    return pairs
