from __future__ import division
from copy import copy as sys_copy

from .constants import CIGAR, COLUMNS, DNA_ALPHABET, ORIENT, reverse_complement, STRAND, SVTYPE
from .error import InvalidRearrangement, NotSpecifiedError
from .interval import Interval


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
        return {
            'chr': self.chr,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'seq': self.seq,
            'orientation': self.orient,
            'type': self.__class__.__name__
        }


class BreakpointPair:
    """
    """

    def __getattr__(self, attr):
        data = object.__getattribute__(self, 'data')
        try:
            return data[COLUMNS[attr]]
        except (KeyError, AttributeError):
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
        temp.data = {}
        temp.data.update(self.data)
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

        if self.break1.orient != ORIENT.NS and self.break2.orient != ORIENT.NS:
            if self.opposing_strands is not None:
                if (self.break1.orient == self.break2.orient and not self.opposing_strands) \
                        or (self.break1.orient != self.break2.orient and self.opposing_strands):
                    raise InvalidRearrangement(
                        'invalid breakpoint pair cannot form a valid combination', b1, b2, self.opposing_strands)
            else:
                self.opposing_strands = self.break1.orient == self.break2.orient
        # between break1 and break2 not in either
        self.untemplated_seq = untemplated_seq
        self.data = {}
        if data is not None:
            self.data.update(data)
            conflicts = set(data.keys()) & set(kwargs.keys())
            if conflicts:
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
        can be written directly as a tab row
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
        for col in temp:
            temp[col] = str(temp[col])
        row.update(temp)
        return row

    @classmethod
    def classify(cls, pair, discriminate=False):
        """
        uses the chr, orientations and strands to determine the
        possible structural_variant types that this pair could support

        Args:
            pair (BreakpointPair): the pair to classify
            discriminate (bool): flag if true will use untemplated sequence to discriminate between insertion/deletion classifications
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
                    if pair.untemplated_seq is not None and discriminate:
                        min_del_dist = max(pair.break2.start - pair.break1.end - 1, 0)
                        max_del_dist = pair.break2.end - pair.break1.start - 1
                        if len(pair.untemplated_seq) > max_del_dist:
                            return [SVTYPE.INS]
                        elif len(pair.untemplated_seq) < min_del_dist:
                            return [SVTYPE.DEL]
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
    def call_breakpoint_pair(cls, read1, read2=None, reference_genome=None):
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
            return cls._call_from_single_contig(read1, reference_genome=reference_genome)
        return cls._call_from_paired_contig(read1, read2, reference_genome=reference_genome)

    @classmethod
    def _call_from_single_contig(cls, read, reference_genome=None):
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
        for i, (state, freq) in enumerate(read.cigar):
            if state not in [CIGAR.I, CIGAR.D, CIGAR.N]:
                continue
            if i > 0 and read.cigar[i - 1][0] in [CIGAR.I, CIGAR.D, CIGAR.N]:
                start, _, size = read_events[-1]
                read_events[-1] = (start, i, size + freq)
            else:
                read_events.append((i, i, freq))

        if not read_events:
            raise UserWarning(
                'Cannot call event breakpoints from single contig. Contig does not contain any ins/del events')
        biggest = max(read_events, key=lambda x: x[2])

        first_breakpoint = read.reference_start - 1
        second_breakpoint = first_breakpoint
        untemplated_seq = ''

        seq_pos = 0
        seq_first_stop = -1
        seq_second_start = -1

        for i, (state, freq) in enumerate(read.cigar):
            if i < biggest[0]:  # before the event
                if state in [CIGAR.I, CIGAR.S]:
                    seq_pos += freq
                elif state in [CIGAR.D, CIGAR.N]:
                    first_breakpoint += freq
                    second_breakpoint += freq
                else:
                    first_breakpoint += freq
                    second_breakpoint += freq
                    seq_pos += freq
                seq_first_stop = seq_pos
            elif i >= biggest[0] and i <= biggest[1]:  # inside the event
                if state in [CIGAR.I, CIGAR.S]:
                    if state == CIGAR.I:
                        untemplated_seq += read.query_sequence[seq_pos:seq_pos + freq]
                    seq_pos += freq
                elif state in [CIGAR.D, CIGAR.N]:
                    second_breakpoint += freq
                else:
                    second_breakpoint += freq
                    seq_pos += freq
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
        if first_breakpoint == second_breakpoint and reference_genome:  # insertion or duplication
            refseq = reference_genome[read.reference_name].seq[
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
    def _call_from_paired_contig(cls, read1, read2, reference_genome=None):
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
            assert read1.query_sequence == read2.query_sequence
        else:
            assert read1.query_sequence == reverse_complement(read2.query_sequence)
            length = len(read1.query_sequence) - 1
            r2_qci = Interval(length - r2_qci.end, length - r2_qci.start)
        break1 = None
        break2 = None

        orient1 = ORIENT.NS
        orient2 = ORIENT.NS
        strand1 = STRAND.NEG if read1.is_reverse else STRAND.POS
        strand2 = STRAND.NEG if read2.is_reverse else STRAND.POS

        # <==============
        r1_st = read1.cigar[0][1] if read1.cigar[0][0] == CIGAR.S else 0
        r1_end = read1.cigar[-1][1] if read1.cigar[-1][0] == CIGAR.S else 0

        if r1_st > r1_end:
            orient1 = ORIENT.RIGHT
        elif r1_end > r1_st:
            orient1 = ORIENT.LEFT

        r2_st = read2.cigar[0][1] if read2.cigar[0][0] == CIGAR.S else 0
        r2_end = read2.cigar[-1][1] if read2.cigar[-1][0] == CIGAR.S else 0

        if r2_st > r2_end:
            orient2 = ORIENT.RIGHT
        elif r2_end > r2_st:
            orient2 = ORIENT.LEFT

        if orient1 == ORIENT.NS or orient2 == ORIENT.NS:
            raise AssertionError(
                'read does not have softclipping on either end and cannot therefore determine orientation',
                read1.cigar, read2.cigar)

        if orient1 == ORIENT.LEFT:  # ========++++
            break1 = Breakpoint(
                read1.reference_name,
                read1.reference_end,  # 1-based from 0-based exclusive end so don't need to -1
                strand=strand1,
                orient=orient1,
                seq=read1.query_alignment_sequence
            )
        else:
            break1 = Breakpoint(
                read1.reference_name,
                read1.reference_start + 1,  # 1-based from 0-based
                strand=strand1,
                orient=orient1,
                seq=read1.query_alignment_sequence
            )

        overlap = 0
        if Interval.overlaps(r1_qci, r2_qci):
            # adjust the second read to remove the overlapping query region
            overlap = len(r1_qci & r2_qci)

        if orient2 == ORIENT.RIGHT:
            break2 = Breakpoint(
                read2.reference_name,
                read2.reference_start + overlap + 1,  # 1-based from 0-based
                strand=strand2,
                orient=orient2,
                seq=read2.query_alignment_sequence[overlap:]
            )
        else:
            seq = read2.query_alignment_sequence
            if overlap > 0:
                seq = read2.query_alignment_sequence[:-1 * overlap]
            break2 = Breakpoint(
                read2.reference_name,
                read2.reference_end - overlap,  # 1-based from 0-based exclusive end so don't need to -1
                strand=strand2,
                orient=orient2,
                seq=seq
            )
        # now check for untemplated sequence
        untemplated_seq = ''
        dist = Interval.dist(r1_qci, r2_qci)
        if dist > 0:  # query coverage for read1 is after query coverage for read2
            untemplated_seq = read1.query_sequence[r2_qci[1] + 1:r1_qci[0]]
        elif dist < 0:  # query coverage for read2 is after query coverage for read1
            untemplated_seq = read1.query_sequence[r1_qci[1] + 1:r2_qci[0]]
        else:  # query coverage overlaps
            pass

        return BreakpointPair(break1, break2, untemplated_seq=untemplated_seq)

    def breakpoint_sequence_homology(self, reference_genome):
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
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by template/chr name

        Returns:
            tuple:
                - :class:`str` - homologous sequence at the first breakpoint
                - :class:`str` - homologous sequence at the second breakpoint

        Raises:
            AttributeError: for non specific breakpoints
        """

        b1_refseq = reference_genome[self.break1.chr].seq
        b2_refseq = reference_genome[self.break2.chr].seq
        if len(self.break1) > 1 or len(self.break2) > 1:
            raise AttributeError('cannot call shared sequence for non-specific breakpoints')

        # find for the first breakpoint
        pos1 = self.break1.start - 1
        pos2 = self.break2.start - 1
        shift1 = -1 if self.break1.orient == ORIENT.LEFT else 1
        shift2 = -1 if self.break2.orient == ORIENT.RIGHT else 1
        pos2 += shift2

        first_seq = []
        while all([
            pos1 >= 0,
            pos1 < len(b1_refseq),
            pos2 >= 0,
            pos2 < len(b2_refseq),
            self.interchromosomal or pos2 > self.break1.start - 1,
            self.interchromosomal or pos1 < self.break2.start - 1
        ]):
            break1_bp = b1_refseq[pos1] if self.break1.strand == STRAND.POS else reverse_complement(b1_refseq[pos1])
            break2_bp = b2_refseq[pos2] if self.break2.strand == STRAND.POS else reverse_complement(b2_refseq[pos2])
            if DNA_ALPHABET.match(break1_bp, break2_bp):
                first_seq.append(b1_refseq[pos1])
            else:
                break
            pos1 += shift1
            pos2 += shift2

        if shift1 < 0:
            first_seq.reverse()
        # now go over the second breakpoint
        pos1 = self.break1.start - 1  # reset the start points
        pos2 = self.break2.start - 1
        shift1 *= -1  # flip the directions
        shift2 *= -1
        pos1 += shift1

        second_seq = []
        while all([
            pos1 >= 0,
            pos1 < len(b1_refseq),
            pos2 >= 0,
            pos2 < len(b2_refseq),
            self.interchromosomal or pos2 > self.break1.start - 1,
            self.interchromosomal or pos1 < self.break2.start - 1
        ]):
            break1_bp = b1_refseq[pos1] if self.break1.strand == STRAND.POS else reverse_complement(b1_refseq[pos1])
            break2_bp = b2_refseq[pos2] if self.break2.strand == STRAND.POS else reverse_complement(b2_refseq[pos2])
            if DNA_ALPHABET.match(break1_bp, break2_bp):
                second_seq.append(b2_refseq[pos2])
            else:
                break
            pos1 += shift1
            pos2 += shift2

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
