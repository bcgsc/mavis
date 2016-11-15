from __future__ import division

from copy import copy as sys_copy
from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.interval import Interval


class Breakpoint(Interval):
    """
    class for storing information about a SV breakpoint
    coordinates are given as 1-indexed
    """
    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand)

    def __init__(self, chr, start, end=None, strand=STRAND.NS, orient=ORIENT.NS, label=None):
        """
        Args:
            chr (str): the chromosome
            start (int): the genomic position of the breakpoint
            end (int, optional): if the breakpoint is uncertain (a range) then specify the end of the range here
            strand (STRAND, default=STRAND.NS): the strand
            orient (ORIENT, default=ORIENT.NS): the orientation (which side is retained at the break)
            label (str): the label for the breakpoint
        """
        Interval.__init__(self, start, end)
        self.orient = ORIENT.enforce(orient)
        self.chr = str(chr)
        self.strand = STRAND.enforce(strand)
        self.label = label

    def __repr__(self):
        temp = '{0}:{1}{2}{3}{4}'.format(
            self.chr, self.start, '-' + str(self.end) if self.end != self.start else '', self.orient, self.strand)
        if self.label is not None:
            temp += '#{0}'.format(self.label)
        return 'Breakpoint(' + temp + ')'

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)


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
        raise IndexError(
            'index input accessor is out of bounds: 1 or 2 only', index)

    @property
    def key(self):
        return self.break1.key, self.break2.key, self.opposing_strands, self.stranded, self.untemplated_sequence

    @property
    def interchromosomal(self):
        if self.break1.chr == self.break2.chr:
            return False
        return True

    def copy(self):
        temp = sys_copy(self)
        temp.break1 = sys_copy(self.break1)
        temp.break2 = sys_copy(self.break2)
        temp.flags = sys_copy(self.flags)
        return temp

    def __init__(self, b1, b2, stranded=False, opposing_strands=None, untemplated_sequence=None, flags=[], label=None):
        """
        Args:
            b1 (Breakpoint): the first breakpoint
            b2 (Breakpoint): the second breakpoint
            stranded (bool, default=False): if not stranded then +/- is equivalent to -/+
            opposing_strands (bool, optional): are the strands at the breakpoint opposite? i.e. +/- instead of +/+
            untemplated_sequence (str, optional): sequence between the breakpoints that is not part of either breakpoint
            flags (list, default=[]):
            label (str):
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
        self.flags = flags
        self.label = label

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
            raise AttributeError('must specify if opposing_strands')

        # try classifying to make sure it's a valid combination
        BreakpointPair.classify(self)

    def __str__(self):
        return 'BPP<{}==>{}; opposing={} seq={}>'.format(
            str(self.break1),
            str(self.break2),
            self.opposing_strands,
            repr(self.untemplated_sequence)
        )

    @classmethod
    def classify(cls, pair):
        """
        uses the chr, orientations and strands to determine the
        possible structural_variant types that this pair could support

        Args:
            pair (BreakpointPair): the pair to classify
        Returns:
            List[SVTYPE]: a list of possible SVTYPEs
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
        if read2 is None:
            return cls._call_from_single_contig(read1)
        return cls._call_from_paired_contig(read1, read2)

    @classmethod
    def _call_from_single_contig(cls, read):
        # first find the major unaligned event
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

        biggest = max(read_events, key=lambda x: x[2])

        first_breakpoint = read.reference_start - 1
        second_breakpoint = first_breakpoint
        untemplated_seq = ''

        seq_pos = 0
        for i, t in enumerate(read.cigar):
            v, f = t
            if i < biggest[0]:
                if v in [CIGAR.I, CIGAR.S]:
                    seq_pos += f
                elif v in [CIGAR.D, CIGAR.N]:
                    first_breakpoint += f
                    second_breakpoint += f
                else:
                    first_breakpoint += f
                    second_breakpoint += f
                    seq_pos += f
            elif i >= biggest[0] and i <= biggest[1]:
                if v in [CIGAR.I, CIGAR.S]:
                    if v == CIGAR.I:
                        untemplated_seq += read.query_sequence[seq_pos:seq_pos + f]
                    seq_pos += f
                elif v in [CIGAR.D, CIGAR.N]:
                    second_breakpoint += f
                else:
                    second_breakpoint += f
                    seq_pos += f
            else:
                break

        break1 = Breakpoint(
            read.reference_name,
            first_breakpoint,
            orient=ORIENT.LEFT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS
        )
        break2 = Breakpoint(
            read.reference_name,
            second_breakpoint + 1,  # need to add to be past the event
            orient=ORIENT.RIGHT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS
        )

        return BreakpointPair(break1, break2, opposing_strands=False, untemplated_sequence=untemplated_seq)

    @classmethod
    def _call_from_paired_contig(cls, read1, read2):
        """
        if the first read and the second read have overlapping range on their mutual query sequence
        then the breakpoint is called as is on the first and shifted on the second
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
        print('qci1', r1_qci)
        print('qci2', r2_qci)

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
            b1 = Breakpoint(read1.reference_name, read1.reference_end - 1, strand=s1, orient=o1)
        else:
            b1 = Breakpoint(read1.reference_name, read1.reference_start, strand=s1, orient=o1)

        overlap = 0
        if Interval.overlaps(r1_qci, r2_qci):
            # adjust the second read to remove the overlapping query region
            overlap = len(r1_qci & r2_qci)

        if o2 == ORIENT.RIGHT:
            b2 = Breakpoint(read2.reference_name, read2.reference_start + overlap, strand=s2, orient=o2)
        else:
            b2 = Breakpoint(read2.reference_name, read2.reference_end - 1 - overlap, strand=s2, orient=o2)

        # now check for untemplated sequence
        untemplated_sequence = ''
        d = Interval.dist(r1_qci, r2_qci)
        if d > 0:  # query coverage for read1 is after query coverage for read2
            untemplated_sequence = read1.query_sequence[r2_qci[1] + 1:r1_qci[0]]
        elif d < 0:  # query coverage for read2 is after query coverage for read1
            untemplated_sequence = read1.query_sequence[r1_qci[1] + 1:r2_qci[0]]
        else:  # query coverage overlaps
            pass

        if read1.is_reverse and untemplated_sequence != '':
            untemplated_sequence = reverse_complement(untemplated_sequence)

        return BreakpointPair(b1, b2, untemplated_sequence=untemplated_sequence)
