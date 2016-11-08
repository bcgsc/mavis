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
        return 'BPP<{0}==>{1}{2}{3}>'.format(
            str(self.break1),
            str(self.break2),
            ('; OPP' if self.opposing_strands else '; EQ') if self.opposing_strands is not None else '',
            '; seq=\'{}\''.format(self.untemplated_sequence) if self.untemplated_sequence else ''
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
