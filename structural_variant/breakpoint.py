from __future__ import division

from vocab import Vocab
from structural_variant.constants import ORIENT
from structural_variant.constants import STRAND
from structural_variant.constants import SVTYPE
from structural_variant.interval import Interval

import scipy.stats as stat
import itertools
import numpy as np
import networkx as nx
import warnings

class Breakpoint:
    """
    class for storing information about a SV breakpoint
    coordinates are given as 1-indexed
    """
    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand, self.left_seq, self.right_seq)

    @property
    def start(self):
        return self.pos.start
    
    @property
    def end(self):
        return self.pos.end
    
    @property
    def freq(self):
        return self.pos.freq
    
    @classmethod
    def weighted_mean(cls, breakpoints):
        if len(breakpoints) == 0:
            raise AttributeError('cannot calculate the weighted mean of an empty list')
        return Interval.weighted_mean([ b.pos for b in breakpoints ])
    
    def softclipped_sequence(self):
        if self.orient == ORIENT.NS:
            raise UserWarning('invalid on non-specific orientations')
        elif self.orient == ORIENT.LEFT:
            return self.right_seq
        else:
            return self.left_seq

    def __init__(self, chr, interval_start, interval_end, **kwargs):
        self.orient = ORIENT.enforce( kwargs.pop('orient', ORIENT.NS) )
        self.chr = str(chr)
        self.pos = Interval(interval_start, interval_end)
        self.strand = STRAND.enforce( kwargs.pop('strand', STRAND.NS) )
        self.label = kwargs.pop('label', None)
        self.left_seq = kwargs.pop('left_seq', None)
        self.right_seq = kwargs.pop('right_seq', None)
    
    def __repr__(self):
        temp = '{0}:{1}{2}{3}{4}'.format(
                self.chr, self.start, '-' + str(self.end) if self.end != self.start else '', self.orient, self. strand)
        if self.label is not None:
            temp += '#{0}'.format(self.label)
        return 'Breakpoint(' + temp + ')'

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        if self.key != other.key:
            return False
        return True
    
    def __hash__(self):
        return hash((self.chr, self.pos, self.strand, self.orient, self.label))
 
class BreakpointPair:    
    @property
    def key(self):
        return self.break1.key, self.break2.key, self.opposing_strands, self.stranded

    def __init__(self, b1, b2, **kwargs):
        if b1.key > b2.key:
            self.break1 = b2
            self.break2 = b1
        else:
            self.break1 = b1
            self.break2 = b2
        self.stranded = kwargs.pop('stranded', False)
        self.opposing_strands = kwargs.pop('opposing_strands', None)
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
                    raise UserWarning('invalid breakpoint pair cannot form a valid combination', b1, b2, self.opposing_strands)

    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return '{0}==>{1}{2}'.format(
                str(self.break1), 
                str(self.break2), 
                ( '[OPP]' if self.opposing_strands else '[EQ]') if self.opposing_strands is not None else '')

    def comparable(self, other):
        """
        determines if two pairs are comparable i.e. if they could possibly support the same event

        >>> a = Breakpoint(1, 50, 51, orient = 'R', strand = '+')
        >>> b = Breakpoint(1, 10, 11, orient = 'L', strand = '-')
        >>> c = Breakpoint('X', 50, 51, orient = 'R', strand = '+')
        >>> d = Breakpoint(1, 10, 11, orient = 'R', strand = '-')
        >>> e = Breakpoint(1, 10, 11, orient = 'L', strand = '+')
        >>> f = Breakpoint(1, 50, 51, orient = 'R', strand = '-')
        >>> g = Breakpoint(1, 50, 51, orient = 'R', strand = '?')
        >>> h = Breakpoint(1, 10, 11, orient = 'L', strand = '?')

        # different chromosome
        >>> BreakpointPair(a, b, stranded = False).comparable(BreakpointPair(c, b, stranded = False))
        False

        # same pair
        >>> BreakpointPair(a, b, stranded = False).comparable(BreakpointPair(a, b, stranded = False))
        True

        # same pair and stranded
        >>> BreakpointPair(a, b, stranded = True).comparable(BreakpointPair(a, b, stranded = True))
        True

        # different orientation
        >>> BreakpointPair(a, b, stranded = False).comparable(BreakpointPair(a, d, stranded = False))
        False

        # different relative strand
        >>> BreakpointPair(a, b, stranded = False).comparable(BreakpointPair(a, e, stranded = False))
        False

        # same relative strands, but opposite strands and stranded not specified
        >>> BreakpointPair(a, b, stranded = False).comparable(BreakpointPair(f, e, stranded = False))
        True

        # same relative strands, but opposite strands and stranded specified
        >>> BreakpointPair(a, b, stranded = True).comparable(BreakpointPair(f, e, stranded = True))
        False

        # adding not specified strands
        >>> BreakpointPair(a, b, stranded = True).comparable(BreakpointPair(g, e, stranded = True))
        False

        # adding not specified strands
        >>> BreakpointPair(a, b, stranded = True).comparable(BreakpointPair(f, h, stranded = True))
        False

        # adding not specified strands
        >>> BreakpointPair(a, b, stranded = True).comparable(BreakpointPair(g, h, stranded = True))
        True
        """
        if self.break1.chr != other.break1.chr \
                or self.break2.chr != other.break2.chr \
                or (self.break1.orient != other.break1.orient 
                        and ORIENT.NS not in [self.break1.orient, other.break1.orient]) \
                or (self.break2.orient != other.break2.orient 
                        and ORIENT.NS not in [self.break2.orient, other.break2.orient]):
                    return False
        elif self.stranded and other.stranded: # both of them cares about the strand
            if (self.break1.strand == other.break1.strand 
                    or STRAND.NS in [self.break1.strand, other.break1.strand]) \
                    and (self.break2.strand == other.break2.strand or 
                            STRAND.NS in [self.break2.strand, other.break2.strand]):
                return True
            else:
                return False
        else:
            if STRAND.NS in [self.break1.strand, self.break2.strand, other.break1.strand, other.break2.strand]:
                return True
            elif (self.break1.strand == self.break2.strand) == (other.break1.strand == other.break2.strand):
                return True
            else:
                return False

class SVAnnotation:
    def __init__(self, breakpoint_pair, event_type, transcript1, transcript2, **kwargs):
        """
        holds the association between a pair of transcripts and an event
        transcript1 and transcript2 can both be None or neither can be None
        """
        self.breakpoint_pair = breakpoint_pair
        self.event_type = SVTYPE.enforce(event_type)
        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.HUMAN_REFERENCE_GENOME = kwargs.pop('HUMAN_REFERENCE_GENOME', None)
    
    def fusion_ref_sequence(self):
        pass

    def fusion_frame(self):
        pass


if __name__ == '__main__':
    import doctest
    doctest.testmod()
