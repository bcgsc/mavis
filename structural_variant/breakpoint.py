from __future__ import division

from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.interval import Interval
import warnings

class Breakpoint:
    """
    class for storing information about a SV breakpoint
    coordinates are given as 1-indexed
    """
    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand, self.label, self.left_seq, self.right_seq)

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

    def __init__(self, chr, interval_start, interval_end=None, **kwargs):
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
        return hash(self.key)
 
class BreakpointPair:    
    @property
    def key(self):
        return self.break1.key, self.break2.key, self.opposing_strands, self.stranded, self.untemplated_sequence

    def __init__(self, b1, b2, **kwargs):
        if b1.key > b2.key:
            self.break1 = b2
            self.break2 = b1
        else:
            self.break1 = b1
            self.break2 = b2
        self.stranded = kwargs.pop('stranded', False)
        self.opposing_strands = kwargs.pop('opposing_strands', None)
        self.untemplated_sequence = kwargs.pop('untemplated_sequence', None) # between break1 and break2 not in either
        self.classification = kwargs.pop('classification', None)
        
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
                    raise InvalidRearrangement('invalid breakpoint pair cannot form a valid combination', b1, b2, self.opposing_strands)

    def __repr__(self):
        return str(self)
    
    def __str__(self):
        return '{0}==>{1}{2}'.format(
                str(self.break1), 
                str(self.break2), 
                ( '[OPP]' if self.opposing_strands else '[EQ]') if self.opposing_strands is not None else '')
    
    @classmethod
    def classify(cls, pair):
        """
        uses the chr, orientations and strands to determine the
        possible structural_variant types that this pair could support
        """
        if pair.break1.chr == pair.break2.chr: # intrachromosomal
            if pair.opposing_strands is None:
                if pair.break1.orient == ORIENT.LEFT:    
                    if pair.break2.orient == ORIENT.LEFT: # LL
                        return [SVTYPE.INV]
                    elif pair.break2.orient == ORIENT.RIGHT: # LR
                        return [SVTYPE.DEL, SVTYPE.INS]
                    else: # L?
                        return [SVTYPE.DEL, SVTYPE.INS, SVTYPE.INV]
                elif pair.break1.orient == ORIENT.RIGHT: 
                    if pair.break2.orient == ORIENT.LEFT: # RL
                        return [SVTYPE.DUP]
                    elif pair.break2.orient == ORIENT.RIGHT: # RR
                        return [SVTYPE.INV]
                    else: # R?
                        return [SVTYPE.DUP, SVTYPE.INV]
                else:
                    if pair.break2.orient == ORIENT.LEFT: #?L
                        return [SVTYPE.DUP, SVTYPE.INV]
                    elif pair.break2.orient == ORIENT.RIGHT: #?R
                        return [SVTYPE.INV, SVTYPE.DEL, SVTYPE.INS]
                    else: # ??
                        return [SVTYPE.DUP, SVTYPE.INV, SVTYPE.DEL, SVTYPE.INS]
            elif pair.opposing_strands:
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
        else: # interchromosomal
            if pair.opposing_strands is None: # can't tell if they're opposite strands or not
                if pair.break1.orient not in [ORIENT.LEFT, ORIENT.RIGHT] \
                        or pair.break2.orient not in [ORIENT.LEFT, ORIENT.RIGHT]:
                            return [SVTYPE.ITRANS, SVTYPE.TRANS]
                elif pair.break1.orient == pair.break2.orient:
                    return [SVTYPE.ITRANS]
                else:
                    return [SVTYPE.TRANS]
            elif pair.opposing_strands:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.RIGHT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.LEFT):
                            raise InvalidRearrangement(pair)
                return [SVTYPE.ITRANS]
            else:
                if (pair.break1.orient == ORIENT.LEFT and pair.break2.orient == ORIENT.LEFT) \
                        or (pair.break1.orient == ORIENT.RIGHT and pair.break2.orient == ORIENT.RIGHT):
                            raise InvalidRearrangement(pair)
                return [SVTYPE.TRANS] 

class SVAnnotation:
    def __init__(self, breakpoint_pair, event_type, transcript1, transcript2, **kwargs):
        """
        holds the association between a pair of transcripts and an event
        transcript1 and transcript2 can both be None or neither can be None
        """
        self.breakpoint_pair = breakpoint_pair
        self.sv_class = SVTYPE.enforce(event_type)

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
