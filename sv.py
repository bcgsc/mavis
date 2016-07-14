from __future__ import division

from vocab import Vocab
from constants import ORIENT
from constants import STRAND

import scipy.stats as stat
import itertools
import numpy as np
import networkx as nx
import warnings

class Interval:
    #TYPE = Vocab(INC='inclusive', EXC='exclusive')

    def __init__(self, start, end = None, **kwargs):
        self.start = int( start )
        #self.start_type = Interval.TYPE.enforce( kwargs.pop('start_type', Interval.TYPE.INC) )
        self.end = int( end ) if end is not None else self.start
        #self.end_type = Interval.TYPE.enforce( kwargs.pop('end_type', Interval.TYPE.INC) )
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed')
        #if self.start == self.end and self.start_type == Interval.TYPE.EXC and self.end_type != Interval.TYPE.EXC
        self.weight = int(kwargs.pop('weight', 1))
        if self.weight < 0:
            raise AttributeError('Intervals do not support negative weights')
    
    def overlap(self, other):
        """
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(0, 3) )
        >>> x.overlap(y)
        False
        >>> x.overlap(z)
        True
        >>> y.overlap(x)
        False
        >>> y.overlap(z)
        True
        """
        if self - other == 0:
            return True
        return False
    
    def __len__(self):
        return self.end - self.start + 1

    def __repr__(self):
        if self.weight != 1:
            return str(self.start) + '_' + str(self.end) + 'x' + str(self.weight)
        return str(self.start) + '_' + str(self.end)
    
    @property
    def center(self):
        """
        >>> x, y, z = ( Interval(1, 5), Interval(-1, 0), Interval(1, 2) )
        >>> x.center
        3.0
        >>> y.center
        -0.5
        >>> z.center
        1.5
        """
        return self.start + (len(self) - 1)/2
    
    @classmethod
    def weighted_mean(cls, *intervals):
        if len(intervals) == 0:
            raise AttributeError('input list cannot be empty')
        first = intervals[0]
        centers = []
        weights = []
        
        for i in intervals:
            for temp in range(0, i.weight):
                centers.append(i.center)
                weights.append(1 / len(i))

        return np.average(centers,weights=weights)

    def combine(self, other):
        """
        adding two intervals returns the minimum interval that covers both input intervals
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(1, 2) )
        >>> x.combine(y)
        [-1, 4]
        >>> x.combine(z)
        [1, 4]
        >>> y.combine(z)
        [-1, 2]
        """
        return Interval(min(self.start, other.start), max(self.end, other.end))
    
    def __eq__(self, other):
        if not hasattr(other, 'start') \
                or not hasattr(other, 'end') \
                or not hasattr(other, 'weight') \
                or self.start != other.start \
                or self.end != other.end \
                or self.weight != other.weight:
            return False
        return True
    
    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start:
            if self.end < other.end:
                return True
            elif self.end == other.end:
                if self.weight < other.weight:
                    return True
        return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start:
            if self.end > other.end:
                return True
            elif self.end == other.end:
                if self.weight > other.weight:
                    return True
        return False

    def __sub__(self, other):
        """ 
        returns the minimum distance between intervals
        
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(0, 3) )
        >>> x - y
        1
        >>> y - x
        -1
        >>> x - z
        0
        >>> z - x
        0
        """
        if self.end < other.start:
            return self.end - other.start
        elif self.start > other.end:
            return self.start - other.end
        else:
            return 0
    
    def __hash__(self):
        return hash((self.start, self.end, self.weight))

class Breakpoint:
    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand)

    @property
    def start(self):
        return self.pos.start
    
    @property
    def end(self):
        return self.pos.end
    
    @property
    def weight(self):
        return self.pos.weight
    
    @classmethod
    def weighted_mean(cls, *breakpoints):
        if len(breakpoints) == 0:
            raise AttributeError('cannot calculate the weighted mean of an empty list')
        return Interval.weighted_mean(*[ b.pos for b in breakpoints ])

    def __init__(self, chr, interval_start, interval_end, orient, strand, **kwargs):
        self.orient = ORIENT.enforce( orient )
        self.chr = chr
        self.pos = Interval(interval_start, interval_end)
        self.strand = STRAND.enforce( strand )
        self.label = kwargs.pop('label', None)
    
    def __repr__(self):
        temp = '{0}:{1}_{2}{3}{4}'.format(self.chr, self.start, self.end, self.orient, self. strand)
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
    
    def expand(self):
        breakpoints = []
        
        strand_options = [self.strand]
        if self.strand == STRAND.NS:
            strand_options = [STRAND.POS, STRAND.NEG]
        
        orient_options = [self.orient]
        if self.orient == ORIENT.NS:
            orient_options = [ORIENT.LEFT, ORIENT.RIGHT]
        
        for s in strand_options:
            for o in orient_options:
                breakpoints.append( 
                        Breakpoint(
                            self.chr, 
                            self.pos.start, 
                            self.pos.end,
                            o,
                            s
                            )
                        )
        return breakpoints

def is_complete(G, N):
    """
    for a given input graph and a set of nodes N in G
    checks if N is a complete subgraph of G
    """
    for node, other in itertools.combinations(N, 2):
        if not G.has_node(node) or not G.has_node(other):
            raise AttributeError('invalid node is not part of the input graph')
        if not G.has_edge(node, other):
            return False
    return True 


def generate_redundant_maximal_kcliques(G, **kwargs):
    """
    for a give graph returns all cliques up to a size k
    any clique which is a proper subset of another clique is removed
    nodes can participate in multiple cliques
    """
#print('generate_redundant_maximal_kcliques N(G) =', len(G.nodes()))
    k = int( kwargs.pop('k', 10) )
    if k >= 20:
        warnings.warn('k >= 20 is not recommended as the number of combinations increases exponentially')
    
    cliques = []
    for component in nx.connected_components(G):
        comp_cliques = []
        # take an exhaustive approach to finding the possible cliques
        for ktemp in range(1, k + 1):
            for putative_kclique in itertools.combinations(component, ktemp):
                if is_complete(G, putative_kclique):
                    cliques.append(set(putative_kclique))
#print('generate_redundant_maximal_kcliques initial cliques:', len(cliques))
    # remove subsets to ensure cliques are maximal (up to k)
    refined_cliques = []
    for i in range(0, len(cliques)):
        is_subset = False
        for j in range(i + 1, len(cliques)):
            if cliques[i] <= cliques[j]:
                is_subset = True
                break
        if not is_subset:
            refined_cliques.append(cliques[i])
#print('generate_redundant_maximal_kcliques final =', len(refined_cliques) )
    return refined_cliques

def cluster_breakpoints(input_pairs, **kwargs):
    # 0. sort the breakpoints by start and then end
    # 1a. split/duplicate breakpoints into sets of things that could possibly support the same event
    # 1b. split breakpoint pairs by chr pair (can be the same chr)
    # 2. set the initial clusters based on overlap
    # 3. iterate over the clusters 
    #   # stop when no clusters improve/change or we hit a maximum number of iterations
    r = int( kwargs.pop('r', 10) )
    k = int( kwargs.pop('k', 10) )
    
    # classify the breakpoints.... by the possible pairs they could support (explicit only)
    classify = {}

    for pair in input_pairs:
        b1i, b2i = pair
        if b1i.key > b2i.key:
            b1i, b2i = (b2i, b1i)
        
        interval_pair = ( Interval(b1i.start, b1i.end), Interval(b2i.start, b2i.end) )
        key = (interval_pair[0].start, interval_pair[0].end, interval_pair[1].start, interval_pair[1].end)

        for b1, b2 in itertools.product(b1i.expand(), b2i.expand()):
            explicit = (b1.chr, b2.chr, b1.orient, b2.orient, b1.strand, b2.strand)
            opposite = (
                    b1.chr, b2.chr, b1.orient, b2.orient,
                    STRAND.POS if b1.strand == STRAND.NEG else STRAND.NEG,
                    STRAND.POS if b2.strand == STRAND.NEG else STRAND.NEG
                    ) # only applicable to non-stranded inputs
            if explicit not in classify:
                classify[explicit] = {}
            if opposite not in classify:
                classify[opposite] = {}
            
            if key not in classify[explicit]:
                classify[explicit][key] = interval_pair
            else:
                classify[explicit][key][0].weight += 1
                classify[explicit][key][1].weight += 1
            if key not in classify[opposite]:
                classify[opposite][key] = interval_pair
            else:
                classify[opposite][key][0].weight += 1
                classify[opposite][key][1].weight += 1
    
    # set the initial clusters
    
    # nodes are the breakpoint pair keys
    # edges are the distance between matched breakpoints (less than 2k)
    events = 0
    for k, group in sorted(classify.items()):
        chr1, chr2, or1, or2, st1, st2 =  k
        #if chr1 == chr2:
        #    continue
        group = group.values()
        G = nx.Graph()
        # set up the edges for this grouping
        for n in group:
            G.add_node(n)
        for curr, other in itertools.combinations(group, 2):
            dist = abs(curr[0] - other[0]) \
                    + abs(curr[1] - other[1])
            if dist < 2 * r:
                G.add_edge(curr, other)
        # get all the connected components
        # for all connected components find all cliques of a given size k
        for clique in generate_redundant_maximal_kcliques(G):
            events += 1
            size = len(clique)
            if size < 2:
                continue
            start = int(round(Interval.weighted_mean(*[b[0] for b in clique]), 0))
            end = int(round(Interval.weighted_mean(*[b[1] for b in clique]), 0))
            WINDOW = 100
            print('clique:', size, '{0}:{1}-{2} {3}:{4}-{5}'.format(
                chr1, start - WINDOW, start + WINDOW, chr2, end - WINDOW, end + WINDOW))
    print('found', events, 'events')
    
    # cluster until no splitting is performed

    return classify

if __name__ == '__main__':
    import doctest
    doctest.testmod()
