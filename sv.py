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

    def __init__(self, start, end = None, **kwargs):
        self.start = int( start )
        self.end = int( end ) if end is not None else self.start
        if self.start > self.end:
            raise AttributeError('interval start > end is not allowed')
        self.freq = int(kwargs.pop('freq', 1))
        if self.freq <= 0:
            raise AttributeError('Interval frequency must be a natural number')
    
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
        temp = str(self.start)
        if self.end != self.start:
            temp += '-' + str(self.end)
        if self.freq != 1:
            temp += 'x' + str(self.freq)
        return self.__class__.__name__ + '(' + temp + ')'

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
    def weighted_mean(cls, intervals):
        """
        returns the weighted mean for a set of intervals
        the weight is the inverse of the size of the interval
        so that broader intervals are weighted less than
        more specific/tighter intervals
        """
        if len(intervals) == 0:
            raise AttributeError('input list cannot be empty')
        first = next(iter(intervals))
        centers = []
        weights = []
        
        for i in intervals:
            for temp in range(0, i.freq):
                centers.append(i.center)
                weights.append(1 / len(i))

        return np.average(centers,weights=weights)

    def combine(self, other):
        """
        adding two intervals returns the minimum interval that covers both input intervals
        >>> x, y, z = ( Interval(1, 4), Interval(-1, 0), Interval(1, 2) )
        >>> x.combine(y)
        Interval(-1-4)
        >>> x.combine(z)
        Interval(1-4)
        >>> y.combine(z)
        Interval(-1-2)
        """
        return Interval(min(self.start, other.start), max(self.end, other.end))
    
    def __eq__(self, other):
        if not hasattr(other, 'start') \
                or not hasattr(other, 'end') \
                or not hasattr(other, 'weight') \
                or self.start != other.start \
                or self.end != other.end \
                or self.freq != other.freq:
            return False
        return True
    
    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start:
            if self.end < other.end:
                return True
            elif self.end == other.end:
                if self.freq < other.freq:
                    return True
        return False

    def __gt__(self, other):
        if self.start > other.start:
            return True
        elif self.start == other.start:
            if self.end > other.end:
                return True
            elif self.end == other.end:
                if self.freq > other.freq:
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
        return hash((self.start, self.end, self.freq))
    
    @classmethod
    def paired_weighted_means(cls, intervals):
        int_a = Interval.weighted_mean( [ x[0] for x in intervals ] )
        int_b = Interval.weighted_mean( [ x[1] for x in intervals ] )
        return int_a, int_b
    
    @classmethod
    def paired_set_distance(cls, intervals, other_intervals):
        """
        for two sets of interval pairs (as tuples) computes the weighted mean
        of each interval set (a total of four) and then returns the 
        distance between the sets of pairs as the cumulative distance
        of the weighted means of their pairs
        """
        int_a, int_b = cls.paired_weighted_means(intervals)
        oint_a, oint_b = cls.paired_weighted_means(other_intervals)
        return abs(int_a - oint_a) + abs(int_b - oint_b)

    @classmethod
    def redundant_ordered_hierarchical_clustering(cls, clusters, **kwargs):
        """
        for an input set of of clusters, do hierarchical clustering
        redundant b/c we allow clusters to be grouped more than once
        into either of their immediate neighbours
        """
        r = int(kwargs.pop('r'))
        if kwargs:
            raise AttributeError('invalid parameter', kwargs)
        if r < 0:
            raise AttributeError('r must be a positive integer')
        # order the clusters by weighted mean
        complete = []
        queue = sorted(clusters, key=lambda x: cls.paired_weighted_means(x) )
    
        
        while len(queue) > 0:
            temp_queue = []
            
            for i in range(0, len(queue)):
                curr = queue[i]
                joined = False
                
                if i > 0:
                    dist = cls.paired_set_distance(curr, clusters[i - 1])
                    if dist <= r:
                        joined = True
                if i < len(queue) - 1:
                    dist = cls.paired_set_distance(curr, clusters[i + 1])
                    if dist <= r:
                        temp_queue.append(curr.union(clusters[i + 1]))
                        joined = True 
                if not joined:
                    complete.append(curr)
            queue = temp_queue
        return complete

    @classmethod
    def union(cls, intervals):
        """
        returns the union of the set of input intervals
        
        >>> l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        >>> Interval.union(l)
        Interval(1-10)
        >>> l.append(Interval(11))
        >>> Interval.union(l)
        Interval(1-11)
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the union of an empty set of intervals')
        curr = next(iter(intervals))
        low = curr.start
        high = curr.end
        
        for i in intervals:
            if i.start < low:
                low = i.start
            if i.end > high:
                high = i.end
        return Interval(low, high)

    @classmethod
    def intersection(cls, intervals):
        """
        returns None if there is no intersection
        
        >>> l = [Interval(1, 10), Interval(5, 7), Interval(7)]
        >>> Interval.intersection(l)
        Interval(7)
        >>> l.append(Interval(11))
        >>> Interval.intersection(l)
        """
        if len(intervals) < 1:
            raise AttributeError('cannot compute the intersection of an empty set of intervals')
        curr = next(iter(intervals))
        low = curr.start
        high = curr.end
        
        for i in intervals:
            if high < i.start or i.end < low:
                return None
            low = max(low, i.start)
            high = min(high, i.end)
        return Interval(low, high)

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
    def freq(self):
        return self.pos.freq
    
    @classmethod
    def weighted_mean(cls, breakpoints):
        if len(breakpoints) == 0:
            raise AttributeError('cannot calculate the weighted mean of an empty list')
        return Interval.weighted_mean([ b.pos for b in breakpoints ])

    def __init__(self, chr, interval_start, interval_end, orient, strand, **kwargs):
        self.orient = ORIENT.enforce( orient )
        self.chr = str(chr)
        self.pos = Interval(interval_start, interval_end)
        self.strand = STRAND.enforce( strand )
        self.label = kwargs.pop('label', None)
    
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
        #self.gtype = kwargs.pop('gtype', 'DNA')
        self.stranded = kwargs.pop('stranded', False)
        self.opposing_strands = kwargs.pop('opposing_strands', None)
        if self.break1.strand != STRAND.NS and self.break2.strand != STRAND.NS:
            opposing = self.break1.strand != self.break2.strand
            if self.opposing_strands is None:
                self.opposing_strands = opposing
            elif self.opposing_strands != opposing:
                raise AttributeError('conflict in input arguments, opposing_strands must agree with input breakpoints'
                    'when the strand has been specified')
    
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

        >>> a = Breakpoint(1, 50, 51, 'R', '+')
        >>> b = Breakpoint(1, 10, 11, 'L', '-')
        >>> c = Breakpoint('X', 50, 51, 'R', '+')
        >>> d = Breakpoint(1, 10, 11, 'R', '-')
        >>> e = Breakpoint(1, 10, 11, 'L', '+')
        >>> f = Breakpoint(1, 50, 51, 'R', '-')
        >>> g = Breakpoint(1, 50, 51, 'R', '?')
        >>> h = Breakpoint(1, 10, 11, 'L', '?')

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


def redundant_maximal_kcliques(G, **kwargs):
    """
    for a give graph returns all cliques up to a size k
    any clique which is a proper subset of another clique is removed
    nodes can participate in multiple cliques if they are equal fit
    """
    k = int( kwargs.pop('k') )
    if kwargs:
        raise AttributeError('invalid parameter', kwargs)
    if k < 1:
        raise AttributeError('k must be greater than 0')
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
    
    # remove subsets to ensure cliques are maximal (up to k)
    refined_cliques = []
    for i in range(0, len(cliques)):
        is_subset = False
        for j in range(i + 1, len(cliques)):
            if cliques[i].issubset(cliques[j]):
                is_subset = True
                break
        if not is_subset:
            refined_cliques.append(cliques[i])
    
    participation = {}
    for c in refined_cliques:
        for node in c:
            participation[node] = participation.get(node, 0) + 1

    for node, count in sorted(participation.items(), reverse=True):
        distances = []
        for cluster in refined_cliques:
            if node not in cluster:
                continue
            d = Interval.paired_set_distance([node], cluster)
            distances.append((d, cluster))
        lowest = min(distances, key = lambda x: x[0])[0]
        for score, cluster in distances:
            if score > lowest:
                cluster.remove(node)

    return refined_cliques


def cluster_breakpoints(input_pairs, **kwargs):
    # 0. sort the breakpoints by start and then end
    # 1a. split/duplicate breakpoints into sets of things that could possibly support the same event
    # 1b. split breakpoint pairs by chr pair (can be the same chr)
    # 2. set the initial clusters based on overlap
    # 3. iterate over the clusters 
    #   # stop when no clusters improve/change or we hit a maximum number of iterations
    r = int( kwargs.pop('r') )
    k = int( kwargs.pop('k') )
    if kwargs:
        raise AttributeError('invalid parameter', kwargs)
    # classify the breakpoints.... by the possible pairs they could support (explicit only)
    classify = {}
    track_breakpoints = {}

    for pair in input_pairs:
        
        interval_pair = ( Interval(pair.break1.start, pair.break1.end), Interval(pair.break2.start, pair.break2.end) )
        interval_key = (interval_pair[0].start, interval_pair[0].end, interval_pair[1].start, interval_pair[1].end)
        
        for chr1, chr2, o1, o2, s1, s2 in itertools.product(
                [pair.break1.chr], 
                [pair.break2.chr],
                [pair.break1.orient] if pair.break1.orient != ORIENT.NS else [ORIENT.LEFT, ORIENT.RIGHT],
                [pair.break2.orient] if pair.break2.orient != ORIENT.NS else [ORIENT.LEFT, ORIENT.RIGHT],
                [pair.break1.strand] if pair.break1.strand != STRAND.NS else [STRAND.POS, STRAND.NEG],
                [pair.break2.strand] if pair.break2.strand != STRAND.NS else [STRAND.POS, STRAND.NEG]
                ):
            classification_key = chr1, chr2, o1, o2, s1, s2

            classification_keys = [classification_key]

            if not pair.stranded:
                skey = (chr1, chr2, o1, o2, 
                        STRAND.NEG if s1 == STRAND.POS else STRAND.POS,
                        STRAND.NEG if s2 == STRAND.POS else STRAND.POS)
                classification_keys.append(skey)
            
            for ckey in classification_keys:
                if ckey not in classify:
                    classify[ckey] = {}
                    track_breakpoints[ckey] = set()
                track_breakpoints[ckey].add(pair)
                
                if interval_key not in classify[ckey]:
                    classify[ckey][interval_key] = interval_pair
                else:
                    classify[ckey][interval_key][0].freq += 1
                    classify[ckey][interval_key][1].freq += 1
    
    # set the initial clusters
    
    # nodes are the breakpoint pair keys
    # edges are the distance between matched breakpoints (less than 2k)
    clusters = {}
    events = 0
    for key, group in sorted(classify.items()):
        chr1, chr2, or1, or2, st1, st2 =  key
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
            if dist <= 2 * r:
                G.add_edge(curr, other)
        # get all the connected components
        # for all connected components find all cliques of a given size k
        complete_subgraphs = redundant_maximal_kcliques(G, k=k) # every node will be present
        # before we group the subgraphs, should choose the best fit for each node in the subgraphs
        new_sets = Interval.redundant_ordered_hierarchical_clustering(complete_subgraphs, r=r)
        # now using the complete subgraphs as our inputs, do ordered hierarchical clustering
        participants = {}
        for clique in new_sets:
            # create the new weighted interval and then add to the clusters set
            for node in clique:
                participants[node] = participants.get(node, 0) + 1
            events += 1
            size = len(clique)
            # add the cluster breakpoint as the key with a list of the breakpoints used to create it
            current_cluster_breakpoint_pairs = set()
            for f, s in clique:
                for pair in track_breakpoints[key]:
                    if pair.break1.start == f.start \
                            and pair.break1.end == f.end \
                            and pair.break2.start == s.start \
                            and pair.break2.end == s.end:
                        current_cluster_breakpoint_pairs.add(pair)

            start, end = Interval.paired_weighted_means(clique)
            cluster_breakpoint_pair = BreakpointPair(
                    Breakpoint(chr1, start, start, or1, st1),
                    Breakpoint(chr2, end, end, or2, st2)
                    )
            clusters[cluster_breakpoint_pair] = current_cluster_breakpoint_pairs
        for node, freq in participants.items():
            if freq > 1:
                print('node', node, 'participates in', freq, 'clusters')
    print('found', events, 'events')
    
    # now for each cluster, merge clusters with the same location and equivalent strand information
    merged_clusters = {}
    for centroid, cluster_set in clusters.items():
        grouping_key = (
                centroid.opposing_strands, # should always be determined for centroids
                centroid.break1.chr, 
                centroid.break2.chr, 
                centroid.break1.orient, 
                centroid.break2.orient,
                centroid.break1.start,
                centroid.break1.end,
                centroid.break2.start,
                centroid.break2.end)
        # anything else with this same grouping key can differ only in strands and should
        # be grouped assuming there are no stranded pairs in either group
        merged = False
        has_stranded = True if sum([1 if pair.stranded else 0 for pair in cluster_set]) > 0 else False
        for cluster, cset in merged_clusters.get(grouping_key, []):
            curr_has_stranded = True if sum([1 if pair.stranded else 0 for pair in cset]) > 0 else False
            if not has_stranded and not curr_has_stranded:
                # can merge when one strand has not been specified
                cset.update(cluster_set)
                cluster.break1.strand = STRAND.NS
                cluster.break2.strand = STRAND.NS
                merged = True
        if not merged:
            if grouping_key not in merged_clusters:
                merged_clusters[grouping_key] = []
            merged_clusters[grouping_key].append((centroid, cluster_set))

    final_clusters = {}
    for group in merged_clusters:
        for centroid, cset in merged_clusters[group]:
            final_clusters[centroid] = cset

    return final_clusters

if __name__ == '__main__':
    import doctest
    doctest.testmod()
