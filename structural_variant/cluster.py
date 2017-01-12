from __future__ import division

from .constants import *
from .error import *
from .interval import Interval
from .breakpoint import BreakpointPair, Breakpoint

import itertools
import networkx as nx
import warnings


class IntervalPair:
    """
    """
    def __init__(self, start, end, **kwargs):
        """
        Args:
            start (Interval): the first interval
            end (Interval): the second interval
        """
        self.data = kwargs
        self.start = start if isinstance(start, Interval) else Interval(start[0], start[1])
        self.end = end if isinstance(end, Interval) else Interval(end[0], end[1])

    def __eq__(self, other):
        if not hasattr(other, 'start') or not hasattr(other, 'end') \
                or self.start != other.start or self.end != other.end:
            return False
        return True

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start and self.end < other.end:
            return True
        return False

    def __hash__(self):
        return hash((self.start, self.end))

    @classmethod
    def weighted_mean(cls, *interval_pairs):
        """
        returns a new IntervalPair where the start interval is the weighted mean of the starts of all
        the input interval pairs, similar for the end

        Args:
            interval_pairs (IntervalPair): interval pairs
        Returns:
            IntervalPair: the new IntervalPair
        """
        start = Interval.weighted_mean(*[i.start for i in interval_pairs])
        end = Interval.weighted_mean(*[i.end for i in interval_pairs])
        return IntervalPair(start, end)

    def dist(self, other):
        """
        computes the distance between IntervalPairs by averaging the distance between
        the first Interval centers of each and the second Interval centers of each

        Returns:
            int: the distance between interval pairs
        """
        d = abs(self.start.center - other.start.center)
        d += abs(self.end.center - other.end.center)
        d /= 2
        return d

    def __repr__(self):
        return '{}<{}, {}, data={}>'.format(self.__class__.__name__, self.start, self.end, self.data)

    @classmethod
    def _redundant_maximal_kcliques(cls, G, k=10):
        """
        for a give graph returns all cliques up to a size k
        any clique which is a proper subset of another clique is removed
        nodes can participate in multiple cliques if they are equal fit
        """
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

        for count, node in sorted([(c, n) for n, c in participation.items() if c > 1], reverse=True):
            distances = []
            for cluster in refined_cliques:
                if node not in cluster:
                    continue
                d = sum([node.dist(x) for x in cluster if x != node]) / (len(cluster) - 1)
                distances.append((d, cluster))
            lowest = min(distances, key=lambda x: x[0])[0]
            for score, cluster in distances:
                if score > lowest:
                    cluster.remove(node)

        for node in G.nodes():
            found = False
            for clique in refined_cliques:
                if node in clique:
                    found = True
                    break
            if not found:
                raise AssertionError(
                    'error, lost a node somehow', node, refined_cliques)
        return refined_cliques

    @classmethod
    def _redundant_ordered_hierarchical_clustering(cls, groups, r):
        """
        given a set of IntervalPair objects group sets that have less than
        a given distance between weighted means of the different groups

        Args:
            groups (list of set of IntervalPair): a list of sets of interval pairs
            r (int): the distance to determine grouping
        """
        queue = sorted(groups, key=lambda x: IntervalPair.weighted_mean(*x))
        complete_groups = []

        while len(queue) > 0:
            temp_queue = []
            for i in range(0, len(queue)):
                merged = False
                curr = queue[i]
                curri = IntervalPair.weighted_mean(*curr)
                if i > 0:
                    prev = queue[i - 1]
                    if IntervalPair.weighted_mean(*prev).dist(curri) <= r:
                        d = curr | prev
                        if d not in temp_queue:
                            temp_queue.append(d)
                        merged = True
                if i < len(queue) - 1:
                    nexxt = queue[i + 1]
                    if IntervalPair.weighted_mean(*nexxt).dist(curri) <= r:
                        d = curr | nexxt
                        if d not in temp_queue:
                            temp_queue.append(d)
                        merged = True
                if not merged:
                    complete_groups.append(curr)
            queue = sorted(temp_queue, key=lambda x: IntervalPair.weighted_mean(*x))
        return complete_groups

    @classmethod
    def cluster(cls, pairs, r, k):
        """
        clusters a list of IntervalPair objects

        Args:
            pairs (list of IntervalPair): list of IntervalPair objects
            r (int): the distance for grouping clusters
            k (int): the clique size to look for

        Returns:
            list of set of IntervalPair: a list of sets of interval pairs representing their clusters/groupings
        """
        # build the initial graph
        G = nx.Graph()
        for p in pairs:
            G.add_node(p)
        for curr, other in itertools.combinations(pairs, 2):
            if curr.dist(other) <= r:
                G.add_edge(curr, other)

        # pull out the highly connected components
        subgraphs = cls._redundant_maximal_kcliques(G, k)
        subgraphs = cls._redundant_ordered_hierarchical_clustering(subgraphs, r)
        return subgraphs


def is_complete(G, N):
    """
    for a given input graph and a set of nodes N in G
    checks if N is a complete subgraph of G

    Args:
        G (nx.Graph): the input supergraph
        N (list): a list of nodes in G
    Returns:
        bool: True if N as a subgraph of G is complete False otherwise
    """
    for node, other in itertools.combinations(N, 2):
        if not G.has_node(node) or not G.has_node(other):
            raise AttributeError('invalid node is not part of the input graph')
        if not G.has_edge(node, other):
            return False
    return True


def cluster_breakpoint_pairs(input_pairs, r, k):
    # 0. sort the breakpoints by start and then end
    # 1a. split/duplicate breakpoints into sets of things that could possibly support the same event
    # 1b. split breakpoint pairs by chr pair (can be the same chr)
    # 2. set the initial clusters based on overlap
    # 3. iterate over the clusters
    #   # stop when no clusters improve/change or we hit a maximum number of iterations
    # classify the breakpoints.... by the possible pairs they could support
    # (explicit only)
    input_pairs = list(input_pairs)
    node_sets = {}
    input_mapping = {}  # new node to input index

    for index, bpp in enumerate(input_pairs):
        added = False
        for o1, o2 in itertools.product(
                ORIENT.expand(bpp.break1.orient),
                ORIENT.expand(bpp.break2.orient),
        ):
            try:
                temp = BreakpointPair.copy(bpp)
                bpp.break1.orient = o1
                bpp.break2.orient = o2
                BreakpointPair.classify(bpp)  # will throw error if invalid combination
                b1 = Interval(bpp.break1.start, bpp.break1.end)
                b2 = Interval(bpp.break2.start, bpp.break2.end)
                new_bpp = IntervalPair(b1, b2)
                classification_key = (
                    bpp.break1.chr, bpp.break2.chr,
                    o1, o2,
                    bpp.break1.strand,
                    bpp.break2.strand,
                    bpp.opposing_strands,
                    bpp.stranded,
                    bpp.untemplated_sequence
                )
                node_sets.setdefault(classification_key, set()).add(new_bpp)
                input_mapping.setdefault(new_bpp, set()).add(index)
                added = True
            except InvalidRearrangement:
                pass
        if not added:
            raise AssertionError('error. did not add to clustering', bpp)

    result = {}
    for ckey, group in sorted(node_sets.items()):
        chr1, chr2, o1, o2, s1, s2, opposing_strands, stranded, seq = ckey
        clusters = IntervalPair.cluster(group, r, k)

        for node in group:
            particpation = sum([1 for c in clusters if node in c])
            if particpation > 1:
                warnings.warn('interval pair participates in multiple clusters')
            elif particpation == 0:
                raise AssertionError('error: dropped input pair did not complete clustering', node)
        for c in clusters:
            ip = IntervalPair.weighted_mean(*c)
            # create the new breakpoint pair that represents the cluster
            bpp = BreakpointPair(
                Breakpoint(chr1, ip.start[0], ip.start[1], strand=s1, orient=o1),
                Breakpoint(chr2, ip.end[0], ip.end[1], strand=s2, orient=o2),
                opposing_strands=opposing_strands,
                untemplated_sequence=seq,
                stranded=stranded
            )
            # gather the original input pairs using the mapping
            original_input_pairs = itertools.chain.from_iterable([input_mapping[node] for node in c])
            result.setdefault(bpp, set()).update(original_input_pairs)

    all_input_indices = set()
    for bpp, inputs in result.items():
        all_input_indices.update(inputs)

    for i in range(0, len(input_pairs)):
        if i not in all_input_indices:
            raise AssertionError('input breakpoint pair was not clustered', i, str(input_pairs[i]))

    for bpp in result:
        result[bpp] = [input_pairs[i] for i in result[bpp]]

    return result
