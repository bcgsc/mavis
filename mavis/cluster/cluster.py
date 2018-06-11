from __future__ import division

from collections import namedtuple
from copy import copy
import itertools

from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import ORIENT, STRAND
from ..interval import Interval
from ..util import LOG, DEVNULL


class BreakpointPairGroupKey(namedtuple('BreakpointPairGroupKey', [
    'chr1', 'chr2', 'orient1', 'orient2', 'strand1', 'strand2', 'opposing_strands', 'explicit_strand'
])):

    def __new__(cls, chr1, chr2, orient1, orient2, strand1, strand2, opposing_strands=None, explicit_strand=False):
        if STRAND.NS in [strand1, strand2] and explicit_strand:
            raise ValueError('cannot have unspecified strand when explicit_strand is set')
        if not explicit_strand and opposing_strands is None:
            raise ValueError('opposing_strands must be specified when explicit_strand is false')
        if explicit_strand:
            opp = (strand1 != strand2)
            if opposing_strands is None:
                opposing_strands = opp
            elif opposing_strands != opp:
                raise ValueError('strand1 v strand2 v opposing_strands conflict.', strand1, strand2, opposing_strands)
        STRAND.enforce(strand1)
        STRAND.enforce(strand2)
        ORIENT.enforce(orient1)
        ORIENT.enforce(orient2)
        self = super(BreakpointPairGroupKey, cls).__new__(
            cls, chr1, chr2, orient1, orient2, strand1, strand2, opposing_strands, explicit_strand)
        return self


def weighted_mean(values, weights=None):
    if weights is None:
        weights = [1 for v in values]
    return sum(x * w for x, w in zip(values, weights)) / sum(weights)


def merge_integer_intervals(*intervals, weight_adjustment=0):
    """
    Merges a set of integer intervals into a single interval where the center is the
    weighted mean of the input intervals. The weight is inversely proportional to the
    length of each interval. The length of the final interval is the average of the lengths
    of the input intervals capped in size so that it never extends beyond the union of the
    input intervals

    Args:
        weight_adjustment (int): add to length to lower weighting differences between small intervals
    """
    float_offset = 0.99999999
    intervals = list(intervals)
    centers = []
    weights = []
    lengths = []

    if not intervals:
        raise AttributeError('cannot compute the weighted mean interval of an empty set of intervals')
    for i in range(0, len(intervals)):
        curr = intervals[i]
        intervals[i] = Interval(curr[0], curr[1] + float_offset)
        for _ in range(0, intervals[i].freq):
            centers.append(intervals[i].center)
            weights.append((weight_adjustment + 1) / (intervals[i].length() + weight_adjustment))
            lengths.append(intervals[i].length())

    center = round(weighted_mean(centers, weights=weights) * 2, 0) / 2
    size = weighted_mean(lengths)  # -1 b/c center counts as one
    start = max([round(center - size / 2, 0), min([i[0] for i in intervals])])
    end = min([round(center + size / 2, 0), max([i[1] for i in intervals])])
    offset = min([center - start, end - center])
    result = Interval(int(round(center - offset, 0)), int(round(center + max(0, offset - float_offset), 0)))
    return result


def pair_key(pair):
    return (
        pair.break1.chr,
        pair.break2.chr,
        pair.break1.start,
        pair.break2.start,
        pair.break1.end,
        pair.break2.end,
        pair.break1.orient,
        pair.break2.orient,
        pair.break1.strand if pair.stranded else STRAND.NS,
        pair.break2.strand if pair.stranded else STRAND.NS,
        pair.stranded,
        pair.opposing_strands)


def all_pair_group_keys(pair, explicit_strand=False):
    opt = [
        [pair.break1.chr],
        [pair.break2.chr],
        ORIENT.expand(pair.break1.orient),
        ORIENT.expand(pair.break2.orient),
        [STRAND.NS] if not explicit_strand else STRAND.expand(pair.break1.strand),
        [STRAND.NS] if not explicit_strand else STRAND.expand(pair.break2.strand),
        [pair.opposing_strands]
    ]
    result = []
    for c1, c2, o1, o2, s1, s2, opp in list(itertools.product(*opt)):
        if explicit_strand and (s1 != s2) != opp:
            continue
        elif opp == (o1 != o2):
            continue
        result.append(BreakpointPairGroupKey(c1, c2, o1, o2, s1, s2, opp, explicit_strand=explicit_strand))
    return result


def merge_by_union(input_pairs, group_key, weight_adjustment=10, cluster_radius=200):
    """
    for a given set of breakpoint pairs, merge the union of all pairs that are
    within the given distance (cluster_radius)
    """
    pairs_by_start = sorted(input_pairs, key=lambda x: x.break1.start)
    pairs_by_end = sorted(input_pairs, key=lambda x: x.break2.start)
    edges = {pair_key(p): set() for p in input_pairs}
    pairs_by_key = {}

    for i in range(0, len(input_pairs)):
        pairs_by_key.setdefault(pair_key(pairs_by_start[i]), []).append(pairs_by_start[i])
        for ordering in [pairs_by_start, pairs_by_end]:
            # try all combinations until start distance alone is too far
            curr = ordering[i]
            ckey = pair_key(curr)
            edges.setdefault(ckey, set())

            for j in range(i + 1, len(input_pairs)):
                other = ordering[j]
                okey = pair_key(other)
                distance = abs(Interval.dist(curr.break1, other.break1))
                if distance > cluster_radius:
                    break
                distance += abs(Interval.dist(curr.break2, other.break2))
                if distance <= cluster_radius:
                    edges[ckey].add(okey)
                    edges[okey].add(ckey)

    merged = set()
    merge_nodes = []
    for node in edges:
        if node in merged:
            continue
        adj = edges[node] | {node}
        merged.add(node)
        unmerged = adj - merged
        # follow edges to merge all connected nodes until all edges have been visited
        # extracts the current connected component
        while unmerged:
            for other in unmerged:
                adj.update(edges[other])
                merged.add(other)
            unmerged = adj - merged
        merge_nodes.append(adj)
    nodes = {}
    for node_keys in merge_nodes:
        pairs = []
        for pkey in node_keys:
            pairs.extend(pairs_by_key[pkey])
        itvl1 = merge_integer_intervals(*[p.break1 for p in pairs], weight_adjustment=weight_adjustment)
        itvl2 = merge_integer_intervals(*[p.break2 for p in pairs], weight_adjustment=weight_adjustment)
        if group_key.chr1 == group_key.chr2:
            itvl1.end = min(itvl2.end, itvl1.end)
            itvl2.start = max(itvl2.start, itvl1.start)
            itvl1.start = min(itvl1.start, itvl1.end)
            itvl2.end = max(itvl2.end, itvl2.start)
        b1 = Breakpoint(group_key.chr1, itvl1.start, itvl1.end, orient=group_key.orient1, strand=group_key.strand1)
        b2 = Breakpoint(group_key.chr2, itvl2.start, itvl2.end, orient=group_key.orient2, strand=group_key.strand2)
        # create the new bpp representing the merge of the input pairs
        new_bpp = BreakpointPair(
            b1, b2, opposing_strands=group_key.opposing_strands, stranded=group_key.explicit_strand)
        nodes.setdefault(new_bpp, []).extend(pairs)
    return nodes


def merge_breakpoint_pairs(input_pairs, cluster_radius=200, cluster_initial_size_limit=25, verbose=False):
    """
    two-step merging process

    1. merges all 'small' (see cluster_initial_size_limit) events as the union of all events that
        fall within the cluster_radius
    2. for all remaining events choose the 'best' merge for any event within cluster_radius of an
        existing node. Otherwise the node is added unmerged. The events in the second phase are
        done in order of smallest total breakpoint interval size to largest

    Args:
        input_pairs (list of BreakpointPair): the pairs to be merged
        cluster_radius (int) maximum distance allowed for a node to merge
        cluster_initial_size_limit (int): maximum size of breakpoint intervals allowed in the first merging phase

    Returns:
        dict of list of BreakpointPair by BreakpointPair: mapping of merged breakpoint pairs to the input pairs used in the merge
    """
    def pair_center_distance(pair1, pair2):
        d = abs(pair1.break1.center - pair2.break1.center)
        d += abs(pair1.break2.center - pair2.break2.center)
        return d
    mapping = {}
    groups = {}  # split the groups by putative pairings
    pair_weight = {}
    explicit_strand = False
    phase2_groups = {}
    for pair in input_pairs:
        if pair.stranded:
            explicit_strand = True
            break

    doubled = 0
    for i, old_pair in enumerate(input_pairs):
        pair = copy(old_pair)
        pair.data['tag'] = i
        k = pair_key(pair)
        pair_weight.setdefault(k, []).append(pair)

        putative_group_keys = all_pair_group_keys(pair, explicit_strand=explicit_strand)
        doubled += len(putative_group_keys)
        if len(putative_group_keys) < 1:
            raise NotImplementedError('bad breakpoint input does not fit any groups', pair)
        for key in putative_group_keys:
            if len(pair.break1) + len(pair.break2) > cluster_initial_size_limit:
                phase2_groups.setdefault(key, []).append(pair)
            else:
                groups.setdefault(key, []).append(pair)
    # now try all pairwise combinations within groups
    for group_key in sorted(set(list(groups) + list(phase2_groups))):
        count = len(groups.get(group_key, [])) + len(phase2_groups.get(group_key, []))
        if verbose:
            LOG(group_key, 'pairs:', count)
        nodes = merge_by_union(
            groups.get(group_key, []), group_key,
            weight_adjustment=cluster_initial_size_limit, cluster_radius=cluster_radius)

        # phase 2. Sort all the breakpoint pairs left by size and merge the smaller ones in first
        # this is be/c we assume that a larger breakpoint interval indicates less certainty in the call
        phase2_pairs = sorted(
            phase2_groups.get(group_key, []), key=lambda p: (len(p.break1) + len(p.break2), pair_key(p)))

        for pair in phase2_pairs:
            distances = sorted([(pair_center_distance(pair, node), node) for node in nodes], key=lambda x: x[0])
            merged = False

            if len(distances) > 0:
                best = min(distances, key=lambda x: x[0])
                for dist, node in distances:
                    if dist > best[0] or dist > cluster_radius:
                        break
                    pairs = nodes[node] + [pair]

                    itvl1 = merge_integer_intervals(
                        *[p.break1 for p in pairs], weight_adjustment=cluster_initial_size_limit)
                    itvl2 = merge_integer_intervals(
                        *[p.break2 for p in pairs], weight_adjustment=cluster_initial_size_limit)
                    if group_key.chr1 == group_key.chr2:
                        itvl1.end = min(itvl2.end, itvl1.end)
                        itvl2.start = max(itvl2.start, itvl1.start)
                        itvl1.start = min(itvl1.start, itvl1.end)
                        itvl2.end = max(itvl2.end, itvl2.start)

                    b1 = Breakpoint(
                        group_key.chr1, itvl1.start, itvl1.end, orient=group_key.orient1, strand=group_key.strand1)
                    b2 = Breakpoint(
                        group_key.chr2, itvl2.start, itvl2.end, orient=group_key.orient2, strand=group_key.strand2)

                    new_bpp = BreakpointPair(
                        b1, b2, opposing_strands=group_key.opposing_strands, stranded=explicit_strand)
                    del nodes[node]
                    nodes.setdefault(new_bpp, []).extend(pairs)
                    merged = True
            if not merged:
                b1 = Breakpoint(
                    group_key.chr1, pair.break1.start, pair.break1.end,
                    orient=group_key.orient1, strand=group_key.strand1)

                b2 = Breakpoint(
                    group_key.chr2, pair.break2.start, pair.break2.end,
                    orient=group_key.orient2, strand=group_key.strand2)

                new_bpp = BreakpointPair(
                    b1, b2, opposing_strands=group_key.opposing_strands, stranded=explicit_strand)
                nodes.setdefault(new_bpp, []).append(pair)
        if verbose:
            LOG('merged', count, 'down to', len(nodes))
        for node, pairs in nodes.items():
            if node in mapping:
                raise KeyError('duplicate merge node', str(node), node, pair_key(node))
            mapping[node] = pairs
    # assertion to check that no nodes were left out of merging
    merge_sources = set()
    for merge_node, sources in mapping.items():
        merge_sources.update([p.data['tag'] for p in sources])
    if len(merge_sources) != len(input_pairs):
        raise AssertionError('merged node inputs ({}) does not equal the number of pairs input ({})'.format(
            len(merge_sources), len(input_pairs)))
    return mapping
