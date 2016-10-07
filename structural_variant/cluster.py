from __future__ import division

from structural_variant.constants import *
from structural_variant.interval import Interval
from structural_variant.breakpoint import BreakpointPair, Breakpoint

import scipy.stats as stat
import itertools
import numpy as np
import networkx as nx
import warnings

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
    
    for node in G.nodes():
        found = False
        for clique in refined_cliques:
            if node in clique:
                found = True
                break
        if not found:
            raise AssertionError('error, lost a node somehow', node, refined_cliques)
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
                    Breakpoint(chr1, start, start, orient = or1, strand = st1),
                    Breakpoint(chr2, end, end, orient = or2, strand = st2)
                    )
            clusters[cluster_breakpoint_pair] = current_cluster_breakpoint_pairs
        for node, freq in participants.items():
            if freq > 1:
                warnings.warn('node ' + node + ' participates in ' + str(freq) + ' clusters')
    
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


