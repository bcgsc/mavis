import itertools
import warnings

import distance
import networkx as nx

from .bam import cigar as _cigar
from .bam.read import calculate_alignment_score, nsb_align, sequence_complexity
from .constants import reverse_complement
from .interval import Interval
from .util import DEVNULL


class Contig:
    """
    """

    def __init__(self, sequence, score):
        self.seq = sequence
        self.remapped_sequences = {}  # alignment score contribution on the contig by read
        self.score = score
        self.alignments = set()
        self.input_reads = set()
        self.strand_specific = False

    def __hash__(self):
        return hash(self.seq)

    def complexity(self):
        return sequence_complexity(self.seq)

    def add_mapped_sequence(self, read, multimap=1):
        self.remapped_sequences[read] = 1 / multimap

    def remap_score(self):
        return sum(self.remapped_sequences.values())

    def remap_coverage(self):
        itvls = Interval.min_nonoverlapping(*[(x.reference_start, x.reference_end - 1) for x in self.remapped_sequences])
        cov = sum([len(i) for i in itvls])
        return cov / len(self.seq)

    def remap_depth(self, query_range=None):
        """
        the average depth of remapped reads over a give range of the contig sequence

        Args:
            query_range (Interval): 1-based inclusive range
        """
        if query_range is None:
            query_range = Interval(1, len(self.seq))
        if query_range.start < 1 or query_range.end > len(self.seq):
            raise ValueError('query range must be within contig seq', query_range, len(self.seq))
        result = 0
        for read, weight in self.remapped_sequences.items():
            read_qrange = Interval(read.reference_start + 1, read.reference_end)
            if Interval.overlaps(query_range, read_qrange):
                result += len(query_range & read_qrange) * weight
        return result / len(query_range)


class DeBruijnGraph(nx.DiGraph):
    """
    wrapper for a basic digraph
    enforces edge weights
    """

    def get_edge_freq(self, n1, n2):
        """
        returns the freq from the data attribute for a specified edge
        """
        if not self.has_edge(n1, n2):
            raise KeyError('missing edge', n1, n2)
        data = self.get_edge_data(n1, n2)
        return data['freq']

    def add_edge(self, n1, n2, freq=1):
        """
        add a given edge to the graph, if it exists add the frequency to the existing frequency count
        """
        if self.has_edge(n1, n2):
            data = self.get_edge_data(n1, n2)
            freq += data['freq']
        nx.DiGraph.add_edge(self, n1, n2, freq=freq)

    def all_edges(self, *nodes, data=False):
        return self.in_edges(*nodes, data=data) + self.out_edges(*nodes, data=data)

    def trim_tails_by_freq(self, min_weight):
        """
        for any paths where all edges are lower than the minimum weight trim

        Args:
            min_weight (int): the minimum weight for an edge to be retained
        """
        ends = {n for n in self.nodes() if self.out_degree(n) == 0 or self.in_degree(n) == 0}
        visited = set()

        while ends:
            curr = ends.pop()
            if not self.has_node(curr) or curr in visited:
                continue
            visited.add(curr)
            # follow until the path forks or we run out of low weigh edges
            if self.out_degree(curr) == 0 or self.in_degree(curr) == 0:
                for src, tgt, data in list(self.all_edges(curr, data=True)):
                    if data['freq'] < min_weight:
                        self.remove_edge(src, tgt)
                    if src not in visited:
                        ends.add(src)
                    if tgt not in visited:
                        ends.add(tgt)

        # remove any resulting singlets
        for node in visited:
            if not self.has_node(node):
                continue
            if self.degree(node) == 0:
                self.remove_node(node)

    def trim_forks_by_freq(self, min_weight):
        """
        for all nodes in the graph, if the node has an out-degree > 1 and one of the outgoing
        edges has freq < min_weight. then that outgoing edge is deleted
        """
        nodes = [n for n in self.nodes() if self.degree(n) > 2]
        for node in nodes:
            if self.out_degree(node) > 1:
                outgoing_edges = self.out_edges(node, data=True)
                best = max([e[2]['freq'] for e in outgoing_edges])
                for src, tgt, data in outgoing_edges:
                    if data['freq'] < min_weight and data['freq'] != best:
                        self.remove_edge(src, tgt)
            if self.in_degree(node) > 1:
                ingoing_edges = self.in_edges(node, data=True)
                best = max([e[2]['freq'] for e in ingoing_edges])
                for src, tgt, data in ingoing_edges:
                    if data['freq'] < min_weight and data['freq'] != best:
                        self.remove_edge(src, tgt)

    def trim_noncutting_paths_by_freq(self, min_weight):
        """
        trim any low weight edges where another path exists between the source and target
        of higher weight
        """
        current_edges = list(self.all_edges(data=True))
        for src, tgt, data in sorted(current_edges, key=lambda x: (x[2]['freq'], x[0], x[1])):
            # come up with the path by extending this edge either direction until the degree > 2
            if not self.has_node(src) or not self.has_node(tgt) or not self.has_edge(src, tgt):
                continue

            if src == tgt and data['freq'] < min_weight:
                self.remove_edge(src, tgt)
            else:
                path = []
                while self.in_degree(src) == 1 and self.out_degree(src) == 1:
                    s, t, data = self.in_edges(src, data=True)[0]
                    if data['freq'] >= min_weight or s in path:
                        break
                    path.insert(0, src)
                    src = s
                path.insert(0, src)

                while self.in_degree(tgt) == 1 and self.out_degree(tgt) == 1:
                    s, t, data = self.out_edges(tgt, data=True)[0]
                    if data['freq'] >= min_weight or t in path:
                        break
                    path.append(tgt)
                    tgt = t
                path.append(tgt)
                start_edge_data = self.get_edge_data(path[0], path[1])
                self.remove_edge(path[0], path[1])

                end_edge_data = None
                if len(path) > 2:
                    end_edge_data = self.get_edge_data(path[-2], path[-1])
                    self.remove_edge(path[-2], path[-1])

                if not nx.has_path(self, src, tgt):
                    self.add_edge(path[0], path[1], **start_edge_data)
                    if len(path) > 2:
                        self.add_edge(path[-2], path[-1], **end_edge_data)
                else:
                    for node in path[1:-1]:
                        self.remove_node(node)

    def get_sinks(self, subgraph=None):
        """
        returns all nodes with an outgoing degree of zero
        """
        nodeset = set()
        if subgraph is None:
            subgraph = self.nodes()
        for node in subgraph:
            if self.out_degree(node) == 0:
                nodeset.add(node)
        return nodeset

    def get_sources(self, subgraph=None):
        """
        returns all nodes with an incoming degree of zero
        """
        nodeset = set()
        if subgraph is None:
            subgraph = self.nodes()
        for node in subgraph:
            if self.in_degree(node) == 0:
                nodeset.add(node)
        return nodeset


def digraph_connected_components(graph, subgraph=None):
    """
    the networkx module does not support deriving connected
    components from digraphs (only simple graphs)
    this function assumes that connection != reachable
    this means there is no difference between connected components
    in a simple graph and a digraph

    Args:
        graph (networkx.DiGraph): the input graph to gather components from

    Returns:
        :class:`list` of :class:`list`: returns a list of compnents which are lists of node names
    """
    if subgraph is None:
        subgraph = set(graph.nodes())
    g = nx.Graph()
    for src, tgt in graph.edges():
        if src in subgraph and tgt in subgraph:
            g.add_edge(src, tgt)
    for node in subgraph:
        if graph.has_node(node):
            g.add_node(node)
    return nx.connected_components(g)


def pull_contigs_from_component(
    assembly, component, min_edge_trim_weight, assembly_max_paths, log=DEVNULL
):
    """
    builds contigs from the a connected component of the assembly DeBruijn graph

    Args:
        assembly (DeBruijnGraph): the assembly graph
        component (list):  list of nodes which make up the connected component
        min_edge_trim_weight (int): the minimum weight to not remove a non cutting edge/path
        assembly_max_paths (int): the maximum number of paths allowed before the graph is further simplified
        log (function): the log function

    Returns:
        :class:`Dict` of :class:`int` by :class:`str`: the paths/contigs and their scores
    """
    path_scores = {}  # path_str => score_int
    w = min_edge_trim_weight
    unresolved_components = [component]

    while unresolved_components:
        # since now we know it's a tree, the assemblies will all be ltd to
        # simple paths
        component = unresolved_components.pop(0)
        paths_est = len(assembly.get_sinks(component)) * len(assembly.get_sources(component))

        if paths_est > assembly_max_paths:
            edge_weights = sorted([e[2]['freq'] for e in assembly.all_edges(
                assembly.get_sources(component) | assembly.get_sinks(component), data=True)])
            w = max([w + 1, edge_weights[0]])

            if w > edge_weights[-1]:
                continue
            log(
                'reducing estimated paths. Current estimate is {}+ from'.format(paths_est),
                len(component), 'nodes', 'filter increase', w)
            assembly.trim_forks_by_freq(w)
            assembly.trim_noncutting_paths_by_freq(w)
            assembly.trim_tails_by_freq(w)

            unresolved_components.extend(digraph_connected_components(assembly, component))
        else:
            for source, sink in itertools.product(assembly.get_sources(component), assembly.get_sinks(component)):
                paths = list(nx.all_simple_paths(assembly, source, sink))
                for path in paths:
                    s = path[0] + ''.join([p[-1] for p in path[1:]])
                    score = 0
                    for i in range(0, len(path) - 1):
                        score += assembly.get_edge_freq(path[i], path[i + 1])
                    path_scores[s] = max(path_scores.get(s, 0), score)
    return path_scores


def filter_contigs(contigs, assembly_min_uniq=0.01):
    """
    given a list of contigs, removes similar contigs to leave the highest (of the similar) scoring contig only
    """
    filtered_contigs = {}
    # ordering: highest scoring, then longest, then aphanumeric
    for contig in sorted(contigs, key=lambda x: (-1 * x.score, -1 * len(x.seq), x.seq)):
        rseq = reverse_complement(contig.seq)
        if contig.seq in filtered_contigs or rseq in filtered_contigs:
            continue
        drop = False
        # drop all contigs that are more than 'x' percent similar to existing contigs
        for other_seq in filtered_contigs:
            kmer_length = min(len(other_seq), len(contig.seq))
            okmer_list = set(kmers(other_seq, kmer_length))

            for okmer, ckmer in itertools.product(okmer_list, set(kmers(contig.seq, kmer_length))):
                if distance.hamming(okmer, ckmer, normalized=True) < assembly_min_uniq:
                    drop = True
                    break
            if not drop:
                for okmer, ckmer in itertools.product(okmer_list, set(kmers(rseq, kmer_length))):
                    if distance.hamming(okmer, ckmer, normalized=True) < assembly_min_uniq:
                        drop = True
                        break
            if drop:
                break

        if not drop:
            filtered_contigs[contig.seq] = contig

    return list(filtered_contigs.values())


def assemble(
    sequences,
    kmer_size,
    min_edge_trim_weight=3,
    assembly_max_paths=20,
    assembly_min_uniq=0.01,
    min_complexity=0,
    log=lambda *pos, **kwargs: None,
    **kwargs
):
    """
    for a set of sequences creates a DeBruijnGraph
    simplifies trailing and leading paths where edges fall
    below a weight threshold and the return all possible unitigs/contigs

    drops any sequences too small to fit the kmer size

    Args:
        sequences (:class:`list` of :class:`str`): a list of strings/sequences to assemble
        kmer_size: see :term:`assembly_kmer_size` the size of the kmer to use
        min_edge_trim_weight: see :term:`assembly_min_edge_trim_weight`
        remap_min_match: Minimum match percentage of the remapped read (based on the exact matches in the cigar)
        remap_min_overlap: defaults to the kmer size. Minimum amount of overlap between the contig and the remapped read
        min_contig_length: Minimum length of contigs assemble to attempt remapping reads to. Shorter contigs will be ignored
        remap_min_exact_match: see :term:`assembly_min_exact_match_to_remap`
        assembly_max_paths: see :term:`assembly_max_paths`
        log (function): the log function

    Returns:
        :class:`list` of :class:`Contig`: a list of putative contigs
    """
    if not sequences:
        return []
    min_seq = min([len(s) for s in sequences])
    kmer_size = int(round(kmer_size, 0))
    min_contig_length = kwargs.pop('min_contig_length', min_seq + 1)
    remap_min_overlap = kwargs.pop('remap_min_overlap', kmer_size)
    remap_min_exact_match = kwargs.pop('remap_min_exact_match', 6)
    remap_min_match = kwargs.pop('remap_min_match', 0.95)

    if kwargs:
        raise TypeError('unrecognized keyword argument(s)', kwargs)
    assembly = DeBruijnGraph()
    for s in sequences:
        if len(s) < kmer_size:
            continue
        kmers_list = kmers(s, kmer_size)
        for kmer in kmers_list:
            assembly.add_edge(kmer[:-1], kmer[1:])
    # use the ab min edge weight to remove all low weight edges first
    nodes = list(assembly.nodes())
    for n in nodes:
        if assembly.in_degree(n) == 0 and assembly.out_degree(n) == 0:
            assembly.remove_node(n)
    # drop all cyclic components
    for component in digraph_connected_components(assembly):
        subgraph = assembly.subgraph(component)
        if not nx.is_directed_acyclic_graph(subgraph):
            log('dropping cyclic component', time_stamp=False)
            for node in subgraph.nodes():
                assembly.remove_node(node)
    # initial data cleaning
    assembly.trim_forks_by_freq(min_edge_trim_weight)
    assembly.trim_tails_by_freq(min_edge_trim_weight)
    assembly.trim_noncutting_paths_by_freq(min_edge_trim_weight)

    path_scores = {}

    for component in digraph_connected_components(assembly):
        # pull the path scores
        path_scores.update(pull_contigs_from_component(
            assembly, component,
            min_edge_trim_weight=min_edge_trim_weight,
            assembly_max_paths=assembly_max_paths,
            log=log
        ))

    # now map the contigs to the possible input sequences
    log('filtering contigs by size and complexity', len(path_scores), time_stamp=False)
    contigs = []
    for seq, score in list(path_scores.items()):
        contig = Contig(seq, score)
        if len(contig.seq) >= min_contig_length and (not min_complexity or contig.complexity() >= min_complexity):
            contigs.append(contig)
    log('filtering similar contigs', len(contigs))
    # remap the input reads
    contigs = filter_contigs(contigs, assembly_min_uniq)
    log('remapping reads to {} contigs'.format(len(contigs)))

    for input_seq in sequences:
        maps_to = {}  # contig, score
        for contig in contigs:
            alignment = nsb_align(
                contig.seq,
                input_seq,
                min_overlap_percent=min(1, remap_min_overlap / len(input_seq)),  # accounts for hardclipped reads which may be short
                min_match=remap_min_match,
                min_consecutive_match=remap_min_exact_match
            )
            if len(alignment) != 1:
                continue
            if _cigar.match_percent(alignment[0].cigar) < remap_min_match:
                continue
            maps_to[contig] = alignment[0]
        if maps_to:
            scores = []
            for contig, read in maps_to.items():
                score = calculate_alignment_score(read)
                scores.append((contig, read, score, read.reference_end - read.reference_start))
            max_score = max([(t[2], t[3]) for t in scores])[0:2]
            best_alignments = []
            for contig, read, score1, score2 in scores:
                if max_score == (score1, score2):
                    best_alignments.append((contig, read))
            assert len(best_alignments) >= 1
            for contig, read in best_alignments:
                contig.add_mapped_sequence(read, len(best_alignments))
    log('assemblies complete')
    return contigs


def kmers(s, size):
    """
    for a sequence, compute and return a list of all kmers of a specified size

    Args:
        s (str): the input sequence
        size (int): the size of the kmers

    Returns:
        :class:`list` of :class:`str`: the list of kmers

    Example:
        >>> kmers('abcdef', 2)
        ['ab', 'bc', 'cd', 'de', 'ef']
    """
    kmers = []
    for i in range(0, len(s)):
        if i + size > len(s):
            break
        kmers.append(s[i:i + size])
    return kmers
