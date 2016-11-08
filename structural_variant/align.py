import pysam
import itertools
import warnings
from structural_variant.constants import *
from Bio.Seq import Seq
import networkx as nx
from datetime import datetime


class BamCache:
    """
    caches reads by name to facilitate getting read mates without jumping around
    the file if we've already read that section
    """
    def __init__(self, bamfile):
        """
        Args:
            bamfile (str): path to the input bam file
        """
        self.cache = {}
        self.fetch_history = {}  # chr => Interval => set()
        self.fh = bamfile
        if not hasattr(bamfile, 'fetch'):
            self.fh = pysam.AlignmentFile(bamfile, 'rb')

    def add_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read to add to the cache
        """
        self.cache.setdefault(read.query_name, set()).add(read)

    def reference_id(self, chrom):
        """
        Args:
            chrom (str): the chromosome/reference name
        Returns:
            int: the reference id corrsponding to input chromosome name
        """
        tid = self.fh.get_tid(chrom)
        if tid == -1:
            raise KeyError('invalid reference name not present in bam file')
        return tid

    def chr(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read we want the chromosome name for
        Returns:
            str: the name of the chromosome
        """
        return self.fh.get_reference_name(read.reference_id)

    @classmethod
    def _generate_fetch_bins(cls, start, stop, sample_bins, bin_gap_size):
        """
        Args:
            start (int): the start if the area to fetch reads from
            stop (int): the end of the region
            sample_bins (int): the number of bins to split the region into
            bin_gap_size (int): the space to skip between bins
        """
        bin_size = int(((stop - start + 1) - bin_gap_size * (sample_bins - 1)) / sample_bins)

        fetch_regions = [(start, start + bin_size)]  # exclusive ranges for fetch
        for i in range(0, sample_bins - 1):
            st = fetch_regions[-1][1] + bin_gap_size
            end = st + bin_size
            fetch_regions.append((st, end))
        fetch_regions[-1] = fetch_regions[-1][0], stop
        return fetch_regions

    def fetch(self, chrom, start, stop, read_limit=10000, cache=False, sample_bins=3,
              cache_if=lambda x: True, bin_gap_size=0):
        """
        wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
        position from within the loop. Also caches reads if requested and can return a limited read number
        """
        # try using the cache to avoid fetching regions more than once
        result = []
        bin_limit = int(read_limit / sample_bins) if read_limit else None
        # split into multiple fetches based on the 'sample_bins'
        for fstart, fend in self.__class__._generate_fetch_bins(start, stop, sample_bins, bin_gap_size):
            count = 0
            for read in self.fh.fetch(chrom, fstart, fend):
                if bin_limit is not None and count >= bin_limit:
                    warnings.warn(
                        'hit read limit. Fetch will not cache all reads for this bin: {}-{}'.format(fstart, fend))
                    break
                result.append(read)
                if cache:
                    if cache_if(read):
                        self.add_read(read)
                count += 1
        return set(result)

    def get_mate(self, read, primary_only=True, allow_file_access=True):
        """
        Args:
            read (pysam.AlignedSegment): the read
            primary_only (boolean, default=True): ignore secondary alignments
            allow_file_access (boolean, default=True): determines if the bam can be accessed to try to find the mate
        Returns:
            List[pysam.AlignedSegment]: list of mates of the input read
        """
        # NOTE: will return all mate alignments that have been cached
        putative_mates = self.cache.get(read.query_name, set())
        # if SUFFIX_DELIM in read.query_name:
        #     prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
        #     putative_mates.update(self.cache.get(temp, set()))
        mates = []
        for mate in putative_mates:
            if any([read.is_read1 == mate.is_read1,
                    read.is_read2 == mate.is_read2,
                    read.next_reference_start != mate.reference_start,
                    read.next_reference_id != mate.reference_id,
                    primary_only and mate.is_secondary]):
                continue
            mates.append(mate)
        if len(mates) == 0:
            if not allow_file_access or read.mate_is_unmapped:
                raise KeyError('mate is not found in the cache')
            else:
                warnings.warn(
                    'looking for uncached mate of {0}. This requires file access and'.format(read.query_name) +
                    ' requests may be slow. This should also not be using in a loop iterating using the file pointer ' +
                    ' as it will change the file pointer position')
                m = self.fh.mate(read)
                self.add_read(m)
                return [m]
        return mates

    def close(self):
        """
        close the bam file handle
        """
        self.fh.close()


class CigarTools:
    """
    holds methods related to processing cigar tuples. Cigar tuples are generally
    an iterable list of tuples where the first element in each tuple is the
    CIGAR value (i.e. 1 for an insertion), and the second value is the frequency
    """
    @classmethod
    def recompute_cigar_mismatch(cls, read, ref):
        """
        for cigar tuples where M is used, recompute to replace with X/= for increased
        utility and specificity
        """
        temp = []
        offset = 0

        ref_pos = read.reference_start
        seq_pos = 0

        for cigar_value, freq in read.cigar:
            if cigar_value in [CIGAR.S, CIGAR.I]:
                temp.append((cigar_value, freq))
                seq_pos += freq
            elif cigar_value in [CIGAR.D, CIGAR.N]:
                temp.append((cigar_value, freq))
                ref_pos += freq
            elif cigar_value in [CIGAR.M, CIGAR.X, CIGAR.EQ]:
                for offset in range(0, freq):
                    if DNA_ALPHABET.match(ref[ref_pos], read.query_sequence[seq_pos]):
                        if len(temp) == 0 or temp[-1][0] != CIGAR.EQ:
                            temp.append((CIGAR.EQ, 1))
                        else:
                            temp[-1] = (CIGAR.EQ, temp[-1][1] + 1)
                    else:
                        if len(temp) == 0 or temp[-1][0] != CIGAR.X:
                            temp.append((CIGAR.X, 1))
                        else:
                            temp[-1] = (CIGAR.X, temp[-1][1] + 1)
                    ref_pos += 1
                    seq_pos += 1
            else:
                raise NotImplementedError('unexpected CIGAR value {0} is not supported currently'.format(cigar_value))
        assert(sum([x[1] for x in temp]) == sum(x[1] for x in read.cigar))
        return temp

    @classmethod
    def longest_fuzzy_match(cls, cigar, max_fuzzy_interupt=1):
        """
        computes the longest sequence of exact matches allowing for 'x' event interrupts

        """
        temp = CigarTools.join(cigar)
        longest_fuzzy_match = 0
        for pos, c in enumerate(temp):
            if c[0] != CIGAR.EQ:
                continue
            current_fuzzy_match = c[1]
            pos += 1
            fuzzy_count = 0
            while pos < len(temp) and fuzzy_count <= max_fuzzy_interupt:
                v, f = temp[pos]
                if v != CIGAR.EQ:
                    fuzzy_count += 1
                else:
                    current_fuzzy_match += f
                pos += 1
            if current_fuzzy_match > longest_fuzzy_match:
                longest_fuzzy_match = current_fuzzy_match
        return longest_fuzzy_match

    @classmethod
    def longest_exact_match(cls, cigar):
        """
        returns the longest consecutive exact match
        """
        return CigarTools.longest_fuzzy_match(cigar, 0)

    @classmethod
    def score(cls, cigar, **kwargs):
        """scoring based on sw alignment properties with gap extension penalties

        Args:
            cigar (List<(CIGAR,int)>): list of cigar tuple values
            MISMATCH (int, default=-1): mismatch penalty
            MATCH (int, default=2): match penalty
            GAP (int, default=-4): initial gap penalty
            GAP_EXTEND (int, default=-1): gap extension penalty

        Returns:
            int: the score value
        """

        MISMATCH = kwargs.pop('MISMATCH', -1)
        MATCH = kwargs.pop('MATCH', 2)
        GAP = kwargs.pop('GAP', -4)
        GAP_EXTEND = kwargs.pop('GAP_EXTEND', -1)

        score = 0
        for v, freq in cigar:
            if v == CIGAR.EQ:
                score += MATCH * freq
            elif v == CIGAR.X:
                score += MISMATCH * freq
            elif v in [CIGAR.I, CIGAR.D]:
                score += GAP + GAP_EXTEND * (freq - 1)
            elif v in [CIGAR.S, CIGAR.N]:
                pass
            else:
                raise AssertionError('unexpected cigar value', v)
        return score

    @classmethod
    def match_percent(cls, cigar):
        """
        calculates the percent of aligned bases (matches or mismatches) that are matches
        """
        matches = 0
        mismatches = 0
        for v, f in cigar:
            if v == CIGAR.EQ:
                matches += f
            elif v == CIGAR.X:
                mismatches += f
        if matches + mismatches == 0:
            return 0
        else:
            return matches / (matches + mismatches)

    @classmethod
    def join(cls, *pos):
        """
        given a number of cigar lists, joins them and merges any consecutive tuples
        with the same cigar value

        Example:
            >>> CigarTools.join([(1, 1), (4, 7)], [(4, 3), (2, 4)])
            [(1, 1), (4, 10), (2, 4)]
        """
        result = []
        for cigar in pos:
            for v, f in cigar:
                if len(result) > 0 and result[-1][0] == v:
                    result[-1] = (v, f + result[-1][1])
                else:
                    result.append((v, f))
        return result

    @classmethod
    def extend_softclipping(cls, original_cigar, min_exact_to_stop_softclipping):
        """
        given some input cigar, extends softclipping if there are mismatches/insertions/deletions
        close to the end of the aligned portion. The stopping point is defined by the
        min_exact_to_stop_softclipping parameter. this function will throw an error if there is no
        exact match aligned portion to signal stop

        Args:
            original_cigar (List[CIGAR,int]): the input cigar
            min_exact_to_stop_softclipping (int): number of exact matches to terminate extension

        Returns:
            (List[CIGAR,int], int): the new cigar string and a number representing the shift from the original start position
        """
        cigar = original_cigar[:]
        preshift = 0
        # determine how far to scoop for the front softclipping
        head = None
        while len(cigar) > 0:
            v, f = cigar[0]
            if v == CIGAR.EQ and f >= min_exact_to_stop_softclipping:
                break
            elif v in [CIGAR.D, CIGAR.N]:
                preshift -= f
                del cigar[0]
                continue
            elif head is None:
                head = (CIGAR.S, f)
            else:
                head = (CIGAR.S, head[1] + f)
            if v == CIGAR.S:
                pass
            else:
                preshift += f
            del cigar[0]
        tail = None
        while len(cigar) > 0:
            v, f = cigar[-1]
            if v == CIGAR.EQ and f >= min_exact_to_stop_softclipping:
                break
            elif v in [CIGAR.D, CIGAR.N]:
                pass
            elif tail is None:
                tail = (CIGAR.S, f)
            else:
                tail = (CIGAR.S, tail[1] + f)
            del cigar[-1]
        if len(cigar) == 0:
            raise AttributeError(
                'input does not have a long enough match sequence')
        if head is not None:
            cigar.insert(0, head)
        if tail is not None:
            cigar.append(tail)

        return cigar, preshift

    @classmethod
    def compute(cls, ref, alt, **kwargs):
        """
        given a ref and alt sequence compute the cigar string representing the alt

        returns the cigar tuples along with the start position of the alt relative to the ref
        """
        force_softclipping = kwargs.pop('force_softclipping', True)
        min_exact_to_stop_softclipping = kwargs.pop('min_exact_to_stop_softclipping', 6)

        if len(ref) != len(alt):
            raise AttributeError('ref and alt must be the same length')
        cigar = []
        for r, a in zip(ref, alt):
            if r == a and r == GAP:
                pass
            elif r == GAP:
                cigar.append((CIGAR.I, 1))
            elif a == GAP:
                cigar.append((CIGAR.D, 1))
            elif DNA_ALPHABET.match(r, a):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))
        cigar = cls.join(cigar)

        has_terminator = False
        for v, f in cigar:
            if v == CIGAR.EQ and f >= min_exact_to_stop_softclipping:
                has_terminator = True
                break

        if not force_softclipping or not has_terminator:
            if not has_terminator and force_softclipping:
                warnings.warn(
                    'could not force softclipping. did not find a long enough match sequence to terminate clipping')
            return cigar, 0

        start_pos = 0
        if cigar[0][0] in [CIGAR.D, CIGAR.N]:
            start_pos += cigar[0][1]
            del cigar[0]
        cigar, shift = cls.extend_softclipping(
            cigar, min_exact_to_stop_softclipping)
        return cigar, start_pos + shift

    @classmethod
    def convert_for_igv(cls, cigar):
        """
        igv does not support the extended CIGAR values for match v mismatch

        Example:
            >>> CigarTools.convert_for_igv([(7, 4), (8, 1), (7, 5)])
            [(0, 10)]
        """
        result = []
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ]:
                v = CIGAR.M
            result.append((v, f))
        return cls.join(result)

    @classmethod
    def alignment_matches(cls, cigar):
        """
        counts the number of aligned bases irrespective of match/mismatch
        this is equivalent to counting all CIGAR.M
        """
        result = 0
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ, CIGAR.M]:
                result += f
        return result


def reverse_complement(s):
    temp = Seq(s, DNA_ALPHABET)
    return str(temp.reverse_complement())


def longest_homopolymer(sequence):
    if len(sequence) == 0:
        raise AttributeError('cannot compute longest_homopolymer on an empty sequence')
    hp = [0]
    last_char = sequence[0]

    for char in sequence:
        if char == last_char:
            hp[-1] += 1
        else:
            hp.append(0)
        last_char = char
    return max(hp)


def breakpoint_pos(read, orient=ORIENT.NS):
    """
    assumes the breakpoint is the position following softclipping on the side with more
    softclipping (unless and orientation has been specified)

    Args:
        read (psyam.AlignedSegment): the read object
        orient (ORIENT): the orientation

    Returns:
        int: the position of the breakpoint in the input read
    """
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping')

    if orient == ORIENT.NS:
        if (typ == CIGAR.S and end_typ == CIGAR.S and freq > end_freq) \
                or typ == CIGAR.S and end_typ != CIGAR.S:
            orient = ORIENT.RIGHT
            # soft clipped to the left
        else:
            # soft clipped to the right
            orient = ORIENT.LEFT

    if orient == ORIENT.RIGHT:
        return read.reference_start
    else:
        return read.reference_end - 1


def kmers(s, size):
    """
    for a sequence, compute and return a list of all kmers of a specified size

    Args:
        s (str): the input sequence
        size (int): the size of the kmers

    Returns:
        List[str]: the list of kmers
    """
    kmers = []
    for i in range(0, len(s)):
        if i + size > len(s):
            break
        kmers.append(s[i:i + size])
    return kmers


class Contig:
    """
    """
    def __init__(self, sequence, score):
        self.seq = sequence
        self.remapped_reads = {}
        self.score = score
        self.alignments = None

    def __hash__(self):
        return hash(self.seq)

    def add_mapped_read(self, read, multimap=1):
        self.remapped_reads[read] = min(self.remapped_reads.get(read, 1), 1 / multimap)

    def remap_score(self):
        return sum(self.remapped_reads.values())


class DeBruijnGraph(nx.DiGraph):
    """
    wrapper for a basic digraph
    enforces edge weights
    """

    def __init__(self, *pos, **kwargs):
        nx.DiGraph.__init__(self, *pos, **kwargs)
        self.edge_freq = {}

    def add_edge(self, n1, n2, freq=1):
        self.edge_freq[(n1, n2)] = self.edge_freq.get((n1, n2), 0) + freq
        nx.DiGraph.add_edge(self, n1, n2)

    def remove_edge(self, n1, n2):
        del self.edge_freq[(n1, n2)]
        nx.DiGraph.remove_edge(self, n1, n2)

    def trim_low_weight_tails(self, min_weight):
        """
        for any paths where all edges are lower than the minimum weight trim
        """
        for n in list(self.nodes()):
            if not self.has_node(n):
                continue
            # follow until the path forks or we run out of low weigh edges
            curr = n
            while self.degree(curr) == 1:
                if self.out_degree(curr) == 1:
                    curr, other = self.out_edges(curr)[0]
                    if self.edge_freq[(curr, other)] < min_weight:
                        self.remove_node(curr)
                        curr = other
                    else:
                        break
                elif self.in_degree(curr) == 1:
                    other, curr = self.in_edges(curr)[0]
                    if self.edge_freq[(other, curr)] < min_weight:
                        self.remove_node(curr)
                        curr = other
                    else:
                        break
                else:
                    break
        for n in list(self.nodes()):
            if not self.has_node(n):
                continue
            if self.degree(n) == 0:
                self.remove_node(n)


def nsb_align(ref, seq, weight_of_score=0.5, min_overlap_percent=100):
    """
    given some reference string and a smaller sequence string computes the best non-space-breaking alignment
    i.e. an alignment that does not allow for indels (straight-match)
    """
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence')

    if min_overlap_percent <= 0 or min_overlap_percent > 100:
        raise AttributeError('percent must be greater than 0 and up to 100', min_overlap_percent)

    min_overlap = int(round(min_overlap_percent * len(seq) / 100, 0))

    # store to improve speed and space (don't need to store all alignments)
    best_score = 0
    results = []

    for ref_start in range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap):
        score = 0
        cigar = []
        for i in range(0, len(seq)):
            r = ref_start + i
            if r < 0 or r >= len(ref):
                cigar.append((CIGAR.S, 1))
                continue
            if DNA_ALPHABET.match(ref[r], seq[i]):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))

        cigar = CigarTools.join(cigar)
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])

        qstart = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]

        score = CigarTools.score(cigar) * weight_of_score + \
            CigarTools.longest_exact_match(cigar) * (1 - weight_of_score)
        score = CigarTools.score(cigar)
        a = pysam.AlignedSegment()
        a.query_sequence = str(seq)
        a.reference_start = ref_start + qstart
        a.cigar = cigar

        if score >= best_score:
            best_score = score
            results.append((a, score))

    filtered = [x for x, y in results if y == best_score]

    return filtered


def digraph_connected_components(graph):
    """
    the networkx module does not support deriving connected
    components from digraphs (only simple graphs)
    this function assumes that connection != reachable
    this means there is no difference between connected components
    in a simple graph and a digraph
    """
    g = nx.Graph()
    for src, tgt in graph.edges():
        g.add_edge(src, tgt)
    for n in graph.nodes():
        g.add_node(n)
    return nx.connected_components(g)


def assemble(sequences, kmer_size=None, min_edge_weight=3, min_match_quality=0.95, min_read_mapping_overlap=None,
             min_contig_length=None):
    """
    for a set of sequences creates a DeBruijnGraph
    simplifies trailing and leading paths where edges fall
    below a weight threshold and the return all possible unitigs/contigs

    Args:
        sequences (List[str]): a list of strings/sequences to assemble
        kmer_size (int): the size of the kmer to use
        min_edge_weight (int): applies to trimming (see desc)
        min_match_quality (float): percent match for re-aligned reads to contigs
        min_read_mapping_overlap (int): the minimum amount of overlap required when aligning reads to contigs

    Returns:
        List[Contig]: a list of putative contigs
    """
    if len(sequences) == 0:
        return []
    print(datetime.now(), 'assemble() start', len(sequences), )
    min_seq = min([len(s) for s in sequences])
    if kmer_size is None:
        temp = int(min_seq * 0.75)
        if temp < 10:
            kmer_size = min(min_seq, 10)
        else:
            kmer_size = temp
    elif kmer_size > min_seq:
        kmer_size = min_seq
        warnings.warn(
            'cannot specify a kmer size larger than one of the input sequences. reset to {0}'.format(min_seq))
    min_read_mapping_overlap = kmer_size if min_read_mapping_overlap is None else min_read_mapping_overlap
    min_contig_length = min_seq + 1 if min_contig_length is None else min_contig_length
    print(datetime.now(), 'assemble() start', len(sequences), 'kmer_size=', kmer_size)
    assembly = DeBruijnGraph()

    for s in sequences:
        for kmer in kmers(s, kmer_size):
            l = kmer[:-1]
            r = kmer[1:]
            assembly.add_edge(l, r)

    if not nx.is_directed_acyclic_graph(assembly):
        NotImplementedError('assembly not supported for cyclic graphs')
    
    print('assembly', len(assembly.nodes()), len(assembly.edges()), min_edge_weight)
    for s, t in sorted(assembly.edges()):
        f = assembly.edge_freq[(s, t)]
    # now just work with connected components
    # trim all paths from sources or to sinks where the edge weight is low
    assembly.trim_low_weight_tails(min_edge_weight)
    print('assembly after trim', len(assembly.nodes()))
    path_scores = {}  # path_str => score_int

    for component in digraph_connected_components(assembly):
        # since now we know it's a tree, the assemblies will all be ltd to
        # simple paths
        sources = set()
        sinks = set()
        for node in component:
            if assembly.degree(node) == 0:  # ignore singlets
                pass
            elif assembly.in_degree(node) == 0:
                sources.add(node)
            elif assembly.out_degree(node) == 0:
                sinks.add(node)
        if len(sources) * len(sinks) > 10:
            warnings.warn('source/sink combinations: {}'.format(len(sources) * len(sinks)))

        for source, sink in itertools.product(sources, sinks):
            for path in nx.all_simple_paths(assembly, source, sink):
                s = path[0] + ''.join([p[-1] for p in path[1:]])
                score = 0
                for i in range(0, len(path) - 1):
                    score += assembly.edge_freq[(path[i], path[i + 1])]
                path_scores[s] = max(path_scores.get(s, 0), score)
    print('path_scores', len(path_scores.items()))
    # now map the contigs to the possible input sequences
    contigs = {}
    for seq, score in list(path_scores.items()):
        if seq not in sequences and len(seq) >= min_contig_length:
            contigs[seq] = Contig(seq, score)

    # remap the input reads
    for input_seq in sequences:
        maps_to = {}  # contig, score
        for contig in contigs.values():
            a = nsb_align(
                contig.seq,
                input_seq,
                min_overlap_percent=min_read_mapping_overlap / len(contig.seq)
            )
            if len(a) != 1:
                continue
            if CigarTools.match_percent(a[0].cigar) < min_match_quality:
                continue
            maps_to[contig] = a[0]
        for contig, read in maps_to.items():
            contig.add_mapped_read(read, len(maps_to.keys()))
    return list(contigs.values())
