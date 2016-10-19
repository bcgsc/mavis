import pysam
import itertools
import warnings
from structural_variant.constants import *
from Bio.Seq import Seq
import networkx as nx


class BamCache:
    """
    caches reads by name to facilitate getting read mates without jumping around
    the file if we've already read that section
    """
    def __init__(self, bamfile):
        self.cache = {}
        self.fh = bamfile
        if not hasattr(bamfile, 'fetch'):
            self.fh = pysam.AlignmentFile(bamfile, 'rb')

    def add_read(self, read):
        if read.query_name not in self.cache:
            self.cache[read.query_name] = set()
        self.cache[read.query_name].add(read)

    def reference_id(self, chrom):
        tid = self.fh.get_tid(chrom)
        if tid == -1:
            raise KeyError('invalid reference name not present in bam file')
        return tid

    def chr(self, read):
        return self.fh.get_reference_name(read.reference_id)

    def fetch(self, chrom, start, stop, **kwargs):
        """
        wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
        position from within the loop. Also caches reads if requested and can return a limited read number
        """
        read_limit = kwargs.pop('read_limit', None)
        cache = kwargs.pop('cache', False)
        cache_if = kwargs.pop('cache_if', lambda x: True)
        result = []
        count = 0
        for read in self.fh.fetch(chrom, start, stop):
            if read_limit is not None and count >= read_limit:
                break
            result.append(read)
            if cache:
                if cache_if(read):
                    self.add_read(read)
        return result

    def get_mate(self, read, primary_only=True, allow_file_access=True):
        # NOTE: will return all mate alignments that have been cached
        putative_mates = self.cache.get(read.query_name, set())
        if SUFFIX_DELIM in read.query_name:
            prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
            putative_mates.update(self.cache.get(temp, set()))
        mates = []
        for mate in putative_mates:
            if any([read.is_read1 == mate.is_read1,
                    read.is_read2 == mate.is_read2,
                    read.next_reference_start != mate.reference_start,
                    read.next_reference_id != mate.reference_id,
                    primary_only and mate.is_secondary and mate.is_mapped]):
                continue
            mates.append(mate)
        if len(mates) == 0:
            if not allow_file_access:
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
        self.fh.close()


class CigarTools:

    @classmethod
    def recompute_cigar_mismatch(cls, read, ref):
        """
        for cigar tuples where M is used, recompute to replace with X/= for increased
        utility and specificity

        CIGAR VALUES
        M 0 alignment match (can be a sequence match or mismatch)
        I 1 insertion to the reference
        D 2 deletion from the reference
        N 3 skipped region from the reference
        S 4 soft clipping (clipped sequences present in SEQ)
        H 5 hard clipping (clipped sequences NOT present in SEQ)
        P 6 padding (silent deletion from padded reference)
        = 7 sequence match
        X 8 sequence mismatch

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
        if i + size >= len(s):
            break
        kmers.append(s[i:i + size])
    return kmers


class Contig:

    def __init__(self, sequence, score):
        self.seq = sequence
        self.remapped_reads = {}
        self.score = score

    def __hash__(self):
        return hash((self.seq, self.score))

    def add_mapped_read(self, read, initial_reads, multimap=1):
        k = (read, tuple(list(initial_reads)))
        if k in self.remapped_reads and self.remapped_reads[k] != 1 / multimap:
            raise AttributeError('cannot specify a read twice with a different value')
        self.remapped_reads[k] = 1 / multimap

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

    def add_edge(self, n1, n2):
        self.edge_freq[(n1, n2)] = self.edge_freq.get((n1, n2), 0) + 1
        nx.DiGraph.add_edge(self, n1, n2)

    def remove_edge(self, n1, n2):
        del self.edge_freq[(n1, n2)]
        nx.DiGraph.remove_edge(self, n1, n2)


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


def assemble(sequences, kmer_size=None, min_edge_weight=3, min_match_quality=0.95, min_read_mapping_overlap=None):
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
    min_seq = min([len(s) for s in sequences])
    kmer_size = int(min_seq * 0.8) if kmer_size is None else kmer_size
    min_read_mapping_overlap = kmer_size if min_read_mapping_overlap is None else min_read_mapping_overlap

    assembly = DeBruijnGraph()

    if kmer_size > min_seq:
        kmer_size = min_seq
        warnings.warn(
            'cannot specify a kmer size larger than one of the input sequences. reset to {0}'.format(min_seq))

    for s in sequences:
        for kmer in kmers(s, kmer_size):
            l = kmer[:-1]
            r = kmer[1:]
            assembly.add_edge(l, r)

    # n = 'tmp_assembly_' + str(int(random.random()*1000000))+'.txt'
    # with open(n, 'w') as fh:
    #    print('writing', n)
    #    fh.write('source\ttarget\tlabel\n')
    #    for src, tgt in assembly.edges():
    #        fh.write('{0}\t{1}\t{2}\n'.format(src, tgt, assembly.edge_freq[(src, tgt)]))
    # if we assume even coverage then we can assume the number of times an edge must be
    # visited is approximated by its weight
    # however this is not necessarily a fair assumption when we are assembling selected reads and not an entire genome
    # so here it would be better to have an abstract way to represent cycles where possible
    # for example ATCGCGCGAATG would become AT{CG}x3AATG or from the graph AT{CG}x?AATG
    # for now we don't care to assemble regions with repeats

    if not nx.is_directed_acyclic_graph(assembly):
        NotImplementedError('assembly not supported for cyclic graphs')

    # now just work with connected components
    results = {}
    # trim all paths from sources or to sinks where the edge weight is low
    for n in list(assembly.nodes()):
        if not assembly.has_node(n):
            continue
        if assembly.in_degree(n) == 0 and assembly.out_degree(n) == 0:
            assembly.remove_node(n)
        elif assembly.in_degree(n) == 0:
            # follow until the path forks or we run out of low weigh edges
            curr = n
            while assembly.out_degree(curr) == 1:
                curr, nextt = assembly.out_edges(curr)[0]
                if assembly.edge_freq[(curr, nextt)] < min_edge_weight:
                    assembly.remove_node(curr)
                    curr = nextt
                else:
                    break
        elif assembly.out_degree(n) == 0:
            curr = n
            while assembly.in_degree(curr) == 1:
                prev, curr = assembly.in_edges(curr)[0]
                if assembly.edge_freq[(prev, curr)] < min_edge_weight:
                    assembly.remove_node(curr)
                    curr = prev
                else:
                    break

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
        if len(sources) * len(sinks) > 1:
            print('source/sink combinations:', len(sources) * len(sinks))

        for source, sink in itertools.product(sources, sinks):
            for path in nx.all_simple_paths(assembly, source, sink):
                s = path[0] + ''.join([p[-1] for p in path[1:]])
                score = 0
                for i in range(0, len(path) - 1):
                    score += assembly.edge_freq[(path[i], path[i + 1])]
                results[s] = max(results.get(s, 0), score)

    for seq, score in list(results.items()):
        results[seq] = Contig(seq, score)
        if seq in sequences:
            del results[seq]

    # remap the input reads
    for seq in sequences:
        maps_to = {}  # contig, score
        for contig in results:
            if len(contig) < len(seq):
                continue
            a = nsb_align(
                contig, seq, min_overlap_percent=min_read_mapping_overlap / len(seq))
            if len(a) != 1:
                continue
            a = a[0]
            if CigarTools.match_percent(a.cigar) < min_match_quality:
                continue
            maps_to[contig] = a
        for contig, read in maps_to.items():
            results[contig].add_mapped_read(
                read, sequences[seq], len(maps_to.keys()))

    return results.values()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
