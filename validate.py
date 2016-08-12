import pysam
import re
from constants import *
from sv import Breakpoint, BreakpointPair
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Alphabet import Gapped
from Bio import SeqIO
import numpy as np
from Bio.Data.IUPACData import ambiguous_dna_values
import warnings

def compute_cigar(ref, seq):
    """
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
    assert(len(ref) == len(seq) and len(ref) > 0)
    seq = Seq(str(seq), DNA_ALPHABET)
    ref = Seq(str(ref), DNA_ALPHABET)
    # ---XXXX gives you 3 as the start pos
    start = len(seq) - len(seq.lstrip(seq.alphabet.gap_char)) # don't need to add one b/c 0-indexed
    end = len(seq.rstrip(seq.alphabet.gap_char)) #----XXXX-- 10 - 2 = 8; exclusive 0-indexed range
    
    cigar_tuples = []
    # doing the first one will be easier than checking every time
    
    for i in range(start, end):
        current_mode = CIGAR.X
        if ref[i] == GAP and seq[i] == GAP:
            continue
        elif ref[i] == GAP:
            current_mode = CIGAR.I
        elif seq[i] == GAP:
            current_mode = CIGAR.D
        elif seq.alphabet.match(ref[i], seq[i]):
            current_mode = CIGAR.EQ
        
        if len(cigar_tuples) > 0:
            mode, freq = cigar_tuples[-1]
            if mode == current_mode:
                cigar_tuples[-1] = (current_mode, freq + 1)
            else:
                cigar_tuples.append((current_mode, 1))
        else:
            cigar_tuples.append((current_mode, 1))
    return cigar_tuples

def sw_pairwise_alignment(input_ref, input_seq, **kwargs):
    """
    aligns a short sequence to a relatively short reference sequence
    basic smith-waterman alignment for a pair of sequences (pairwise2 in biopython appears broken)

    arr[i, j] = max of the following
    - arr[i-1,j-1] + MATCH/MISMATCH PENALTY 
    - arr[i-1,j] + GAP_PENALTY
    - arr[i,j-1] + GAP_PENALTY

    deviates from the normal after computing the matrix, for all the best alignments recomputes their 
    scores ignoring consecutive gaps on either end (not penalized) but still penalizing interior gaps

    >>> a = sw_pairwise_alignment('ATGGACTCGGTAAA', 'CGGTAA')
    >>> a.reference_start
    7
    >>> a.query_sequence
    'CGGTAA'
    >>> a.cigar
    [(7, 6)]
    """
    size = (len(input_ref) + 1) * (len(input_seq) + 1) * 64 / 1000000000
    
    if size > 1 and size < 3:
        warnings.warn('warning the pairwise alignment matrix is large: {0} GB'.format(size))
    elif size > 1:
        raise UserWarning('warning the pairwise alignment matrix is too large: {0} GB'.format(size))
    
    ref = '-' + (input_ref.seq if hasattr(input_ref, 'seq') else input_ref)
    seq = '-' + (input_seq.seq if hasattr(input_seq, 'seq') else input_seq)
    
    arr = np.zeros(shape=(len(ref),len(seq)), dtype=np.int)

    MISMATCH = -1
    MATCH = 2
    GAP_PENALTY = -4
    highest_score = 0
    pointers = {} # key by location tuple
    # compute the dynamic programming matrix for the alignment
    for i in range(1, len(ref)):
        for j in range(1, len(seq)):
            pointers[(i, j)] = []
            matched = DNA_ALPHABET.match(ref[i], seq[j])
            diag = arr[i-1][j-1] + MATCH if matched else MISMATCH
            right = arr[i][j-1] + GAP_PENALTY
            down = arr[i-1][j] + GAP_PENALTY

            m = max([diag, right, down])
            arr[i][j] = m
            
            if diag == m:
                pointers[(i, j)].append((i-1, j-1))
            if right == m:
                pointers[(i, j)].append((i, j-1))
            if down == m:
                pointers[(i, j)].append((i-1, j))
            if i == len(ref) -1 or j == len(seq) - 1:
                highest_score = max(arr[i][j], highest_score)
    #with open('alignment_matrix.txt', 'w') as fh:
    #    fh.write(' '.join(seq) + '\n')
    #    for i in range(1, len(ref)):
    #        fh.write(ref[i] + ' ' + ' '.join([str(x) for x in arr[i]]) + '\n')
    # pull out the putative paths from the matrix
    paths = []
    for i in range(0, len(ref)):
        j = len(seq) - 1
        if arr[i][j] == highest_score:
            paths.append([(i, j)])
    for j in range(0, len(seq)):
        i = len(ref) - 1
        if arr[i][j] == highest_score:
            paths.append([(i, j)])
    
    complete_paths = []
    while len(paths) > 0:
        current_path = paths.pop(0)
        s, t = current_path[-1] # last element of the path
        if s == 0 or t == 0:
            complete_paths.append(current_path)
            continue
        for next_step in pointers[(s, t)]:
            paths.append(current_path + [next_step])
        if len(pointers[(s, t)]) == 0:
            print('ERROR: no back pointers')
            print('ref', ref)
            print('seq', seq)
            print('left', arr[s-1][t])
            print('up', arr[s][t-1])
            print('diag', arr[s-1][t-1])
            print('current', arr[s][t])
            print('coord', s, t, ref[s], seq[t])
            assert(len(next_step) != 0)
    
    # re-score the paths to give lower penalties to start and end gaps
    path_scores = []
    for path in complete_paths:
        ref_path = ''
        seq_path = ''
        path.reverse()
        path.append((len(ref) - 1, len(seq) - 1)) # bottom corner
        r0 = 0
        c0 = 0
        for r, c in path:
            down = r - r0
            right = c - c0
            if right == 0: # moving on ref and not seq
                seq_path += '-'*down
                temp = ref[r0+1:r+1]
                assert(down==len(temp))
                ref_path += temp
            elif down == 0:
                ref_path += '-'*right
                seq_path += seq[c0+1:c+1]
            else:
                assert(down == 1 and right == 1)
                ref_path += ref[r]
                seq_path += seq[c]
            r0 = r
            c0 = c
        alignment_score = 0
        for i in range(0, len(ref_path)):
            if ref_path[i] == GAP or seq_path[i] == GAP:
                alignment_score += GAP_PENALTY
            else:
                alignment_score += MATCH if DNA_ALPHABET.match(ref_path[i], seq_path[i]) else MISMATCH
        temp = max([len(ref_path) - len(ref_path.lstrip(GAP)),
                    len(seq_path) - len(seq_path.lstrip(GAP))])
        alignment_score -= (temp - 1)*GAP_PENALTY
        temp = max([len(ref_path) - len(ref_path.rstrip(GAP)),
                    len(seq_path) - len(seq_path.rstrip(GAP))])
        alignment_score -= (temp - 1)*GAP_PENALTY
        path_scores.append((alignment_score, ref_path, seq_path, path))
    
    # choose the best alignment by score and then rightmost, and then sorted alphanumerically
    best_score = max([p[0] for p in path_scores])
    alignment = sorted([(ref, seq) for score, ref, seq, path in path_scores if score == best_score],
            key=lambda x: (-1*len(x[1].rstrip(GAP)), x[1], x[0]))
    rseq, sseq = alignment[0]
    start = len(sseq) - len(sseq.lstrip(GAP)) # number of gaps, and 0-index
    a = pysam.AlignedSegment()
    a.query_sequence = sseq.strip(GAP)
    a.reference_start = start 
    a.cigar = compute_cigar(rseq, sseq)
    return a

class Evidence:
    @property
    def break1(self):
        return self.breakpoint_pair.break1
    
    @property
    def break2(self):
        return self.breakpoint_pair.break2

    def _window(self, breakpoint):
        temp = self.read_length*2 + self.average_insert_size
        start = breakpoint.start - temp - self.call_error - self.read_length - 1
        end = breakpoint.end + temp + self.call_error + self.read_length - 1
        
        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + self.call_error + self.read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - self.call_error - self.read_length - 1
        return (start, end)
    
    def __init__(self, breakpoint_pair, **kwargs):
        """
        @param =average_insert_size \a optional (type: int; default: 450)
        @param =min_anchor_size \a optional (type: int; default: 5)
        @param =min_mapping_quality \a optional (type: int; default: 20)
        @param =read_length \a optional (type: int; default: 125)
        """
        self.read_length = kwargs.pop('read_length', 125)
        self.average_insert_size = kwargs.pop('average_insert_size', 450)
        self.call_error = kwargs.pop('call_error', 10)
        self.min_anchor_size = kwargs.pop('min_anchor_size', 5)
        self.min_splits_reads_resolution = kwargs.pop('min_splits_reads_resolution', 2)
        self.convert_chr_to_index = kwargs.pop('convert_chr_to_index')
        
        self.breakpoint_pair = breakpoint_pair
        self.split_reads = {
                self.breakpoint_pair.break1: set(), 
                self.breakpoint_pair.break2: set()
                }
        self.flanking_reads = {
                self.breakpoint_pair.break1: set(), 
                self.breakpoint_pair.break2: set()
                }
    
    def add_flanking_read(self, read):
        if read.is_unmapped or read.mate_is_unmapped:
            raise UserWarning('input read (and its mate) must be mapped')
        w1 = self._window(self.break1)
        w1 = (w1[0] + 1, w1[1] + 1) # correct for psyam using 0-based coordinates
        w2 = self._window(self.break2)
        w2 = (w2[0] + 1, w2[1] + 1) # correct for psyam using 0-based coordinates
        
        # check if the strands are compatible
        if self.breakpoint_pair.opposing_strands is not None:
            opp = read.is_reverse == read.mate_is_reverse # same b/c natural orientation is opposite for paired reads
            if opp != self.breakpoint_pair.opposing_strands:
                raise UserWarning('strands are not compatible with expected strands')
        # check if this read falls in the first breakpoint window
        added_flanking = False
        if read.reference_start >= w1[0] and read.reference_end <= w1[1] and read.reference_name == self.break1.chr:
            if read.next_reference_start >= w2[0] and read.next_reference_start <= w2[1] \
                    and read.next_reference_name == self.break2.chr:
                self.flanking_reads[self.break1].add(read)
                added_flanking = True
        if read.reference_start >= w2[0] and read.reference_end <= w2[1] and read.reference_name == self.break2.chr:
            if read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                    and read.next_reference_name == self.break1.chr:
                self.flanking_reads[self.break2].add(read)
                added_flanking = True
        
        if not added_flanking:
            raise UserWarning('does not map to the expected regions. does not support the current breakpoint pair')

    def add_split_read(self, read, first_breakpoint=True):
        
        breakpoint = self.break1 if first_breakpoint else self.break2
        opposite_breakpoint = self.break2 if first_breakpoint else self.break1
        s, t = self._window(breakpoint)
        s += 1 # correct for pysam using 0-based coordinates
        t += 1 # correct for pysam using 0-based coordinates
        if read.reference_start > t or read.reference_end < s \
                or read.reference_name != breakpoint.chr:
            raise UserWarning('read does not map within the breakpoint evidence window')
        if len(read.query_sequence) - (read.query_alignment_end + 2) < self.min_anchor_size \
                and (read.query_alignment_start + 1) < self.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria')
        if self.breakpoint_pair.stranded:
            if ( read.is_reverse and breakpoint.strand == STRAND.POS ) \
                    or ( not read.is_reverse and breakpoint.strand == STRAND.NEG ):
                        raise UserWarning('split read not on the appropriate strand')
        left = ''
        right = ''
        if breakpoint.orient == ORIENT.LEFT:
            left = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            right = read.query_sequence[read.query_alignment_end:] # end is exclusive in pysam
            if len(left) < self.min_anchor_size or len(right) < self.min_anchor_size:
                raise UserWarning('split read does not meet the minimum anchor criteria')
            # now grab the soft-clipped portion on the right and try to align it to the partner region
            s, t = self._window(opposite_breakpoint)
            s -= 1
            t -= 1
            ref = HUMAN_REFERENCE_GENOME[opposite_breakpoint.chr]
            ref = ref.seq[s:t+1]
            print((s, t, ref[:10]))
            a = sw_pairwise_alignment(ref, right)
            a.query_sequence = left + a.query_sequence
            #a.query_alignment_start = len(left)
            a.cigar = [(CIGAR.S, len(left))] + a.cigar
            a.reference_start = s + a.reference_start # TODO figure out why the hell this is off by one 
            a.reference_id = self.convert_chr_to_index[opposite_breakpoint.chr]
            a.query_name = read.query_name + '--right'
            a.next_reference_start = read.next_reference_start
            a.next_reference_id = read.next_reference_id
            a.flag = read.flag
           
            print(a)
            self.split_reads[breakpoint].add(read)
            self.split_reads[opposite_breakpoint].add(a)
        elif breakpoint.orient == ORIENT.RIGHT:
            left = read.query_sequence[:read.query_alignment_start]
            right = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            if len(left) < self.min_anchor_size or len(right) < self.min_anchor_size:
                raise UserWarning('split read does not meet the minimum anchor criteria')
        else:
            raise AttributeError('cannot assign split reads to a breakpoint where the orientation has not been '
                    'specified')
        self.split_reads[breakpoint].add(read)
    
    def resolve_breakpoint(self, first_breakpoint=True):
        breakpoint = self.break1 if first_breakpoint else self.break2
        
        if breakpoint not in self.split_reads:
            raise AttributeError('invalid breakpoint. must be already associated with the evidence object through the '
                'breakpoint_pair')
        seqs = {}

        print('split reads for breakpoint', breakpoint.chr + ':' 
                + str((breakpoint.start + breakpoint.end) / 2), len(self.split_reads[breakpoint]))
        
        for read in self.split_reads[breakpoint]:
            left = ''
            right = ''
            pos = None
            
            if breakpoint.orient == ORIENT.LEFT:
                left = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
                right = read.query_sequence[read.query_alignment_end - 1:] # end is exclusive in pysam
                pos = read.reference_end - 1# pysam read coordinates are 0-based and exclusive
            elif breakpoint.orient == ORIENT.RIGHT:
                left = read.query_sequence[:read.query_alignment_start]
                right = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
                pos = read.reference_start # pysam read coordinates are 0-based
            else:
                raise AttributeError('input breakpoint cannot have an ambiguous orientation')
            if pos not in seqs:
                seqs[pos] = []
            seqs[pos].append((left, right, read))
        
        resolved_breakpoints = []
        print('invesitgating possible breakpoints', seqs.keys())
        for pos in seqs:
            if len(seqs[pos]) < self.min_splits_reads_resolution:
                continue
            print('consensus for breakpoint position', pos, breakpoint)
            mleft = max([len(l[0]) for l in seqs[pos]])
            mright = max([len(l[1]) for l in seqs[pos]])
            
            for i, s in enumerate(seqs[pos]):
                l, r, read = s
                l = str(l)
                r = str(r)
                seqs[pos][i] = SeqRecord(Seq(l.rjust(mleft, '-'), DNA_ALPHABET)) \
                        , SeqRecord(Seq(r.ljust(mright, '-'), DNA_ALPHABET)), read
            
            for l, r, read in sorted(seqs[pos], key =lambda x: x[2].query_name):
                print(l.seq, r.seq, read.query_name)
            
            align = MultipleSeqAlignment([s[0] for s in seqs[pos]], alphabet=DNA_ALPHABET)
            summary = AlignInfo.SummaryInfo(align)
            lcons = summary.dumb_consensus(consensus_alpha = DNA_ALPHABET, ambiguous='N')
            align = MultipleSeqAlignment([s[1] for s in seqs[pos]], alphabet=DNA_ALPHABET)
            summary = AlignInfo.SummaryInfo(align)
            rcons = summary.dumb_consensus(consensus_alpha = DNA_ALPHABET, ambiguous='N')
            #print(lcons, rcons)
            b = Breakpoint(
                    breakpoint.chr, 
                    pos, 
                    pos, 
                    strand = breakpoint.strand, 
                    orient = breakpoint.orient,
                    left_seq = lcons,
                    right_seq = rcons
                    ) 
            resolved_breakpoints.append((b, len(seqs[pos])))
        print('resolved_breakpoints', resolved_breakpoints)
        return resolved_breakpoints

def gather_evidence(bamfilename, breakpoint_pair, **kwargs):
    """
    @param =putative_annotations \a optional (type: List<Tuple<Transcript, Transcript>>; default: []) 
    @param =average_insert_size \a optional (type: int; default: 450)
    @param =read_cap \a optional (type: int; default: 100000)
    @param =min_anchor_size \a optional (type: int; default: 5)
    @param =min_mapping_quality \a optional (type: int; default: 20)
    @param =read_length \a optional (type: int; default: 125)
    """
    putative_annotations = kwargs.pop('putative_annotations', [])
    average_insert_size = kwargs.pop('average_insert_size', 450)
    insert_size_stdev = kwargs.pop('insert_size_stdev', 25) # deletions/ins < 50 can be captured with single reads theoretically
    read_cap = kwargs.pop('read_cap', 100000)
    min_anchor_size = kwargs.pop('min_anchor_size', 5)
    min_mapping_quality = kwargs.pop('min_mapping_quality', 20)
    read_length = kwargs.pop('read_length', 125)
    call_error = kwargs.pop('call_error', 20)
    

    or_iterable = lambda x, y: [a | b for a, b in zip(x, y)]
    bamfile = pysam.AlignmentFile(bamfilename, 'rb')
    convert_chr_to_index = {}
    for name in bamfile.references:
        convert_chr_to_index[name] = bamfile.gettid(name)
    
    # determine the region to look for evidence
    # naively classify the event if possible (will help in determining the type of support to look for)
    event_type = breakpoint_pair.naive_classification()
    # iterate over the reads in the target regions
    # look for flanking reads and split reads
    # store and summarize the evidence
    # return the 'new' breakpoint pair and the associated evidence
    # TODO adjust if annotations are given: assumes transcriptome
    
    windows = []
   
    # -1 b/c pysam uses 0-based coordinates
    for breakpoint in [breakpoint_pair.break1, breakpoint_pair.break2]:
        temp = read_length*2 + average_insert_size
        start = breakpoint.start - temp - call_error - read_length - 1
        end = breakpoint.end + temp + call_error + read_length - 1
        # narrow the evidence region if possible
        if breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - read_length - 1
        elif breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + read_length - 1
        windows.append((start, end))
    # while iterating over the reads we'll need to check for
    # support for all possible event types left based on the breakpoint orientations
    # look for soft-clipped reads

    # determine the minimum event size/isize for read pairs to support the current set of breakpoints
    evidence = Evidence(breakpoint_pair, convert_chr_to_index = convert_chr_to_index)

    for read in bamfile.fetch(
            '{0}'.format(breakpoint_pair.break1.chr), 
            windows[0][0],
            windows[0][1]):
        if read.is_unmapped:
            continue

        try:
            evidence.add_split_read(read)
        except UserWarning:
            pass
        try:
            evidence.add_flanking_read(read)
        except UserWarning:
            pass
    for read in bamfile.fetch(
            '{0}'.format(breakpoint_pair.break2.chr), 
            windows[1][0],
            windows[1][1]):
        if read.is_unmapped:
            continue

        try:
            evidence.add_split_read(read, False)
        except UserWarning:
            pass
        try:
            evidence.add_flanking_read(read)
        except UserWarning:
            pass
    b1, b1count = evidence.resolve_breakpoint()[0]
    b2, b2count = evidence.resolve_breakpoint(False)[0]

    s, t = evidence._window(evidence.break2)
    
    # try to map the soft clipping to the other window
    seq = b1.softclipped_sequence()
    ref = HUMAN_REFERENCE_GENOME[evidence.break2.chr].seq
    ref = ref[s:t+1]
    print('reference sequence', ref)
    print('softclipped sequence', seq)
    print('')
    exit()
    #alignments = pairwise2.align.localxx(ref, seq)
    #for align in alignments:
    #    print(align)
    #print('evidence', b)
    #print('evidence', evidence.resolve_breakpoint(False))
    print('flanking pairs', [ len(x) for x in evidence.flanking_reads.values()])

HUMAN_REFERENCE_GENOME = {}


def main():
    # open up the reference sequence
    global HUMAN_REFERENCE_GENOME
    #f = '/projects/seqref/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'
    f = 'chr11_chr22.fa'
    print('loading the human reference genome', f)
    with open(f, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
        print(HUMAN_REFERENCE_GENOME.keys())
    #with open('chr11_chr22.fa', 'w') as fh:
    #    SeqIO.write([HUMAN_REFERENCE_GENOME['11'], HUMAN_REFERENCE_GENOME['22']], fh, 'fasta')
    
    print('finished loading:', f)
    
    fh = open('regions.bed', 'w')
    # Breakpoint(11:128664209_128664209R?)==>Breakpoint(22:29684365_29684365L?)[EQ]
    bp = BreakpointPair(
                Breakpoint('11', 128664209, 128664209, orient = ORIENT.RIGHT),
                Breakpoint('22', 29684365, 29684365, orient = ORIENT.LEFT),
                stranded = False, opposing_strands = False)
    
    fh.write('chr{0}\t{1}\t{2}\t{3}\n'.format(bp.break1.chr, bp.break1.start, bp.break1.end, 'breakpoint1'))
    fh.write('chr{0}\t{1}\t{2}\t{3}\n'.format(bp.break2.chr, bp.break2.start, bp.break2.end, 'breakpoint2'))
    
    
    gather_evidence('/projects/analysis/analysis14/P00159/merge_bwa/125nt/hg19a/P00159_5_lanes_dupsFlagged.bam', bp)
    
    bp = BreakpointPair(
                Breakpoint('11', 128663909, 128664109, orient = ORIENT.LEFT),
                Breakpoint('22', 296843689, 296846389, orient = ORIENT.RIGHT),
                stranded = False, opposing_strands = False)
    
    fh.write('chr{0}\t{1}\t{2}\t{3}\n'.format(bp.break1.chr, bp.break1.start, bp.break1.end, 'breakpoint1'))
    fh.write('chr{0}\t{1}\t{2}\t{3}\n'.format(bp.break2.chr, bp.break2.start, bp.break2.end, 'breakpoint2'))
    
    
    gather_evidence('/projects/analysis/analysis14/P00159/merge_bwa/125nt/hg19a/P00159_5_lanes_dupsFlagged.bam', bp)
    fh.close()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
