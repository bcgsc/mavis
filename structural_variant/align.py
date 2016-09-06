import pysam
import itertools
from copy import copy
import re
from structural_variant.constants import *
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
import time
import structural_variant.annotate as ann

def score_cigar(cigar_tuples):
    """
    for an input set of cigar tuples computes scoring minimums and returns a measure describing
    the alignment
    """
    putative_consecutive_matches = [[0, False]]
    total_events = 0

    for cigar, freq in cigar_tuples:
        if cigar == CIGAR.M:
            raise AttributeError('input cigar cannot have ambiguous values \'M\' is not allowed')
        if cigar != CIGAR.EQ and cigar != CIGAR.S:
            total_events += 1
        if cigar == CIGAR.EQ:
           putative_consecutive_matches[-1][0] += freq
        elif putative_consecutive_matches[-1][1]:
            pass
        else:
            putative_consecutive_matches.append([0, False])
    exact = [0]
    fuzzy = [0]
    for l, event in putative_consecutive_matches:
        if event:
            fuzzy.append(l)
        else:
            exact.append(l)
    return max(exact), max(fuzzy), total_events # longest exact, longest fuzzy, number of events

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

def str_cigar(read):
    """
    >>> read = pysam.AlignedSegment()
    >>> read.query_sequence = 'ATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGTGTTCCAAATGTGTTGGTTCTGCCTAGTATACTGATTTAGAAATTCCACCACATGGAGGTTATTTGAG'
    >>> read.cigar = [(7, 55), (4, 70)]
    >>> str_cigar(read)
    '-----------------------------------------------------------------ATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTagtgttccaaatgtgttggttctgcctagtatactgatttagaaattccaccacatggaggttatttgag--------------------------------------------------'
    """
    index = 0
    s = ''
    clipped = ''
    clipped_left = False
    max_str_len = 120
    primary_tuples = read.cigar

    typ, freq = read.cigar[0]
    typ_end, freq_end = read.cigar[-1]

    if (typ == CIGAR.S and typ_end != CIGAR.S) \
        or (typ == CIGAR.S and freq > freq_end): # soft clipped on the left side
        freq = read.cigar[0][1]
        clipped = read.query_sequence[0:freq]
        primary_tuples = read.cigar[1:]
        clipped_left = True
        index += freq
    elif typ_end == CIGAR.S: # soft clipped at the end
        freq = len(read.query_sequence) - read.cigar[-1][1]
        clipped = read.query_sequence[freq:]
        primary_tuples = read.cigar[:-1]
    
    for cigar, freq in primary_tuples:
        end = freq + index
        temp = read.query_sequence[index:end]
        if cigar == CIGAR.S:
            temp = temp.lower()
        index += freq
        s += temp
    
    if clipped_left:
        return clipped.rjust(max_str_len, '-').lower() + s.ljust(max_str_len, '-')
    return s.rjust(max_str_len, '-') + clipped.ljust(max_str_len, '-').lower()

def compute_cigar(ref, seq, **kwargs):
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

    ASSUMES if force_end or force_start are given that all bases outside are SOFT CLIPPED
    # GTGAGTAAATTCAACATCGTTTTT
    # AACTTAGAATTCAAC---------
    >>> compute_cigar('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
    [(4, 7), (7, 8)]
    >>> compute_cigar('TTTGCCATTGTTTGCTTTTGGAATCCAGGACACATCTTTAGGGCATGGCATTCCAGCTAACATCTTGGATGCCCCAGTGAAGATTTCAGTGCTTTTTTCCAATTTGTTTTTAGTATGCCTTGGCAAAAAAGAAAAGAGATACTAATTGATCAGAAAAAATGTTTTAGACTGCTAGCCCTGCTGTCTTTGGGGAAGTTGTATGCAGTGAGTAAATTCAACATCGTTTTTGGCCTCCCTATCAGTCATTAATGTAAAGTGGGGAGGCAGCTATTGCAGGCCACTATGATTTTGATAGTCAAGTAAAAGCTATGTTTTTTTGTTGCTGTTTGTTTATATCCATTAAGGGGAAAAATGGCCAGGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGGGGGAGGCCAAGGCAGGAGGATCACTTGAGACCAGGAGTTTGAGACCACCCTGGGCAACATAGTGAGACCCTGTCTTTTTAAAGATAAATAAAGATTCTTTAAGAAAATAGAAAAGGAAAGGGGGAAAAATAATCCTCCTCACAGAATATTTGCAGTTCTTCTGTATGGAGAGAGGTAACCAAAACACATAATTATCTTCAGAATTTAGGATCATTGCTCTTTTTAAATTACATTCTATCCACGGGGATTTCTCAACAACTCAAAGGAGTAGTGTTGATTTATGGAGCTCAGCCCCTAGCAGTGTGCTAAAGCCCTAAGAAACCTGGCTCTGTTTTCATTTCTTGAGTTTAGCTCAAAAGAAGCTGCTTTTGGAGACATCTTAGGATATAGAACCCTATAAACTATATGTTAGCCAGAGGGCTTTTTTTGTTGTTGTTATTTAGAAGAGGTTTATTATAGGCTAACTATATAAGGTTATTTCAGCTGTGATGATCCAGCAGTTTTGTTGAACTTAAAAGAGCCTACCTATTAAGGATGCTTTATCGTGATGTAAAGAAGATGGTGCCTTGGTTAGTG', '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AGGGCTTTTTTTGTTGTTGTTATTT---------------------------------------------------------------------------------------------------------------------------------------')
    [(7, 25)]
    """
    if len(ref) != len(seq):
        raise AssertionError('input strings must be of equal length', ref, seq, len(ref), len(seq))
    assert(len(ref) > 0)
    
    force_start = kwargs.pop('force_start', None)
    force_end  = kwargs.pop('force_end', None)
    
    seq = Seq(str(seq), DNA_ALPHABET)
    ref = Seq(str(ref), DNA_ALPHABET)
    # ---XXXX gives you 3 as the start pos
    start = len(seq) - len(seq.lstrip(seq.alphabet.gap_char)) # don't need to add one b/c 0-indexed
    end = len(seq.rstrip(seq.alphabet.gap_char)) #----XXXX-- 10 - 2 = 8; exclusive 0-indexed range

    if force_start is not None and (force_start < start or force_start >= end):
        raise AttributeError('cannot force a start position unless it if >= to the aligned start position')
    if force_end is not None and (force_end > end or force_end <= start): 
        raise AttributeError('cannot force an end position past the length of the putative alignment')
    
    cigar_tuples = []
    if force_start is not None and force_start != start:
        cigar_tuples.append((CIGAR.S, force_start - start))
        start = force_start

    tail_sc = 0
    if force_end is not None and force_end != end:
        tail_sc = end - force_end
        end = force_end

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
    if tail_sc != 0:
        cigar_tuples.append((CIGAR.S, tail_sc))

    if force_start is None:
        # turn any X on the front into soft-clipping
        if cigar_tuples[0][0] == CIGAR.X:
            cigar_tuples[0] = (CIGAR.S, cigar_tuples[0][1])
    if force_end is None:
        if cigar_tuples[-1][0] == CIGAR.X:
            cigar_tuples[-1] = (CIGAR.S, cigar_tuples[-1][1])
    temp = sum([x[1] for x in cigar_tuples])
    if temp != len(seq.strip(GAP)):
        raise AssertionError('cigar tuple must be same length as input string', temp, len(ref), ref, cigar_tuples)
    return cigar_tuples

def build_string_from_reverse_path(ref, seq, path):
    """
    >>> ref = "-mxabdce"
    >>> seq = "-abc"
    >>> build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3)])
    ('e', '-')
    >>> build_string_from_reverse_path(ref, seq, [(7, 3), (6, 3), (5,2), (4,2)])
    ('dce', '-c-')
    >>> build_string_from_reverse_path(ref, seq, [(6, 3), (5,2), (4,2), (3, 1), (2, 0)])
    ('mxabdce', '--ab-c-')
    """
    ref = str(ref)
    seq = str(seq)
    ref_string = []
    seq_string = []
    transitions = []

    end = (len(ref) -1, len(seq) - 1)
    if path[0] != end:
        path.insert(0, end)

    if path[-1][0] == 0 or path[-1][1] == 0:
        if path[-1] != (0, 0):
            path.append((0, 0))
    
    for i in range(0, len(path)):
        if i == 0:
            continue
        r, c = path[i]
        r0, c0 = path[i-1]
        up = r0 - r
        left = c0 - c
        if left == 0: # moving on ref and not seq
            seq_string.append('-'*up)
            temp = ref[r+1:r0+1]
            if up != len(temp):
                raise AssertionError('must be equal length', temp, up, path)
            ref_string.append(temp)
        elif up == 0:
            ref_string.append('-'*left)
            temp = seq[c+1:c0+1] 
            seq_string.append(temp)
            if left != len(temp):
                raise AssertionError('must be equal length', temp, left, path)
        else:
            if up != 1 or left != 1:
                raise AssertionError('diagonal moves must be up ({4}) one and left ({5}) one: {0},{1} => {2},{3}'.format(
                    r0, c0, r, c, up, left), path, i)
            ref_string.append(ref[r0])
            seq_string.append(seq[c0])
    ref_string.reverse()
    seq_string.reverse()
    return ''.join(ref_string), ''.join(seq_string)

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

    """
    max_inner_events = kwargs.pop('max_inner_events', 3)
    
    time_start = time.clock()
    size = (len(input_ref) + 1) * (len(input_seq) + 1) * 64 / 1000000
    
    if size > 10 and size < 2000:
        warnings.warn('warning the pairwise alignment matrix is large: {0} MB'.format(size))
    elif size > 2000:
        raise UserWarning('warning the pairwise alignment matrix is too large: {0} GB'.format(size/1000))
    
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
            if j == len(seq) - 1:
                highest_score = max(arr[i][j], highest_score)

    # pull out the putative paths from the matrix
    paths = []
    for i in range(0, len(ref)): # must consume the smaller sequence
        j = len(seq) - 1
        if arr[i][j] == highest_score:
            paths.append([(i, j)])
    
    complete_paths = []
    while len(paths) > 0:
        current_path = paths.pop(0)
        s, t = current_path[-1] # last element of the path
        if s == 0 or t == 0:
            complete_paths.append(current_path)
            continue
        # check if we should continue with the current path based on the input
        if max_inner_events is not None and len(current_path) > max_inner_events + 1:
            # calculate the cigar string for the current path
            temp = build_string_from_reverse_path(ref, seq, current_path)
            c = compute_cigar(*temp)
            exact, fuzzy, events = score_cigar(c)
            if events > max_inner_events:
                continue
        for next_step in pointers[(s, t)]:
            paths.append(current_path + [next_step])

        # check
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
        #path.reverse()
        ref_path, seq_path = build_string_from_reverse_path(ref, seq, path)

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
    
    if len(path_scores) == 0:
        return []
    # choose the best alignment by score and then rightmost, and then sorted alphanumerically
    best_score = max([p[0] for p in path_scores])
    
    result = []
    for score, ref, seq, path in path_scores:
        if score != best_score:
            continue
        start = len(seq) - len(seq.lstrip(GAP)) # number of gaps, and 0-index
        
        a = pysam.AlignedSegment()
        a.query_sequence = seq.strip(GAP)
        a.reference_start = start 
        a.cigar = compute_cigar(ref, seq)
        result.append(a)
    t = (time.clock() - time_start)
    if t > 5:
        print('sw_pairwise_alignment:', t)
        print('# of paths investigated', len(path_scores))
        print('matrix size', (len(input_ref)+ 1) * (len(input_seq) + 1)/1000, 'K')
    return result

def recompute_cigar(ev, read):
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
    ref = ann.HUMAN_REFERENCE_GENOME[ev.convert_index_to_chr[read.reference_id]].seq
    
    ref_pos = read.reference_start
    seq_pos = 0
    
    for cigar_value, freq in read.cigar:
        #print('current', cigar_value, freq)
        #print('ref_pos, seq_pos', ref_pos, seq_pos)
        #print(temp)
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
    return temp

def breakpoint_pos(read, orient=ORIENT.NS):
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
