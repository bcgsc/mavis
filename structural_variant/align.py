import pysam
import random
import math
import os
import subprocess
import io
import itertools
from copy import copy
import warnings
import re
import TSV
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
import networkx as nx
from subprocess import Popen, PIPE, STDOUT
import tempfile
import codecs
from structural_variant.interval import Interval

class CigarTools:
    @classmethod
    def longest_fuzzy_match(cls, cigar, max_fuzzy_interupt=1):
        """
        computes the longest sequence of exact matches allowing for 'x' event interrupts

        >>> c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        >>> CigarTools.longest_fuzzy_match(c, 1)
        15
        >>> CigarTools.longest_fuzzy_match(c, 0)
        10
        >>> CigarTools.longest_fuzzy_match(c, 2)
        16
        >>> CigarTools.longest_fuzzy_match(c, 4)
        16
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
        """
        scoring based on sw alignment properties with gap extension penalties

        @param cigar [required] (type: List<(CIGAR,int)>) list of cigar tuple values
        @param =MISMATCH [optional] (type: int; default: -1) mismatch penalty
        @param =MATCH [optional] (type: int; default: 2) match penalty
        @param =GAP [optional] (type: int; -4) initial gap penalty
        @param =GAP_EXTEND [optional] (type: int; -1) gap extension penalty
        
        @return (type: int) the score value

        >>> c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        >>> CigarTools.score(c)
        22
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
        >>> c = [(CIGAR.S, 10), (CIGAR.EQ, 1), (CIGAR.X, 4), (CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.EQ, 5)]
        >>> CigarTools.match_percent(c)
        0.8
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
        cigar = original_cigar[:]
        preshift =  0
        # determine how far to scoop for the front softclipping
        head = None
        while True:
            if len(cigar) == 0:
                break
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
        while True:
            if len(cigar) == 0:
                break
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
            raise AttributeError('input does not have a long enough match sequence')
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
        
        >>> CigarTools.compute('GTGAGTAAATTCAACATCGTTTTT', 'AACTTAGAATTCAAC---------')
        ([(4, 7), (7, 8)], 7)
        >>> CigarTools.compute('GTGAGTAAATTCAACATCGTTTTT', '--CTTAGAATTCAAC---------')
        ([(4, 5), (7, 8)], 7)
        """
        force_softclipping = kwargs.pop('force_softclipping', True)
        min_exact_to_stop_softclipping = kwargs.pop('min_exact_to_stop_softclipping', 6)
        start_pos = 0

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
        
        if cigar[0][0] in [CIGAR.D, CIGAR.N]:
            start_pos += cigar[0][1]
            del cigar[0]
        
        has_terminator = False
        for v, f in cigar:
            if v == CIGAR.EQ and f >= min_exact_to_stop_softclipping:
                has_terminator = True
                break

        if not force_softclipping or not has_terminator:
            if not has_terminator and force_softclipping:
                warnings.warn('could not force softclipping. did not find a long enough match sequence to terminate clipping')
            return cigar
        cigar, shift = cls.extend_softclipping(cigar, min_exact_to_stop_softclipping)
        return cigar, start_pos + shift
    
    @classmethod
    def convert_for_igv(cls, cigar):
        result = []
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ]:
                v = CIGAR.M
            result.append((v, f))
        return cls.join(result)
    
    @classmethod
    def alignment_matches(cls, cigar):
        result = 0
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ, CIGAR.M]:
                result += f
        return result

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
        if cigar != CIGAR.EQ and cigar != CIGAR.S and cigar != CIGAR.N:
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
        if cigar_tuples[0][0] in [CIGAR.X, CIGAR.I]:
            cigar_tuples[0] = (CIGAR.S, cigar_tuples[0][1])
    if force_end is None:
        if cigar_tuples[-1][0] in [CIGAR.X, CIGAR.I]:
            cigar_tuples[-1] = (CIGAR.S, cigar_tuples[-1][1])
    temp = sum([x[1] for x in cigar_tuples])
    if temp != len(seq.strip(GAP)):
        raise AssertionError('cigar tuple must be same length as input string', temp, len(ref), ref, cigar_tuples)

    return cigar_tuples

def pairwise_nsb_align(ref, seq):
    """
    test for the best alignment; non-space breaking (i.e. don't allow for indels)
    """
    best_score = 0 # store to improve speed and space (don't need to store all alignments)
    results = []

    for ref_start in range(0, len(ref)):
        cigar = []
        
        if ref_start + len(seq) > len(ref):
            break

        r = ref[ref_start:ref_start + len(seq)]

        cigar = compute_cigar(r, seq)

        score = nsb_cigar_score(cigar)
        
        a = pysam.AlignedSegment()
        a.query_sequence = seq
        a.reference_start = ref_start 
        a.cigar = cigar
        
        if score >= best_score:
            best_score = score
            results.append((a, score))
    
    filtered = [x for x, y in results if y == best_score ]
    
    return filtered

def nsb_align_with_softclipping(ref, seq, **kwargs):
    weight_of_score = kwargs.pop('weight_of_score', 0.5)

    results = []

    for ref_start in range(0, len(ref)):
        cigar = []
        
        if ref_start + len(seq) > len(ref):
            break

        r = ref[ref_start:ref_start + len(seq)]

        cigar = CigarTools.compute(r, seq)

        score = CigarTools.score(cigar) * weight_of_score
        score += CigarTools.longest_exact_match(cigar) / ( len(seq) * ( 1 - weight_of_score) )
        
        preshift = 0
        for v, f in cigar:
            if v not in [CIGAR.M, CIGAR.EQ, CIGAR.X, CIGAR.D]:
                preshift += f
            else:
                break
        a = pysam.AlignedSegment()
        a.query_sequence = seq
        a.reference_start = ref_start + preshift
        a.cigar = cigar
        
        results.append((a, score))
    
    max_score = max([s for read, s in results])
    
    return [ read for read, s in results if s == max_score ] 

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
    ref = ann.HUMAN_REFERENCE_GENOME[ev.settings.convert_index_to_chr[read.reference_id]].seq
    
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

def kmers(s, size):
    """
    for a sequence, compute and return a list of all kmers of a specified size

    @param s (type: str) the input sequence
    @param size (type: int) the size of the kmers

    @return (type: List<str>) the list of kmers
    """
    kmers = []
    for i in range(0, len(s)):
        if i + size >= len(s):
            break
        kmers.append(s[i:i+size])
    return kmers

class BlatAlignedSegment(pysam.AlignedSegment):
    def __init__(self, row):
        pysam.AlignedSegment.__init__(self)
        self.blat = row
    
    def blat_score(self): # convenience
        return self.blat['score']
    
    def query_coverage_interval(self):
        query_ranges = [ (x, x + y - 1) for x, y in zip(self.blat['qstarts'], blat['block_sizes']) ]
        return Interval.union(*query_ranges)

class Blat:
    @staticmethod
    def millibad(row, is_protein=False, is_mrna=False):
        """ 
        this function is used in calculating percent identity
        direct translation of the perl code
        # https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
        """
        size_mul = 1 if not is_protein else 3
        if is_protein and is_mrna:
            raise AttributeError('cannot be a protein AND and mRNA')
        q_ali_size = size_mul * (row['qend'] - row['qstart'])
        t_ali_size = row['tend'] - row['tstart']
        ali_size = min(q_ali_size, t_ali_size)
        if ali_size <= 0:
            return 0
        
        size_dif = q_ali_size - t_ali_size
        if size_dif < 0:
            if is_mrna:
                size_dif = 0
            else:
                size_dif = abs(size_dif)
        
        insert_factor = row['qgap_count']
        if not is_mrna:
            insert_factor += row['tgap_count']

        total = size_mul * (row['match'] + row['repmatch'] + row['mismatch'])
        if total != 0:
            millibad = 0
            round_away_from_zero = 3*math.log(1 + size_dif)
            if round_away_from_zero < 0:
                round_away_from_zero = int(round_away_from_zero - 0.5)
            else:
                round_away_from_zero = int(round_away_from_zero + 0.5)
            millibad = (1000 * (row['mismatch'] * size_mul + insert_factor + round_away_from_zero)) / total
            return millibad
        else:
            return 0
    
    @staticmethod
    def score(row, is_protein=False):
        """
        direct translation from ucsc guidelines on replicating the web blat score
        https://genome.ucsc.edu/FAQ/FAQblat.html#blat4

        below are lines from the perl code i've re-written in python
        my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
        sizmul = 1 for DNA
        my $pslScore = $sizeMul * ($matches + ( $repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumInsert)
        """

        size_mul = 1 if not is_protein else 3
        score = size_mul * (row['match'] + (row['repmatch'] >> 1) ) \
                - size_mul * row['mismatch'] \
                - row['qgap_count'] \
                - row['tgap_count']
        return score
    
    @staticmethod
    def percent_identity(row, is_protein=False, is_mrna=False):
        return 100 - Blat.millibad(row, is_protein, is_mrna) * 0.1

    @staticmethod
    def read_pslx(filename, seqid_to_sequence_mapping, **kwargs):
        is_protein = kwargs.pop('is_protein', False)
        pslx_header = [
                'match', 'mismatch', 'repmatch', 'ncount', 
                'qgap_count', 'qgap_bases',
                'tgap_count', 'tgap_bases',
                'strand',
                'qname', 'qsize', 'qstart', 'qend',
                'tname', 'tsize', 'tstart', 'tend',
                'block_count', 'block_sizes',
                'qstarts', 'tstarts',
                'qseqs', 'tseqs'
                ] 
        
        split_csv_trailing_seq = lambda x: [ s.upper() for s in re.sub(',$', '', x).split(',')]
        split_csv_trailing_ints = lambda x: [ int(s) for  s in re.sub(',$', '', x).split(',')]
        
        header, rows = TSV.read_file(
                filename, 
                header=pslx_header,
                cast = {
                    'match': 'int',
                    'mismatch': 'int',
                    'repmatch': 'int',
                    'ncount': 'int',
                    'qgap_count': 'int', 
                    'qgap_bases': 'int',
                    'tgap_count': 'int', 
                    'tgap_bases': 'int',
                    'qsize': 'int', 
                    'qstart': 'int', 
                    'qend': 'int',
                    'tsize': 'int', 
                    'tstart': 'int', 
                    'tend': 'int',
                    'block_count': 'int'
                    },
                validate = {
                    'strand': '^[\+-]$'
                    },
                transform = {
                    'tname': lambda x: re.sub('^chr', '', x),
                    'block_sizes': split_csv_trailing_ints,
                    'qstarts': split_csv_trailing_ints,
                    'tstarts': split_csv_trailing_ints,
                    'qseqs': split_csv_trailing_seq,
                    'tseqs': split_csv_trailing_seq
                    }
                )
        
        for row in rows:
            row['score'] = Blat.score(row)
            row['percent_ident'] = Blat.percent_identity(row)
            qseq = seqid_to_sequence_mapping[row['qname']]
            row['qseq_full'] = qseq
        return header, rows
    
    @staticmethod
    def pslx_row_to_pysam(row, **kwargs):
        """
        given a 'row' from reading a pslx file. converts the row to a BlatAlignedSegment object

        @param row [required] (type: dict<str,*>) a row object from the 'read_pslx' method
        @param =chr_to_index [required] (type: dict<chr,int>) a mapping of chromosome names to reference_id's for the pysam file

        """
        chr_to_index = kwargs.pop('chr_to_index')
        
        chrom = chr_to_index[row['tname']]
        qseq = row['qseq_full']
        if row['strand'] == STRAND.NEG:
            temp = Seq(qseq, DNA_ALPHABET)
            temp = temp.reverse_complement()
            qseq = str(temp)
        
        # note: converting to inclusive range [] vs end-exclusive [)
        query_ranges = [ (x, x + y - 1) for x, y in zip(row['qstarts'], row['block_sizes']) ]
        ref_ranges = [ (x, x + y - 1) for x, y in zip(row['tstarts'], row['block_sizes']) ] 
        
        cigar = []
        seq = ''

        # add initial soft-clipping
        if query_ranges[0][0] > 0:
            temp = qseq[0:query_ranges[0][0]]
            seq += temp
            cigar.append((CIGAR.S, len(temp)))
        for i in range(0, len(query_ranges)):
            rcurr = ref_ranges[i]
            qcurr = query_ranges[i]
            if i > 0:
                # append based on the prev range
                rprev = ref_ranges[i - 1]
                qprev = query_ranges[i - 1]
                qjump = qcurr[0] - qprev[1]
                rjump = rcurr[0] - rprev[1]
                
                if rjump == 1: # reference is consecutive
                    if qjump > 1: # query range skipped. insertion to the reference sequence
                        cigar.append((CIGAR.I, qjump - 1))
                        seq += qseq[qprev[1] + 1:qcurr[0]] # adds the inserted seq for the pysam read
                elif qjump == 1: # query is consecutive
                    if rjump > 1: # reference range skipped. deletion of the reference sequence
                        cigar.append((CIGAR.D, rjump - 1))
                else: # indel
                    seq += qseq[qprev[1] + 1:qcurr[0]]
                    cigar.append((CIGAR.I, qjump - 1))
                    cigar.append((CIGAR.D, rjump - 1))
            # add the current range of matches
            temp = compute_cigar(row['tseqs'][i], row['qseqs'][i])
            if temp[0][0] == CIGAR.S:
                temp[0] = (CIGAR.X, temp[0][1])
            if temp[-1][0] == CIGAR.S:
                temp[-1] = (CIGAR.X, temp[-1][1])
            cigar = CigarTools.join(cigar, temp)
            seq += qseq[qcurr[0]:qcurr[1] + 1]
        
        if query_ranges[-1][1] < len(qseq) - 1:
            temp = qseq[query_ranges[-1][1] + 1:]
            seq += temp
            cigar.append((CIGAR.S, len(temp)))
        read = BlatAlignedSegment(row)
        read.query_sequence = seq
        read.reference_start = row['tstart']
        read.reference_id = chrom
        read.cigar = cigar
        read.query_name = row['qname']
        if row['strand'] == STRAND.NEG:
            read.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE
        return read

class Contig:
    def __init__(self, sequence, score):
        self.seq = sequence
        self.remapped_reads = {}
        self.score = score

    def __hash__(self):
        return hash((self.seq, self.score))

    def add_mapped_read(self, read, multimap=1):
        if read in self.remapped_reads and self.remapped_reads[read] != 1 / multimap:
            raise AttributeError('cannot specify a read twice with a different value')
        self.remapped_reads[read] = 1 / multimap

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

def percent_query_match(cigar):
    score = 0
    total = 0
    for c, v in cigar:
        if c == CIGAR.S:
            continue
        total += v
        if c == CIGAR.EQ:
            score += v
    return score/total

def nsb_overhang_align(ref, seq, **kwargs):
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence')
    min_overlap = min(kwargs.pop('min_overlap', len(seq)), len(seq))

    if len(seq) < min_overlap or len(ref) < min_overlap:
        raise AttributeError('minimum overlap cannot be greater than the length of either of the input sequences')
    
    best_score = 0 # store to improve speed and space (don't need to store all alignments)
    results = []
    
    for ref_start in range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap):
        score = 0
        cigar = []
        for i in range(0, len(seq)):
            r = ref_start + i
            if r < 0 or r >= len(ref):
                if len(cigar) == 0 or cigar[-1][0] != CIGAR.X:
                    cigar.append((CIGAR.X, 1))
                else:
                    cigar[-1] = (cigar[-1][0], cigar[-1][1] + 1)
                continue
            r = ref[r]
            q = seq[i]
            if DNA_ALPHABET.match(r, q):
                if len(cigar) == 0 or cigar[-1][0] != CIGAR.EQ:
                    cigar.append((CIGAR.EQ, 1))
                else:
                    cigar[-1] = (cigar[-1][0], cigar[-1][1] + 1)
            else:
                if len(cigar) == 0 or cigar[-1][0] != CIGAR.X:
                    cigar.append((CIGAR.X, 1))
                else:
                    cigar[-1] = (cigar[-1][0], cigar[-1][1] + 1)
        
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])
        
        score = nsb_cigar_score(cigar)
        
        a = pysam.AlignedSegment()
        a.query_sequence = str(seq)
        a.reference_start = ref_start 
        a.cigar = cigar
        
        if score >= best_score:
            best_score = score
            results.append((a, score))
    
    filtered = [x for x, y in results if y == best_score ]
    
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

def assemble(sequences, **kwargs):
    """
    for a set of sequences creates a DeBruijnGraph
    simplifies trailing and leading paths where edges fall
    below a weight threshold and the return all possible unitigs/contigs

    @param sequences (type: List<string>) a list of strings/sequences to assemble
    @param =kmer_size (type: int) the size of the kmer to use
    @param =min_edge_weight (type: int) applies to trimming (see desc)
    @param =min_match_quality (type: float) percent match for re-aligned reads to contigs
    @param =min_read_mapping_overlap (type: int) the minimum amount of overlap required when aligning reads to contigs

    @return (type: List<Contig>) a list of putative contigs
    """
    min_seq = min([len(s) for s in sequences])
    kmer_size = kwargs.pop('kmer_size', int(min_seq * 0.8) )
    min_edge_weight = kwargs.pop('min_edge_weight', 3)
    min_match_quality = kwargs.pop('min_match_quality', 0.95)
    min_read_mapping_overlap = int(kwargs.pop('min_read_mapping_overlap', min_seq * 0.8))
    assembly = DeBruijnGraph()
    
    if kmer_size > min_seq:
        kmer_size = min_seq
        warnings.warn('cannot specify a kmer size larger than one of the input sequences. reset to {0}'.format(min_seq))

    for s in sequences:
        for kmer in kmers(s, kmer_size):
            l = kmer[:-1]
            r = kmer[1:]
            assembly.add_edge(l, r)
    
    n = 'tmp_assembly_' + str(int(random.random()*1000000))+'.txt'
    with open(n, 'w') as fh:
        print('writing', n)
        fh.write('source\ttarget\tlabel\n')
        for src, tgt in assembly.edges():
            fh.write('{0}\t{1}\t{2}\n'.format(src, tgt, assembly.edge_freq[(src, tgt)]))
    # if we assume even coverage then we can assume the number of times an edge must be visited is approximated by its weight
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
    
    n = 'tmp_assembly_' + str(int(random.random()*1000000))+'.txt'
    with open(n, 'w') as fh:
        print('writing', n)
        fh.write('source\ttarget\tlabel\n')
        for src, tgt in assembly.edges():
            fh.write('{0}\t{1}\t{2}\n'.format(src, tgt, assembly.edge_freq[(src, tgt)]))
    
    for component in digraph_connected_components(assembly):
        # since now we know it's a tree, the assemblies will all be ltd to simple paths
        sources = set()
        sinks = set()
        for node in component:
            if assembly.degree(node) == 0: # ignore singlets
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
                    score += assembly.edge_freq[(path[i], path[i+1])]
                results[s] = max(results.get(s, 0), score)
     
    for seq, score in list(results.items()):
        results[seq] = Contig(seq, score)
        if seq in sequences:
            del results[seq]
    
    # remap the input reads
    for seq in sequences:
        maps_to = {} # contig, score
        for contig in results:
            if len(contig) < len(seq):
                continue
            a = nsb_overhang_align(contig, seq, min_overlap = min_read_mapping_overlap)
            if len(a) != 1:
                continue
            a = a[0]
            score = percent_query_match(a.cigar)
            if score < min_match_quality:
                continue
            maps_to[contig] = a
        for contig, read in maps_to.items():
            results[contig].add_mapped_read(read, len(maps_to.keys()))
            #remapped_sequences[contig] += 1 / len(maps_to.keys())

    return results.values()

def blat_contigs(sequences, **kwargs):
    """
    given a set of contigs, call blat from the commandline and return the results
    """
    ref = kwargs.pop('ref', '/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit')
    chr_to_index = kwargs.pop('chr_to_index')
    min_query_anchor = kwargs.pop('min_query_anchor', 10)
    min_query_coverage = kwargs.pop('min_query_coverage', 0.25)
    max_indel_size = kwargs.pop('max_indel_size', 50)
    min_percent_of_max_score = kwargs.pop('min_percent_of_max_score', 0.75)
    align_from_end_of_contig = kwargs.pop('align_from_end_of_contig', 10)
    min_identity = kwargs.pop('min_identity', 95)
    min_nonoverlap = kwargs.pop('min_nonoverlap', 10)
    blat_options = kwargs.pop('blat_options', 
            ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(min_identity)])
    is_protein = kwargs.pop('is_protein', False)
    
    # write the input sequences to a fasta file
    fasta = tempfile.NamedTemporaryFile(mode='w', delete=False)
    fasta_name = fasta.name
    query_id_mapping = {}
    count = 1
    for seq in sequences:
        n = 'seq{0}'.format(count)
        query_id_mapping[n] = seq
        fasta.write('>' + n + '\n' + seq + '\n')
        count += 1
    fasta.close()
    
    # call the blat subprocess
    psl = tempfile.NamedTemporaryFile(delete=False)
    # will raise subprocess.CalledProcessError if non-zero exit status
    # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    print(["blat", ref, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options)
    subprocess.check_output(["blat", ref, fasta_name, psl.name, '-out=pslx', '-noHead'] + blat_options) 
    psl.close()
    
    header, rows = Blat.read_pslx(psl.name, query_id_mapping)
    
    # split the rows by query id
    rows_by_query = {}
    for row in rows:
        if row['qname'] not in rows_by_query:
            rows_by_query[row['qname']] = []
        rows_by_query[row['qname']].append(row)
    
    reads_by_query = {}
    for query_id, rows in rows_by_query.items():
        query_seq = query_id_mapping[query_id]
        reads = []
        max_score = max([r['score'] for r in rows])

        acceptable_alignments = sum([1 for x in rows if x['score'] >= max_score * min_percent_of_max_score ])
        for rank, row in enumerate(sorted(rows, key = lambda x: x['score'], reverse=True)):
            if row['score'] < max_score * min_percent_of_max_score or row['percent_ident'] < min_identity:
                continue
            try:
                read = Blat.pslx_row_to_pysam(row, chr_to_index = chr_to_index)
                read.set_tag('bs', row['score'], value_type='i')
                read.set_tag('ba', acceptable_alignments, value_type='i')
                read.set_tag('bp', min_percent_of_max_score, value_type='f')
                read.set_tag('br', rank, value_type='i')
                read.set_tag('bi', row['percent_ident'], value_type='f')
                reads.append(read)
            except KeyError as e:
                warnings.warn('warning: reference template name not recognized {0}'.format(e))
        reads_by_query[query_seq] = reads
    return reads_by_query

if __name__ == '__main__':
    import doctest
    doctest.testmod()
