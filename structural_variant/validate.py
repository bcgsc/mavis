import pysam
import itertools
from copy import copy as sys_copy
import re
from structural_variant.constants import *
#from structural_variant.align import *
from structural_variant.align import CigarTools, blat_contigs, assemble, nsb_align
from structural_variant.interval import Interval
import structural_variant.annotate as ann
from structural_variant.breakpoint import Breakpoint, BreakpointPair
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
import datetime
import time

class EvidenceSettings:
    """
    holds all the user input settings associated with evidence gathering
    separate class to allow for easy transfer and sharing of settings
    between evidence objects
    """
    def __init__(self, **kwargs):
        self.read_length = kwargs.pop('read_length', 125)
        self.average_insert_size = kwargs.pop('average_insert_size', 450)
        self.stdev_insert_size = kwargs.pop('stdev_insert_size', 25)
        self.call_error = kwargs.pop('call_error', 10)
        self.min_splits_reads_resolution = kwargs.pop('min_splits_reads_resolution', 3)
        self.convert_chr_to_index = kwargs.pop('convert_chr_to_index', {})
        self.bamfile = kwargs.pop('bamfile')
        self.min_anchor_exact = kwargs.pop('min_anchor_exact', 6)
        self.min_anchor_fuzzy = kwargs.pop('min_anchor_fuzzy', 10) # allow a single event to interrupt the sequence
        self.min_anchor_size = min(self.min_anchor_exact, self.min_anchor_fuzzy)
        self.min_anchor_match = kwargs.pop('min_anchor_match', 0.75) # match / (miss + match)
        self.min_mapping_quality = kwargs.pop('min_mapping_quality', 20)
        self.max_sc_preceeding_anchor = kwargs.pop('max_sc_preceeding_anchor', self.min_anchor_size)
        self.max_reads_limit = kwargs.pop('max_reads_limit', 1000000)
        self.max_anchor_events = kwargs.pop('max_anchor_events', 5)
        self.filter_secondary_alignments = kwargs.pop('filter_secondary_alignments', True)
        self.convert_chr_to_index = kwargs.pop('convert_chr_to_index', {})
        self.convert_index_to_chr = {}
        self.update_chr_to_index(self.convert_chr_to_index)
        self.consensus_req = kwargs.pop('consensus_req', 3)
        
    
    def update_chr_to_index(self, d):
        temp = {}
        for k, v in d.items():
            if v in temp:
                raise AttributeError('indices must be unique', v)
            temp[v] = k
        self.convert_chr_to_index.update(d)
        self.convert_index_to_chr.update(temp)

class Evidence:
    @property
    def window1(self):
        if self._window1 is None:
            return self._window(self.break1)
        else:
            return self._window1
    
    @property
    def window2(self):
        if self._window2 is None:
            return self._window(self.break2)
        else:
            return self._window1
    
    @property
    def break1(self):
        return self.breakpoint_pair.break1
    
    @property
    def break2(self):
        return self.breakpoint_pair.break2

    def _window(self, breakpoint):
        """
        given some input breakpoint uses the current evidence settting to determine an
        appropriate window/range of where one should search for supporting reads
        """
        temp = self.settings.read_length*2 + self.settings.average_insert_size
        start = breakpoint.start - temp - self.settings.call_error - self.settings.read_length - 1
        end = breakpoint.end + temp + self.settings.call_error + self.settings.read_length - 1
        
        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + self.settings.call_error + self.settings.read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - self.settings.call_error - self.settings.read_length - 1
        return Interval(start, end)
    
    def __init__(self, breakpoint_pair, **kwargs):
        """
        @param =average_insert_size \a optional (type: int; default: 450)
        @param =min_anchor_size \a optional (type: int; default: 5)
        @param =min_mapping_quality \a optional (type: int; default: 20)
        @param =read_length \a optional (type: int; default: 125)
        """
        self.settings = EvidenceSettings(**kwargs)
        
        self.labels = kwargs.pop('labels', {})
        self.classification = kwargs.pop('classification', None)
        if self.classification is not None and self.classification not in BreakpointPair.classify(breakpoint_pair):
            raise AttributeError('breakpoint pair improper classification', 
                    BreakpointPair.classify(breakpoint_pair), self.classification)
        
        self.breakpoint_pair = breakpoint_pair
        # split reads are a read that covers at least one breakpoint
        # to avoid duplicating with spanning should try adding as spanning first
        self.split_reads = {
                self.break1: set(), 
                self.break2: set()
                }
        # flanking reads are read pairs that have a mate within a given breakpoint window and 
        # their pair mapped to the opposite breakpoint window
        self.flanking_reads = {
                self.break1: set(), 
                self.break2: set()
                }         

        # spanning reads are reads spanning BOTH breakpoints
        self.spanning_reads = set() 
        # for each breakpoint stores the number of reads that were read from the associated
        # bamfile for the window surrounding the breakpoint
        self.read_counts = {}
        self._window1 = None
        self._window2 = None
    
    def supporting_reads(self):
        """
        convenience method to return all flanking, split and spanning reads associated with an evidence object
        """
        result = set()
        for s in self.flanking_reads.values(): 
            result.update(s)
        for s in self.split_reads.values():
            result.update(s)
        result.update(self.spanning_reads)
        return result

    def linking_evidence(self):
        # can link a pair of breakpoints if the read prefix is in the evidence for both
        # or if the mate id is in the other breakpoint
        first_prefixes = set()
        second_prefixes = set()
        for read in self.split_reads[self.break1]:
            prefix = read.query_name
            if SUFFIX_DELIM in read.query_name:
                prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
            first_prefixes.add(prefix)
        for read in self.split_reads[self.break2]:
            prefix = read.query_name
            if SUFFIX_DELIM in read.query_name:
                prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
            second_prefixes.add(prefix)
        return len(first_prefixes.intersection(second_prefixes))

    def add_spanning_read(self, read): # TODO
        """
        spanning read: a read covering BOTH breakpoints
        
        this is only applicable to small events
        """
        # check that the read fully covers BOTH breakpoints
        read_start = read.reference_start + 1 - self.settings.call_error # adjust b/c pysam is 0-indexed
        read_end = read.reference_end + self.settings.call_error        # don't adjust b/c pysam is 0-indexed but end coord are one past
        if self.break1.start >= read_start  and self.break1.end <= read_end  \
                and self.break2.start >= read_start and self.break2.end <= read_end:
                    pass # now check if this supports the putative event types
                    raise NotImplementedError('have not added support for indels yet')
        else:
            raise UserWarning('this does not cover/span both breakpoints and cannot be added as spanning evidence')
    
    @staticmethod
    def read_pair_type(read):
        # check if the read pair is in the expected orientation
        """
        assumptions based on illumina pairs: only 4 possible combinations
        ++++> <---- is LR same-strand
        ++++> ++++> is LL opposite
        <---- <---- is RR opposite
        <---- ++++> is RL same-strand
        """
        reverse = False
        if read.reference_id == read.next_reference_id and read.reference_start > read.next_reference_start:
            reverse = True

        if not read.is_reverse and read.mate_is_reverse: # LR
            return READ_PAIR_TYPE.RL if reverse else READ_PAIR_TYPE.LR
        elif not read.is_reverse and not read.mate_is_reverse: # LL opp
            return READ_PAIR_TYPE.LL
        elif read.is_reverse and read.mate_is_reverse: # RR opp
            return READ_PAIR_TYPE.RR
        elif read.is_reverse and not read.mate_is_reverse: # RL
            return READ_PAIR_TYPE.LR if reverse else READ_PAIR_TYPE.RL
        else:
            raise NotImplementedError('unexpected orientation for pair')
    
    def add_flanking_read(self, read):
        """
        checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
        this event
        """
        if read.is_unmapped or read.mate_is_unmapped:
            raise UserWarning('input read (and its mate) must be mapped')
        
        
        # filter by putative event classifications
        classifications = [self.classification] if self.classification else BreakpointPair.classify(self.breakpoint_pair)
        classifications = sorted(classifications)
        
        insert_size = abs(read.template_length)
        
        if classifications == sorted([SVTYPE.DEL, SVTYPE.INS]):
            if insert_size <= self.settings.average_insert_size + self.settings.stdev_insert_size \
                    and insert_size >= self.settings.average_insert_size - self.settings.stdev_insert_size:
                        raise UserWarning('insert size is not abnormal. does not support del/ins', insert_size)
        elif classifications == [SVTYPE.DEL]:
            if insert_size <= self.settings.average_insert_size + self.settings.stdev_insert_size:
                raise UserWarning('insert size is smaller than expected for a deletion type event', insert_size)
        elif classifications == [SVTYPE.INS]:
            if insert_size >= self.settings.average_insert_size - self.settings.stdev_insert_size:
                raise UserWarning('insert size is larger than expected for an insertion type event', insert_size)
        
        # check if the read orientation makes sense with the event type
        rt = Evidence.read_pair_type(read)
        if rt == READ_PAIR_TYPE.LR:
            if len(set([SVTYPE.INS, SVTYPE.DEL, SVTYPE.TRANS]).intersection(set(classifications))) == 0:
                raise UserWarning('read pair orientation does not match event type', rt, classifications)
        elif rt == READ_PAIR_TYPE.RL:
            if SVTYPE.DUP not in classifications and SVTYPE.TRANS not in classifications:
                raise UserWarning('read pair orientation does not match event type', rt, classifications)
        else:
            if SVTYPE.INV not in classifications and SVTYPE.ITRANS not in classifications:
                raise UserWarning('read pair orientation does not match event type', rt, classifications)
        
        # check if this read falls in the first breakpoint window
        w1 = self.window1
        w1 = (w1[0] - 1, w1[1] - 1) # correct for psyam using 0-based coordinates
        w2 = self.window2
        w2 = (w2[0] - 1, w2[1] - 1) # correct for psyam using 0-based coordinates
        
        if read.reference_start >= w1[0] and read.reference_end <= w1[1] \
                and read.reference_id == self.settings.convert_chr_to_index[self.break1.chr] \
                and read.next_reference_start >= w2[0] and read.next_reference_start <= w2[1] \
                and read.next_reference_id == self.settings.convert_chr_to_index[self.break2.chr]:
            # current read falls in the first breakpoint window, mate in the second
            self.flanking_reads[self.break1].add(read)
        elif read.reference_start >= w2[0] and read.reference_end <= w2[1] \
                and self.settings.convert_chr_to_index[self.break2.chr] == read.reference_id \
                and read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                and self.settings.convert_chr_to_index[self.break1.chr] == read.next_reference_id:
            # current read falls in the second breakpoint window, mate in the first
            self.flanking_reads[self.break2].add(read)
        else:
            raise UserWarning('does not map to the expected regions. does not support the current breakpoint pair')
    
    def add_split_read(self, read, first_breakpoint=True):
        """
        adds a split read if it passes the criteria filters and raises a warning if it does not
        """
        breakpoint = self.break1 if first_breakpoint else self.break2
        window = self.window1 if first_breakpoint else self.window2
        opposite_breakpoint = self.break2 if first_breakpoint else self.break1
        opposite_window = self.window2 if first_breakpoint else self.window1
        
        if read.cigar[0][0] != CIGAR.S and read.cigar[-1][0] != CIGAR.S:
            raise UserWarning('split read is not softclipped')
        elif breakpoint.orient == ORIENT.LEFT and read.cigar[-1][0] != CIGAR.S:
            raise UserWarning('split read is not softclipped')
        elif breakpoint.orient == ORIENT.RIGHT and read.cigar[0][0] != CIGAR.S:
            raise UserWarning('split read is not softclipped')

        # the first breakpoint of a BreakpointPair is always the lower breakpoint
        # if this is being added to the second breakpoint then we'll need to check if the 
        # read soft-clipping needs to be adjusted
        
        # need to do this after shifting? assume shifting amount is insignificant
        s, t = (window[0], window[1])
        s -= 1 # correct for pysam using 0-based coordinates
        t -= 1 # correct for pysam using 0-based coordinates
        
        if read.reference_start > t or read.reference_end < s \
                or self.settings.convert_index_to_chr[read.reference_id] != breakpoint.chr:
            raise UserWarning('read does not map within the breakpoint evidence window')
        if self.breakpoint_pair.stranded:
            if ( read.is_reverse and breakpoint.strand == STRAND.POS ) \
                    or ( not read.is_reverse and breakpoint.strand == STRAND.NEG ):
                        raise UserWarning('split read not on the appropriate strand')
        primary = ''
        clipped = ''
        if breakpoint.orient == ORIENT.LEFT:
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            clipped = read.query_sequence[read.query_alignment_end:] # end is exclusive in pysam
        elif breakpoint.orient == ORIENT.RIGHT:
            clipped = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
        else:
            raise AttributeError('cannot assign split reads to a breakpoint where the orientation has not been '
                    'specified')
        if len(primary) < self.settings.min_anchor_size or len(clipped) < self.settings.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria', primary, clipped)
        elif len(read.query_sequence) - (read.query_alignment_end + 2) < self.settings.min_anchor_size \
                and (read.query_alignment_start + 1) < self.settings.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria')
        elif len(primary) < self.settings.min_anchor_size or len(clipped) < self.settings.min_anchor_size:
            raise UserWarning('split read does not meet the minimum anchor criteria')
        
        read = sys_copy(read)
        # recalculate the read cigar string to ensure M is replaced with = or X
        c = CigarTools.recompute_cigar_mismatch(
                read, 
                ann.HUMAN_REFERENCE_GENOME[self.settings.convert_index_to_chr[read.reference_id]].seq
                )
        prefix = 0
        try:
            c, prefix = CigarTools.extend_softclipping(c, self.settings.min_anchor_size)
        except AttributeError:
            pass
        read.cigar = c
        read.reference_start = read.reference_start + prefix
        # data quality filters

        if CigarTools.alignment_matches(read.cigar) >= 10 \
                and CigarTools.match_percent(read.cigar) < self.settings.min_anchor_match:
            raise UserWarning('alignment of too poor quality') #pass #print('bad quality read', read.cigar, read.query_name)
        if CigarTools.longest_exact_match(read.cigar) < self.settings.min_anchor_exact \
                and CigarTools.longest_fuzzy_match(read.cigar, 1) < self.settings.min_anchor_fuzzy:
            raise UserWarning('alignment of too poor quality')
        else:
            self.split_reads[breakpoint].add(read)
            
        # try mapping the soft-clipped portion to the other breakpoint
        w = (opposite_window[0], opposite_window[1])
        opposite_breakpoint_ref = ann.HUMAN_REFERENCE_GENOME[opposite_breakpoint.chr].seq[w[0] - 1: w[1]]
        
        putative_alignments = None
        revcomp_primary = str(Seq(primary, DNA_ALPHABET).reverse_complement())

        if not self.breakpoint_pair.opposing_strands:
            sc_align = nsb_align(opposite_breakpoint_ref, read.query_sequence)
            
            for a in sc_align:
                a.flag = read.flag
            putative_alignments = sc_align
        else:
            # check if the revcomp will give us sc_align better alignment
            revcomp_sc_align = Seq(read.query_sequence, DNA_ALPHABET).reverse_complement()
            revcomp_sc_align = nsb_align(opposite_breakpoint_ref, str(revcomp_sc_align))
            
            for a in revcomp_sc_align:
                a.flag = read.flag ^ 16
            putative_alignments = revcomp_sc_align
        
        scores = []
        
        for a in putative_alignments:
            a.flag = a.flag ^ 64 ^ 128
            
            # add information from the original read
            a.reference_start = w[0] - 1 + a.reference_start  
            a.reference_id = self.settings.convert_chr_to_index[opposite_breakpoint.chr]
            a.query_name = read.query_name + SUFFIX_DELIM + 'clipped-realign'
            a.next_reference_start = read.next_reference_start
            a.next_reference_id = read.next_reference_id
            a.mapping_quality = NA_MAPPING_QUALITY
            try:
                cigar, offset = CigarTools.extend_softclipping(a.cigar, self.settings.min_anchor_size)
                a.cigar = cigar
                a.reference_start = a.reference_start + offset
            except AttributeError:
                # if the matches section is too small you can't extend the softclipping
                pass
            s = CigarTools.score(a.cigar)
            
            if (CigarTools.alignment_matches(a.cigar) >= 10 \
                    and CigarTools.match_percent(a.cigar) < self.settings.min_anchor_match):
                continue 
            if CigarTools.longest_exact_match(a.cigar) < self.settings.min_anchor_exact \
                    and CigarTools.longest_fuzzy_match(a.cigar, 1) < self.settings.min_anchor_fuzzy:
                continue
            if opposite_breakpoint.orient == ORIENT.LEFT:
                if a.cigar[0][0] == CIGAR.S and a.cigar[0][1] > self.settings.max_sc_preceeding_anchor:
                    continue
            elif opposite_breakpoint.orient == ORIENT.RIGHT:
                if a.cigar[-1][0] == CIGAR.S and a.cigar[-1][1] > self.settings.max_sc_preceeding_anchor:
                    continue
            scores.append((s, CigarTools.match_percent(a.cigar), a))
        
        scores = sorted(scores, reverse = True) if len(scores) > 0 else []

        if len(scores) > 1:
            if scores[0][0] != scores[1][0] and scores[0][1] != scores[1][1]: # not multimap
                clipped = scores[0][2]
                self.split_reads[opposite_breakpoint].add(clipped)
        elif len(scores) == 1:
            clipped = scores[0][2]
            self.split_reads[opposite_breakpoint].add(clipped)
    
    def resolve_breakpoints(self):
        """
        use split read evidence to resolve bp-level calls for breakpoint pairs (where possible)
        if a bp level call is not possible for one of the breakpoints then returns None
        if no breakpoints can be resolved returns an the original event only with NO split read evidence
        also sets the SV type call if multiple are input
        """
        pos1 = {}
        pos2 = {}
        
        for breakpoint, d in [(self.break1, pos1), (self.break2, pos2)]:
            for read in self.split_reads[breakpoint]:
                pos = breakpoint_pos(read, breakpoint.orient)
                if pos not in d:
                    d[pos] = []
                d[pos].append(read)
            putative_positions = list(d.keys())
            print('; '.join(['{0}[{1}]'.format(p, len(d[p])) for p in putative_positions]))
            for pos in putative_positions:
                if len(d[pos]) < self.settings.min_splits_reads_resolution:
                    del d[pos]

        linked_pairings = []
        print('putative breakpoint positions:', pos1.keys(), pos2.keys())
        # now pair up the breakpoints with their putative partners
        for first, second in itertools.product(pos1, pos2):
            # can link a pair of breakpoints if the read prefix is in the evidence for both
            # or if the mate id is in the other breakpoint
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                bp.break1.pos = Interval(first)
                bp.break2.pos = Interval(second)
                e = Evidence(bp, bamfile=self.settings.bamfile)
                e.settings = self.settings
                e.classification = cls
                e.split_reads[e.break1] = set(pos1[first])
                e.split_reads[e.break2] = set(pos2[second])
            
                if e.linking_evidence() == 0:
                    continue
                for read in self.flanking_reads[self.break1].union(self.flanking_reads[self.break2]):
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass
                
                linked_pairings.append(e)
        
        psuedo_pairings = set()
        print('pos1:', pos1.keys())
        print('linked pairings', [ temp.break1.start for temp in linked_pairings ])
        print('pos2:', pos2.keys())
        print('linked pairings', [ temp.break2.start for temp in linked_pairings ])
        for first in pos1:
            if first not in [ temp.break1.start for temp in linked_pairings ]:
                for second in pos2:
                    psuedo_pairings.add((first, second))
        for second in pos2:
            if second not in [ temp.break2.start for temp in linked_pairings ]:
                for first in pos1:
                    psuedo_pairings.add((first, second))
        
        for first, second in psuedo_pairings:
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                bp.break1.pos = Interval(first)
                bp.break2.pos = Interval(second)
                e = Evidence(bp, bamfile=self.settings.bamfile)
                e.settings = self.settings
                e.classification = cls
                e.split_reads[e.break1] = set(pos1[first])
                e.split_reads[e.break2] = set(pos2[second])

                for read in self.flanking_reads[self.break1].union(self.flanking_reads[self.break2]):
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass
                
                linked_pairings.append(e)
        
        if len(linked_pairings) == 0:
            assert(len(pos1.keys()) == 0 or len(pos2.keys() == 0))
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                e = Evidence(bp, bamfile=self.settings.bamfile)
                e.settings = self.settings
                e.classification = cls
                
                for read in self.flanking_reads[self.break1].union(self.flanking_reads[self.break2]):
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass
                linked_pairings.append(e)
        return linked_pairings

    def load_evidence(self, open_bam=None):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
        bamfile = None
        if open_bam is None:
            bamfile = pysam.AlignmentFile(self.settings.bamfile, 'rb')
        else:
            bamfile = open_bam
        
        convert_chr_to_index = {}
        for name in bamfile.references:
            convert_chr_to_index[name] = bamfile.gettid(name)
        
        self.settings.update_chr_to_index(convert_chr_to_index)
        # TODO transcriptome window gathering
        
        count = 0
        for read in bamfile.fetch(
                '{0}'.format(self.break1.chr), 
                self.window1[0],
                self.window1[1]):
            
            count += 1
            if count > self.settings.max_reads_limit:
                break
            if read.is_unmapped \
                    or read.mapping_quality < self.settings.min_mapping_quality \
                    or (read.is_secondary and self.settings.filter_secondary_alignments):
                continue
            try:
                self.add_spanning_read(read)
            except (UserWarning, NotImplementedError): # TODO: ADD SUPPORT FOR INDELS
                try:
                    self.add_split_read(read)
                except UserWarning:
                    pass
            try:
                self.add_flanking_read(read)
            except UserWarning:
                pass
        self.read_counts[self.break1] = count
        count = 0
        for read in bamfile.fetch(
                '{0}'.format(self.break2.chr), 
                self.window2[0],
                self.window2[1]):
            
            count += 1
            if count > self.settings.max_reads_limit:
                break
            if read.is_unmapped \
                    or read.mapping_quality < self.settings.min_mapping_quality \
                    or (read.is_secondary and self.settings.filter_secondary_alignments):
                continue
            try:
                self.add_split_read(read, False)
            except UserWarning:
                pass
            try:
                self.add_flanking_read(read)
            except UserWarning:
                pass
        self.read_counts[self.break2] = count
        
        if open_bam is None:
            bamfile.close()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
