"""This module is primarily responsible for collecting read evidence from bams
"""

import itertools
from copy import copy as sys_copy
from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.align import CigarTools, nsb_align, breakpoint_pos
from structural_variant.interval import Interval
from structural_variant.annotate import overlapping_transcripts
from structural_variant.breakpoint import BreakpointPair
from Bio.Seq import Seq
import tools.profile_bam as profile_bam
from functools import partial

CALL_METHOD = Vocab(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads')


class Event:
    def __init__(
            self,
            evidence,
            breakpoint_pair,
            classification,
            call_method,
            contig=None,
            contig_alignment=None,
            flanking_reads=set()
    ):
        self.breakpoint_pair = breakpoint_pair
        self.classification = classification
        self.evidence = evidence
        self.flanking_reads = flanking_reads
        self.call_method = CALL_METHOD.enforce(call_method)
        self.contig = contig
        self.contig_alignment = contig_alignment

    def flanking_read_stdev(self):
        ihist = {}
        for read in self.flanking_reads:
            isize = abs(read.template_length)
            ihist[isize] = ihist.get(isize, 0) + 1
        
        if len(self.flanking_reads) > 1:
            median = profile_bam.histogram_median(ihist)
            stdev = math.sqrt(profile_bam.histogram_stderr(ihist, median))
            return stdev
        else:
            raise AttributeError('cannot calculate the stdev of less than 2')


class EvidenceSettings:
    """
    holds all the user input settings associated with evidence gathering
    separate class to allow for easy transfer and sharing of settings
    between evidence objects
    """

    def __init__(
            self,
            read_length=125,
            median_insert_size=380,
            mm_isize_error=100,
            call_error=10,
            min_splits_reads_resolution=3,
            min_anchor_exact=6,
            min_anchor_fuzzy=10,
            min_anchor_match=0.9,
            min_mapping_quality=20,
            fetch_reads_limit=10000,
            fetch_reads_bins=3,
            filter_secondary_alignments=True,
            consensus_req=3,
            sc_extension_stop=5,
            max_sc_preceeding_anchor=None):
        """
        Args:
            read_length (int, default=125): length of individual reads
            median_insert_size (int, default=380): expected average insert size for paired-end reads
            mm_isize_error (int, default=100): the expected deviation in the insert size
            call_error (int, default=10): buffer for calculating the evidence window
            min_splits_reads_resolution (int, default=3):
                minimum number of reads required to call the same breakpoint for it to be a valid breakpoint
            min_anchor_exact (int, default=6):
                minimum number of consecutive exact matches to satisfy the anchor constraint
            min_anchor_fuzzy (int, default=10):
                minimum number of consecutive exact matches (allowing one mismatch/indel) to satisfy the anchor
                constraint
            min_mapping_quality (int, default=20):
                mapping quality to filter on
            fetch_reads_limit (int, default=100000):
                maximum number of reads to loop over for a given event
            filter_secondary_alignments (bool, default=True):
                don't use secondary alignments when reading evidence from the bam file
            sc_extension_stop (int, default=5):
                when extending softclipped, stop given this number of exact consecutive matches
        """
        self.read_length = read_length
        self.median_insert_size = median_insert_size
        self.mm_isize_error = mm_isize_error
        self.call_error = call_error
        self.min_splits_reads_resolution = min_splits_reads_resolution
        self.min_anchor_exact = min_anchor_exact
        self.min_anchor_fuzzy = min_anchor_fuzzy
        self.min_anchor_size = min(self.min_anchor_exact, self.min_anchor_fuzzy)
        self.min_anchor_match = min_anchor_match
        self.min_mapping_quality = min_mapping_quality
        if max_sc_preceeding_anchor is None:
            self.max_sc_preceeding_anchor = self.min_anchor_size
        else:
            self.max_sc_preceeding_anchor = max_sc_preceeding_anchor
        self.fetch_reads_limit = fetch_reads_limit
        self.fetch_reads_bins = fetch_reads_bins
        self.filter_secondary_alignments = filter_secondary_alignments
        self.consensus_req = consensus_req
        self.sc_extension_stop = sc_extension_stop


class Evidence:
    """
    """
    @property
    def window1(self):
        return self.windows[self.break1]

    @property
    def window2(self):
        return self.windows[self.break2]

    @property
    def break1(self):
        return self.breakpoint_pair.break1

    @property
    def break2(self):
        return self.breakpoint_pair.break2

    @classmethod
    def generate_window(cls, breakpoint, read_length, median_insert_size, call_error, mm_isize_error):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            read_length (int): the read length
            median_insert_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            mm_isize_error:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        fragment_size = read_length * 2 + median_insert_size + mm_isize_error * 2
        start = breakpoint.start - fragment_size - call_error
        end = breakpoint.end + fragment_size + call_error

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + read_length
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - read_length
        return Interval(start, end)

    @classmethod
    def generate_transcriptome_window(cls, breakpoint, annotations, read_length, median_insert_size, call_error,
                                      mm_isize_error):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            annotations (Dict[str,List[Gene]]): the set of reference annotations: genes, transcripts, etc
            read_length (int): the read length
            median_insert_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            mm_isize_error:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        transcripts = overlapping_transcripts(annotations, breakpoint)
        window = cls.generate_window(breakpoint, read_length, median_insert_size, call_error, mm_isize_error)

        tgt_left = breakpoint.start - window.start + 1
        tgt_right = window.end - breakpoint.end + 1

        if len(transcripts) == 0:  # case 1. no overlapping transcripts
            return window

        for t in transcripts:
            current_length = 0
            exons = t.get_exons()
            current_interval = Interval(breakpoint.start, breakpoint.end)

            # first going left
            epos, in_prev_intron = Interval.position_in_range(exons, (breakpoint.start, breakpoint.start))
            if epos == 0 and in_prev_intron:
                continue
            elif in_prev_intron:
                epos -= 1
                current_length += breakpoint.start - exons[epos].end + 1
                current_interval.start = exons[epos].end - 1
                if current_length >= tgt_left:
                    continue
            while epos >= 0:
                if current_length + len(exons[epos]) >= tgt_left:
                    eshift = tgt_left - current_length
                    current_interval.start = exons[epos].end - eshift
                    current_length += eshift
                    assert(current_length == tgt_left)
                    break
                else:
                    current_length += len(exons[epos])
                    current_length.start = exons[epos].start
                    epos -= 1
            if current_length < tgt_left:
                assert(epos == -1)
                eshift = tgt_left - current_length
                current_interval.start = exons[0].start - eshift

            current_length = 0
            # next going right
            epos, in_prev_intron = Interval.position_in_range(exons, (breakpoint.end, breakpoint.end))
            if epos == len(exons):  # after the last exon
                continue
            elif in_prev_intron:
                current_length += exons[epos].start - breakpoint.end + 1
                current_interval.end = exons[epos].start - 1
                if current_length >= tgt_right:
                    continue
            while epos < len(exons):
                if current_length + len(exons[epos]) >= tgt_right:
                    eshift = tgt_right - current_length
                    current_interval.end = exons[epos].start + eshift
                    current_length += eshift
                    assert(current_length == tgt_right)
                    break
                else:
                    current_length += len(exons[epos])
                    current_interval.end = exons[epos].end
                    epos += 1
            if current_length < tgt_right:
                assert(epos == len(exons))
                eshift = tgt_right - current_length
                current_interval.end = exons[-1].end + eshift
            window = window | current_interval
        return window

    def __init__(
            self,
            breakpoint_pair,
            bam_cache,
            human_reference_genome,
            labels={},
            classification=None,
            protocol=PROTOCOL.GENOME,
            annotations={},
            **kwargs):
        """
        Args:
            breakpoint_pair (BreakpointPair): the breakpoint pair to collect evidence for
            bam_cache (BamCache): the bam cache (and assc file) to collect evidence from
            human_reference_genome (SeqIO.iterator[SeqIO.SeqRecord]): the human reference genome read as fasta
            labels (Dict, optional): a dictionary of labels to associate with the evidence object
            classification (SVTYPE, optional): the event type
            protocol (PROTOCOL, default=PROTOCOL.GENOME): genome or transcriptome
            **kwargs: named arguments to be passed to EvidenceSettings
        """
        self.settings = EvidenceSettings(**kwargs)
        self.bam_cache = bam_cache
        self.labels = labels
        self.classification = classification
        self.protocol = PROTOCOL.enforce(protocol)
        self.human_reference_genome = human_reference_genome
        self.annotations = annotations

        if self.classification is not None and self.classification not in BreakpointPair.classify(breakpoint_pair):
            raise AttributeError('breakpoint pair improper classification',
                                 BreakpointPair.classify(breakpoint_pair), self.classification)

        if not self.annotations and self.protocol == PROTOCOL.TRANS:
            raise AttributeError('must specify the reference annotations for transcriptome evidence objects')

        self.breakpoint_pair = breakpoint_pair

        if self.break1.orient == ORIENT.NS or self.break2.orient == ORIENT.NS:
            raise AttributeError(
                'input breakpoint pair must specify strand and orientation. Cannot be \'not specified'
                '\' for evidence gathering')

        # split reads are a read that covers at least one breakpoint
        # to avoid duplicating with spanning should try adding as spanning
        # first
        self.split_reads = {
            self.break1: set(),
            self.break2: set()
        }
        # flanking reads are read pairs that have a mate within a given breakpoint window and
        # their pair mapped to the opposite breakpoint window
        self.flanking_reads = set()

        # spanning reads are reads spanning BOTH breakpoints
        self.spanning_reads = set()
        # for each breakpoint stores the number of reads that were read from the associated
        # bamfile for the window surrounding the breakpoint
        self.read_counts = {}
        self.windows = {}
        self.contigs = []

        if self.protocol == PROTOCOL.GENOME:
            self.windows[self.break1] = Evidence.generate_window(
                self.break1,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                mm_isize_error=self.settings.mm_isize_error
            )
            self.windows[self.break2] = Evidence.generate_window(
                self.break2,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                mm_isize_error=self.settings.mm_isize_error
            )
        elif self.protocol == PROTOCOL.TRANS:
            self.windows[self.break1] = Evidence.generate_transcriptome_window(
                self.break1,
                self.annotations,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                mm_isize_error=self.settings.mm_isize_error
            )
            self.windows[self.break2] = Evidence.generate_transcriptome_window(
                self.break2,
                self.annotations,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                mm_isize_error=self.settings.mm_isize_error
            )
        else:
            raise AttributeError('invalid protocol', self.protocol)
        self.half_mapped = {
            self.break1: set(),
            self.break2: set()
        }

    def supporting_reads(self):
        """
        convenience method to return all flanking, split and spanning reads associated with an evidence object
        """
        result = set()
        result.update(self.flanking_reads)
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

    def add_spanning_read(self, read):  # TODO
        """
        spanning read: a read covering BOTH breakpoints

        this is only applicable to small events
        """
        # check that the read fully covers BOTH breakpoints
        read_start = read.reference_start + 1 - \
            self.settings.call_error  # adjust b/c pysam is 0-indexed
        # don't adjust b/c pysam is 0-indexed but end coord are one past
        read_end = read.reference_end + self.settings.call_error
        if self.break1.start >= read_start and self.break1.end <= read_end  \
                and self.break2.start >= read_start and self.break2.end <= read_end:
            pass  # now check if this supports the putative event types
            raise NotImplementedError('have not added support for indels yet')
        else:
            raise UserWarning(
                'this does not cover/span both breakpoints and cannot be added as spanning evidence')

    @staticmethod
    def read_pair_type(read):
        # check if the read pair is in the expected orientation
        """
        assumptions based on illumina pairs: only 4 possible combinations

        ::

            ++++> <---- is LR same-strand
            ++++> ++++> is LL opposite
            <---- <---- is RR opposite
            <---- ++++> is RL same-strand
        """
        reverse = False
        if read.reference_id == read.next_reference_id and read.reference_start > read.next_reference_start:
            reverse = True

        if not read.is_reverse and read.mate_is_reverse:  # LR
            return READ_PAIR_TYPE.RL if reverse else READ_PAIR_TYPE.LR
        elif not read.is_reverse and not read.mate_is_reverse:  # LL opp
            return READ_PAIR_TYPE.LL
        elif read.is_reverse and read.mate_is_reverse:  # RR opp
            return READ_PAIR_TYPE.RR
        elif read.is_reverse and not read.mate_is_reverse:  # RL
            return READ_PAIR_TYPE.LR if reverse else READ_PAIR_TYPE.RL
        else:
            raise NotImplementedError('unexpected orientation for pair')

    def add_flanking_read(self, read):
        """
        checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
        this event

        Args:
            read (pysam.AlignedSegment): the read to add
        Raises:
            UserWarning: the read does not support this event or does not pass quality filters
        """
        if read.is_unmapped or read.mate_is_unmapped:
            raise UserWarning('input read (and its mate) must be mapped')

        # filter by putative event classifications
        classifications = [self.classification] if self.classification else BreakpointPair.classify(
            self.breakpoint_pair)
        classifications = sorted(classifications)

        insert_size = abs(read.template_length)

        if classifications == sorted([SVTYPE.DEL, SVTYPE.INS]):
            if insert_size <= self.settings.median_insert_size + self.settings.mm_isize_error * 2 \
                    and insert_size >= self.settings.median_insert_size - self.settings.mm_isize_error * 2:
                raise UserWarning(
                    'insert size is not abnormal. does not support del/ins', insert_size)
        elif classifications == [SVTYPE.DEL]:
            if insert_size <= self.settings.median_insert_size + self.settings.mm_isize_error * 2:
                raise UserWarning(
                    'insert size is smaller than expected for a deletion type event', insert_size)
        elif classifications == [SVTYPE.INS]:
            if insert_size >= self.settings.median_insert_size - self.settings.mm_isize_error * 2:
                raise UserWarning(
                    'insert size is larger than expected for an insertion type event', insert_size)

        # check if the read orientation makes sense with the event type
        rt = Evidence.read_pair_type(read)
        if rt == READ_PAIR_TYPE.LR:
            if len(set([SVTYPE.INS, SVTYPE.DEL, SVTYPE.TRANS]).intersection(set(classifications))) == 0:
                raise UserWarning(
                    'read pair orientation does not match event type', rt, classifications)
        elif rt == READ_PAIR_TYPE.RL:
            if SVTYPE.DUP not in classifications and SVTYPE.TRANS not in classifications:
                raise UserWarning(
                    'read pair orientation does not match event type', rt, classifications)
        else:
            if SVTYPE.INV not in classifications and SVTYPE.ITRANS not in classifications:
                raise UserWarning(
                    'read pair orientation does not match event type', rt, classifications)

        # check if this read falls in the first breakpoint window
        w1 = self.window1
        # correct for psyam using 0-based coordinates
        w1 = (w1[0] - 1, w1[1] - 1)
        w2 = self.window2
        # correct for psyam using 0-based coordinates
        w2 = (w2[0] - 1, w2[1] - 1)

        if read.reference_start >= w1[0] and read.reference_end <= w1[1] \
                and read.reference_id == self.bam_cache.reference_id(self.break1.chr) \
                and read.next_reference_start >= w2[0] and read.next_reference_start <= w2[1] \
                and read.next_reference_id == self.bam_cache.reference_id(self.break2.chr) \
                and (not self.breakpoint_pair.stranded or not read.is_read1
                    or read.is_reverse == (self.break1.strand == STRAND.NEG)):
            # current read falls in the first breakpoint window, mate in the
            # second
            self.flanking_reads.add(read)
        elif read.reference_start >= w2[0] and read.reference_end <= w2[1] \
                and self.bam_cache.reference_id(self.break2.chr) == read.reference_id \
                and read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                and self.bam_cache.reference_id(self.break1.chr) == read.next_reference_id \
                and (not self.breakpoint_pair.stranded or not read.is_read2
                    or read.is_reverse != (self.break2.strand == STRAND.NEG)):
            # current read falls in the second breakpoint window, mate in the
            # first
            self.flanking_reads.add(read)
        else:
            raise UserWarning(
                'does not map to the expected regions. does not support the current breakpoint pair')

    def add_split_read(self, read, first_breakpoint):
        """
        adds a split read if it passes the criteria filters and raises a warning if it does not

        Args:
            read (pysam.AlignedSegment): the read to add
            first_breakpoint (boolean): add to the first breakpoint (or second if false)
        Raises:
            UserWarning: the read does not support this breakpoint or does not pass quality filters
            AttributeError: orientation wasn't specified for the breakpoint
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

        # need to do this after shifting? assume shifting amount is
        # insignificant
        s, t = (window[0], window[1])
        s -= 1  # correct for pysam using 0-based coordinates
        t -= 1  # correct for pysam using 0-based coordinates

        if read.reference_start > t or read.reference_end < s \
                or self.bam_cache.chr(read) != breakpoint.chr:
            raise UserWarning(
                'read does not map within the breakpoint evidence window')
        if self.breakpoint_pair.stranded:
            if (read.is_reverse and breakpoint.strand == STRAND.POS) \
                    or (not read.is_reverse and breakpoint.strand == STRAND.NEG):
                raise UserWarning('split read not on the appropriate strand')
        primary = ''
        clipped = ''
        if breakpoint.orient == ORIENT.LEFT:
            primary = read.query_sequence[
                read.query_alignment_start:read.query_alignment_end]
            # end is exclusive in pysam
            clipped = read.query_sequence[read.query_alignment_end:]
        elif breakpoint.orient == ORIENT.RIGHT:
            clipped = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[
                read.query_alignment_start:read.query_alignment_end]
        else:
            raise AttributeError('cannot assign split reads to a breakpoint where the orientation has not been '
                                 'specified')
        if len(primary) < self.settings.min_anchor_size or len(clipped) < self.settings.min_anchor_size:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria', primary, clipped)
        elif len(read.query_sequence) - (read.query_alignment_end + 2) < self.settings.min_anchor_size \
                and (read.query_alignment_start + 1) < self.settings.min_anchor_size:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria')
        elif len(primary) < self.settings.min_anchor_size or len(clipped) < self.settings.min_anchor_size:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria')

        read = sys_copy(read)
        # recalculate the read cigar string to ensure M is replaced with = or X
        c = CigarTools.recompute_cigar_mismatch(
            read,
            self.human_reference_genome[self.bam_cache.chr(read)].seq
        )
        prefix = 0
        try:
            c, prefix = CigarTools.extend_softclipping(
                c, self.settings.sc_extension_stop)
        except AttributeError:
            pass
        read.cigar = c
        read.reference_start = read.reference_start + prefix
        # data quality filters

        if CigarTools.alignment_matches(read.cigar) >= 10 \
                and CigarTools.match_percent(read.cigar) < self.settings.min_anchor_match:
            # pass #print('bad quality read', read.cigar, read.query_name)
            raise UserWarning('alignment of too poor quality')
        if CigarTools.longest_exact_match(read.cigar) < self.settings.min_anchor_exact \
                and CigarTools.longest_fuzzy_match(read.cigar, 1) < self.settings.min_anchor_fuzzy:
            raise UserWarning('alignment of too poor quality')
        else:
            self.split_reads[breakpoint].add(read)

        # try mapping the soft-clipped portion to the other breakpoint
        w = (opposite_window[0], opposite_window[1])
        opposite_breakpoint_ref = self.human_reference_genome[
            opposite_breakpoint.chr].seq[w[0] - 1: w[1]]

        putative_alignments = None

        if not self.breakpoint_pair.opposing_strands:
            sc_align = nsb_align(opposite_breakpoint_ref, read.query_sequence)

            for a in sc_align:
                a.flag = read.flag
            putative_alignments = sc_align
        else:
            # check if the revcomp will give us sc_align better alignment
            revcomp_sc_align = Seq(read.query_sequence,
                                   DNA_ALPHABET).reverse_complement()
            revcomp_sc_align = nsb_align(
                opposite_breakpoint_ref, str(revcomp_sc_align))

            for a in revcomp_sc_align:
                a.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
            putative_alignments = revcomp_sc_align

        scores = []

        for a in putative_alignments:
            # a.flag = a.flag ^ 64 ^ 128
            a.flag = a.flag | PYSAM_READ_FLAGS.SECONDARY
            # add information from the original read
            a.reference_start = w[0] - 1 + a.reference_start
            a.reference_id = self.bam_cache.reference_id(opposite_breakpoint.chr)
            # a.query_name = read.query_name + SUFFIX_DELIM + 'clipped-realign'
            a.query_name = read.query_name
            a.set_tag('cr', 1, value_type='i')
            a.next_reference_start = read.next_reference_start
            a.next_reference_id = read.next_reference_id
            a.mapping_quality = NA_MAPPING_QUALITY
            try:
                cigar, offset = CigarTools.extend_softclipping(
                    a.cigar, self.settings.sc_extension_stop)
                a.cigar = cigar
                a.reference_start = a.reference_start + offset
            except AttributeError:
                # if the matches section is too small you can't extend the
                # softclipping
                pass
            s = CigarTools.score(a.cigar)

            if CigarTools.alignment_matches(a.cigar) >= 10 \
                    and CigarTools.match_percent(a.cigar) < self.settings.min_anchor_match:
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

        scores = sorted(scores, key=lambda x: (x[0], x[1]), reverse=True) if len(scores) > 0 else []

        if len(scores) > 1:
            if scores[0][0] != scores[1][0] and scores[0][1] != scores[1][1]:  # not multimap
                clipped = scores[0][2]
                self.split_reads[opposite_breakpoint].add(clipped)
        elif len(scores) == 1:
            clipped = scores[0][2]
            self.split_reads[opposite_breakpoint].add(clipped)

    def call_events(self):
        """
        use the associated evidence and classifications and split the current evidence object
        into more specific objects.
        """
        results = []
        for cls in self.classifications:
            # event type = cls
            new_calls = []
            # try calling by contigs
            # try calling by split reads
            # try calling by flanking reads

    @classmethod
    def _call_by_contigs(cls, ev, classification):
        # resolve the overlap if multi-read alignment
        events = []
        for ctg in ev.contigs:
            for read1, read2 in ctg.alignments:
                if read2 is None:
                    read_events = []
                    for i, t in enumerate(read1.cigar):
                        v, f = t
                        if v not in [CIGAR.I, CIGAR.D, CIGAR.N]:
                            continue
                        if i > 0 and v == read.cigar[i - 1]:
                            start, last, size = read_events[-1]
                            read_events[-1] = (start, i, size + f)
                        else:
                            read_events.append((i, i, f))
                    
                    biggest = max(read_events, key=lambda x: x[2])

                    first_breakpoint = read1.reference_start - 1
                    second_breakpoint = first_breakpoint
                    untemplated_seq = ''
                    
                    seq_pos = 0
                    for i, t in enumerate(read1.cigar):
                        v, f = t
                        if v in [CIGAR.I, CIGAR.S]:
                            seq_pos += f
                            continue
                        elif v == CIGAR.H:
                            continue
                        
                        if i < biggest[0]:
                            first_breakpoint += f
                            second_breakpoint += f
                        elif i >= biggest[0] and i <= biggest[1]:
                            second_breakpoint += f
                            if v == CIGAR.INS:
                                untemplated_seq += read1.query_sequence[seq_pos:seq_pos + f]
                        
                        if v != CIGAR.D:
                            seq_pos += f
                    
                    break1 = Breakpoint(
                        read1.reference_name,
                        first_breakpoint,
                        orient=ORIENT.LEFT,
                        strand=STRAND.NEG if read1.is_reverse else STRAND.POS
                    )
                    break2 = Breakpoint(
                        read1.reference_name,
                        second_breakpoint,
                        orient=ORIENT.LEFT,
                        strand=STRAND.NEG if read1.is_reverse else STRAND.POS
                    )
                    
                    bpp = BreakpointPair(break1, break2, opposing_strands=False, untemplated_sequence=untemplated_seq)
                    
                    new_event = Event(
                        e,
                        bpp,
                        classification,
                        call_method = CALL_METHOD.CONTIG,
                        contig=contig,
                        contig_alignment=(read1, read2)
                    )
                    events.append(new_event)
                    if classification not in [SVTYPE.INS, SVTYPE.DEL, SVTYPE.DUP]:
                        raise NotImplementedError('unsupported classification type', classification)
                else:
                    pass
                    # is there overlap on the query seq?
                    #if classification == SVTYPE.TRANS:
                    #elif classification == SVTYPE.ITRANS:
                    #elif classification == SVTYPE.INV:
                    #elif classification == SVTYPE.INS:
                    #elif classification == SVTYPE.DEL:
                    #elif classification == SVTYPE.DUP:
                    #else:
                    #    raise NotImplementedError('unsupported classification type', classification)
        # filter flanking reads from ev and add to each


    @classmethod
    def _call_by_flanking_reads(cls, ev, classification):
        pass

    @classmethod
    def _call_by_split_reads(cls, ev, classification):
        """
        use split read evidence to resolve bp-level calls for breakpoint pairs (where possible)
        if a bp level call is not possible for one of the breakpoints then returns None
        if no breakpoints can be resolved returns the original event only with NO split read evidence
        also sets the SV type call if multiple are input
        """
        pos1 = {}
        pos2 = {}

        for breakpoint, d in [(self.break1, pos1), (self.break2, pos2)]:
            for read in self.split_reads[breakpoint]:
                try:
                    pos = breakpoint_pos(read, breakpoint.orient)
                    if pos not in d:
                        d[pos] = []
                    d[pos].append(read)
                except AttributeError:
                    pass
            putative_positions = list(d.keys())
            for pos in putative_positions:
                if len(d[pos]) < self.settings.min_splits_reads_resolution:
                    del d[pos]

        linked_pairings = []
        # now pair up the breakpoints with their putative partners
        for first, second in itertools.product(pos1, pos2):
            # can link a pair of breakpoints if the read prefix is in the evidence for both
            # or if the mate id is in the other breakpoint
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(
                self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                bp.break1.pos = Interval(first)
                bp.break2.pos = Interval(second)
                e = Evidence(bp, bam_cache=self.bam_cache, human_reference_genome=self.human_reference_genome)
                e.settings = self.settings
                e.classification = cls
                e.split_reads[e.break1] = set(pos1[first])
                e.split_reads[e.break2] = set(pos2[second])

                if e.linking_evidence() == 0:
                    continue
                for read in self.flanking_reads:
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass

                linked_pairings.append(e)

        psuedo_pairings = set()
        for first in pos1:
            if first not in [temp.break1.start for temp in linked_pairings]:
                for second in pos2:
                    psuedo_pairings.add((first, second))
        for second in pos2:
            if second not in [temp.break2.start for temp in linked_pairings]:
                for first in pos1:
                    psuedo_pairings.add((first, second))

        for first, second in psuedo_pairings:
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(
                self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                bp.break1.pos = Interval(first)
                bp.break2.pos = Interval(second)
                e = Evidence(bp, bam_cache=self.bam_cache, human_reference_genome=self.human_reference_genome)
                e.settings = self.settings
                e.classification = cls
                e.split_reads[e.break1] = set(pos1[first])
                e.split_reads[e.break2] = set(pos2[second])

                for read in self.flanking_reads:
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass

                linked_pairings.append(e)

        if len(linked_pairings) == 0:
            assert(len(pos1.keys()) == 0 or len(pos2.keys() == 0))
            temp = [self.classification] if self.classification is not None else BreakpointPair.classify(
                self.breakpoint_pair)
            for cls in temp:
                bp = self.breakpoint_pair.copy()
                e = Evidence(bp, bam_cache=self.bam_cache, human_reference_genome=self.human_reference_genome)
                e.settings = self.settings
                e.classification = cls

                for read in self.flanking_reads:
                    try:
                        e.add_flanking_read(read)
                    except UserWarning:
                        pass
                linked_pairings.append(e)
        return linked_pairings

    def load_evidence(self, grab_unmapped_partners=True):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
        # TODO transcriptome window gathering

        reads_with_unmapped_partners = set()
        count = 0
        bin_gap_size = self.settings.read_length // 2

        for read in self.bam_cache.fetch(
                '{0}'.format(self.break1.chr),
                self.window1[0],
                self.window1[1],
                read_limit=self.settings.fetch_reads_limit,
                sample_bins=self.settings.fetch_reads_bins,
                bin_gap_size=bin_gap_size,
                cache=True,
                cache_if=partial(
                    lambda filt_sec, r: not filt_sec or (filt_sec and not r.is_secondary),
                    self.settings.filter_secondary_alignments)):
            count += 1
            if read.is_secondary and self.settings.filter_secondary_alignments:
                continue

            if read.is_unmapped or read.mapping_quality < self.settings.min_mapping_quality:
                continue
            try:
                self.add_spanning_read(read)
            except (UserWarning, NotImplementedError):  # TODO: ADD SUPPORT FOR INDELS
                try:
                    self.add_split_read(read, True)
                except UserWarning:
                    pass
            if read.mate_is_unmapped:
                self.half_mapped[self.break1].add(read)
            else:
                try:
                    self.add_flanking_read(read)
                except UserWarning:
                    pass
        count = 0
        for read in self.bam_cache.fetch(
                '{0}'.format(self.break2.chr),
                self.window2[0],
                self.window2[1],
                read_limit=self.settings.fetch_reads_limit,
                sample_bins=self.settings.fetch_reads_bins,
                bin_gap_size=bin_gap_size,
                cache=True,
                cache_if=partial(
                    lambda filt_sec, r: not filt_sec or (filt_sec and not r.is_secondary),
                    self.settings.filter_secondary_alignments)):
            count += 1
            if read.is_secondary and self.settings.filter_secondary_alignments:
                continue

            if read.is_unmapped or read.mapping_quality < self.settings.min_mapping_quality:
                continue
            try:
                self.add_spanning_read(read)
            except (UserWarning, NotImplementedError):  # TODO: ADD SUPPORT FOR INDELS
                try:
                    self.add_split_read(read, False)
                except UserWarning:
                    pass
            if read.mate_is_unmapped:
                self.half_mapped[self.break2].add(read)
            else:
                try:
                    self.add_flanking_read(read)
                except UserWarning:
                    pass


if __name__ == '__main__':
    import doctest
    doctest.testmod()
