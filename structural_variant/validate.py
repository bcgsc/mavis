"""This module is primarily responsible for collecting read evidence from bams
"""

import itertools
from copy import copy as sys_copy
from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.align import CigarTools, nsb_align, breakpoint_pos, assemble
from structural_variant.interval import Interval
from structural_variant.annotate import overlapping_transcripts
from structural_variant.breakpoint import BreakpointPair
import tools.profile_bam as profile_bam
from functools import partial
import math


CALL_METHOD = Vocab(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads', MIXED='split and flanking')


class EventCall:
    """
    class for holding evidence and the related calls since we can't freeze the evidence object
    directly without a lot of copying. Instead we use call objects which are basically
    just a reference to the evidence object and decisions on class, exact breakpoints, etc
    """
    def __init__(self, evidence, classification, breakpoint_pair, call_method, contig=None, alignment=None):
        """
        Args:
            evidence (Evidence): the evidence object we are calling based on
            classification (SVTYPE): the type of structural variant
            breakpoint_pair (BreakpointPair): the breakpoint pair representing the exact breakpoints
            call_method (CALL_METHOD): the way the breakpoints were called
            contig (Contig, optional): the contig used to call the breakpoints (if applicable)
        """
        self.evidence = evidence
        self.classification = classification
        self.breakpoint_pair = breakpoint_pair
        self.contig = contig
        self.call_method = CALL_METHOD.enforce(call_method)
        self.alignment = alignment
        if contig is None and self.call_method == CALL_METHOD.CONTIG:
            raise AttributeError('if the call is by contig then the contig must be provided')

    def count_flanking_support(self):
        """
        counts the flanking read-pair support for the event called

        Returns:
            tuple[int, int, int]:
            * (*int*) - the number of flanking read pairs
            * (*int*) - the median insert size
            * (*int*) - the standard deviation (from the median) of the insert size
        """
        s = self.evidence.settings
        support = set()
        upper_limit = s.median_insert_size + s.stdev_count_abnormal * s.stdev_isize
        lower_limit = s.median_insert_size - s.stdev_count_abnormal * s.stdev_isize

        ihist = {}
        for read in itertools.chain.from_iterable(self.evidence.flanking_reads):
            isize = abs(read.template_length)
            if (self.classification == SVTYPE.INS and isize <= lower_limit) \
                    or (self.classification == SVTYPE.DEL and isize >= upper_limit) \
                    or self.classification not in [SVTYPE.DEL, SVTYPE.INS]:
                support.add(read)
                ihist[isize] = ihist.get(isize, 0) + 1

        stdev = None
        median = None

        if len(support) > 1:
            median = profile_bam.histogram_median(ihist)
            stdev = math.sqrt(profile_bam.histogram_stderr(ihist, median))
        elif len(support) == 1:
            median = list(support)[0]
            stdev = 0

        return len(support), median, stdev

    def count_split_read_support(self):
        """
        counts the split read support for the event called. split reads are only considered to
        be supporting the current call if they exactly match the breakpoint pair associated
        with this call

        Returns:
            tuple[int, int, int]:
            * (*int*) - the number of split reads supporting the first breakpoint
            * (*int*) - the number of split reads supporting the second breakpoint
            * (*int*) - the number of split reads supporting the pairing of these breakpoints
        """
        support1 = set()
        support2 = set()

        for read in self.evidence.split_reads[0]:
            bpos = breakpoint_pos(read, self.breakpoint_pair.break1.orient)
            if Interval.overlaps((bpos, bpos), self.breakpoint_pair.break1):
                support1.add(read)

        for read in self.evidence.split_reads[1]:
            bpos = breakpoint_pos(read, self.breakpoint_pair.break2.orient)
            if Interval.overlaps((bpos, bpos), self.breakpoint_pair.break2):
                support2.add(read)

        links = 0
        query2_names = set([r.query_name for r in support2])
        for read in support1:
            if read.query_name in query2_names:
                links += 1

        return len(support1), len(support2), links


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
            stdev_isize=100,
            stdev_count_abnormal=2,
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
            max_sc_preceeding_anchor=None,
            assembly_min_edge_weight=3,
            assembly_min_remap=3,
    ):
        """
        Args:
            read_length (int, default=125): length of individual reads
            median_insert_size (int, default=380): expected average insert size for paired-end reads
            stdev_isize (int, default=100): the expected deviation in the insert size
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
            fetch_reads_limit (int, default=10000):
                maximum number of reads to loop over for a given event
            filter_secondary_alignments (bool, default=True):
                don't use secondary alignments when reading evidence from the bam file
            sc_extension_stop (int, default=5):
                when extending softclipped, stop given this number of exact consecutive matches
        """
        self.read_length = read_length
        self.median_insert_size = median_insert_size
        self.stdev_isize = stdev_isize
        self.call_error = call_error
        self.min_splits_reads_resolution = min_splits_reads_resolution
        self.min_anchor_exact = min_anchor_exact
        self.min_anchor_fuzzy = min_anchor_fuzzy
        self.min_anchor_size = min(self.min_anchor_exact, self.min_anchor_fuzzy)
        self.min_anchor_match = min_anchor_match
        if self.min_anchor_match > 1 or self.min_anchor_match < 0:
            raise AttributeError('min_anchor_match must be a number between 0 and 1')
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
        self.assembly_min_edge_weight = assembly_min_edge_weight
        self.stdev_count_abnormal = stdev_count_abnormal
        self.assembly_min_remap = assembly_min_remap


class Evidence:
    @property
    def window1(self):
        """(:class:`~structural_variant.interval.Interval`): the window where evidence will be gathered for the first
        breakpoint
        """
        return self.windows[0]

    @property
    def window2(self):
        """(:class:`~structural_variant.interval.Interval`): the window where evidence will be gathered for the second
        breakpoint
        """
        return self.windows[1]

    @property
    def break1(self):
        """(:class:`~structural_variant.breakpoint.Breakpoint`): the first breakpoint"""
        return self.breakpoint_pair.break1

    @property
    def break2(self):
        """(:class:`~structural_variant.breakpoint.Breakpoint`): the second breakpoint"""
        return self.breakpoint_pair.break2

    @property
    def untemplated_sequence(self):
        return self.breakpoint_pair.untemplated_sequence

    @classmethod
    def generate_window(cls, breakpoint, read_length, median_insert_size, call_error, stdev_isize, stdev_count_abnormal):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            read_length (int): the read length
            median_insert_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            stdev_isize:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        fragment_size = read_length * 2 + median_insert_size + stdev_isize * stdev_count_abnormal
        start = breakpoint.start - fragment_size - call_error
        end = breakpoint.end + fragment_size + call_error

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + read_length
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - read_length
        return Interval(start, end)

    @classmethod
    def generate_transcriptome_window(cls, breakpoint, annotations, read_length, median_insert_size, call_error,
                                      stdev_isize, stdev_count_abnormal):
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
            stdev_isize:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        transcripts = overlapping_transcripts(annotations, breakpoint)
        window = cls.generate_window(
            breakpoint, read_length, median_insert_size, call_error, stdev_isize, stdev_count_abnormal)

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
        self.split_reads = (set(), set())
        # flanking reads are read pairs that have a mate within a given breakpoint window and
        # their pair mapped to the opposite breakpoint window
        self.flanking_reads = (set(), set())

        # spanning reads are reads spanning BOTH breakpoints
        self.spanning_reads = set()
        # for each breakpoint stores the number of reads that were read from the associated
        # bamfile for the window surrounding the breakpoint
        self.read_counts = {}
        self.contigs = []

        if self.protocol == PROTOCOL.GENOME:
            w1 = Evidence.generate_window(
                self.break1,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_isize=self.settings.stdev_isize,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            w2 = Evidence.generate_window(
                self.break2,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_isize=self.settings.stdev_isize,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            self.windows = (w1, w2)
        elif self.protocol == PROTOCOL.TRANS:
            w1 = Evidence.generate_transcriptome_window(
                self.break1,
                self.annotations,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_isize=self.settings.stdev_isize,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            w2 = Evidence.generate_transcriptome_window(
                self.break2,
                self.annotations,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_isize=self.settings.stdev_isize,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            self.windows = (w1, w2)
        else:
            raise AttributeError('invalid protocol', self.protocol)
        self.half_mapped = (set(), set())

    def supporting_reads(self):
        """
        convenience method to return all flanking, split and spanning reads associated with an evidence object
        """
        result = set()
        for s in self.flanking_reads:
            result.update(s)
        for s in self.split_reads:
            result.update(s)
        result.update(self.spanning_reads)
        return result

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
        if read.reference_id > read.next_reference_id or \
                (read.reference_id == read.next_reference_id and read.reference_start > read.next_reference_start):
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
        s = self.settings
        if read.is_unmapped or read.mate_is_unmapped:
            raise UserWarning('input read (and its mate) must be mapped')

        # filter by putative event classifications
        classifications = [self.classification] if self.classification else BreakpointPair.classify(
            self.breakpoint_pair)
        classifications = sorted(classifications)

        insert_size = abs(read.template_length)
        upper_limit = s.median_insert_size + s.stdev_isize * s.stdev_count_abnormal
        lower_limit = s.median_insert_size - s.stdev_isize * s.stdev_count_abnormal

        if classifications == sorted([SVTYPE.DEL, SVTYPE.INS]):
            if insert_size <= upper_limit and insert_size >= lower_limit:
                raise UserWarning('insert size is not abnormal. does not support del/ins', insert_size)
        elif classifications == [SVTYPE.DEL]:
            if insert_size <= upper_limit:
                raise UserWarning('insert size is smaller than expected for a deletion type event', insert_size)
        elif classifications == [SVTYPE.INS]:
            if insert_size >= lower_limit:
                raise UserWarning('insert size is larger than expected for an insertion type event', insert_size)

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

        added = False
        if read.reference_start >= w1[0] and read.reference_end <= w1[1] \
                and read.reference_id == self.bam_cache.reference_id(self.break1.chr) \
                and read.next_reference_start >= w2[0] and read.next_reference_start <= w2[1] \
                and read.next_reference_id == self.bam_cache.reference_id(self.break2.chr) \
                and (not self.breakpoint_pair.stranded or not read.is_read1
                    or read.is_reverse == (self.break1.strand == STRAND.NEG)):
            # current read falls in the first breakpoint window, mate in the
            # second
            self.flanking_reads[0].add(read)
            added = True
        if read.reference_start >= w2[0] and read.reference_end <= w2[1] \
                and self.bam_cache.reference_id(self.break2.chr) == read.reference_id \
                and read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                and self.bam_cache.reference_id(self.break1.chr) == read.next_reference_id \
                and (not self.breakpoint_pair.stranded or not read.is_read2
                    or read.is_reverse != (self.break2.strand == STRAND.NEG)):
            # current read falls in the second breakpoint window, mate in the
            # first
            self.flanking_reads[1].add(read)
            added = True
        if not added:
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

        if not read.has_tag('ca') or read.get_tag('ca') != 1:
            read.set_tag('ca', 1, value_type='i')
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
            self.split_reads[0 if first_breakpoint else 1].add(read)

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
            revcomp_sc_align = reverse_complement(read.query_sequence)
            revcomp_sc_align = nsb_align(opposite_breakpoint_ref, revcomp_sc_align)

            for a in revcomp_sc_align:
                a.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
            putative_alignments = revcomp_sc_align

        scores = []

        for a in putative_alignments:
            # a.flag = a.flag ^ 64 ^ 128
            a.flag = a.flag | PYSAM_READ_FLAGS.SECONDARY
            a.set_tag('ca', 1, value_type='i')
            # add information from the original read
            a.reference_start = w[0] - 1 + a.reference_start
            a.reference_id = self.bam_cache.reference_id(opposite_breakpoint.chr)
            # a.query_name = read.query_name + SUFFIX_DELIM + 'clipped-realign'
            a.query_name = read.query_name
            a.set_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN, 1, value_type='i')
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

            if a.reference_id == a.next_reference_id:
                if a.is_read1:
                    a.template_length = (a.next_reference_start + self.settings.read_length) - a.reference_start
                else:
                    a.template_length = a.reference_start - (a.next_reference_start + self.settings.read_length)

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
                self.split_reads[1 if first_breakpoint else 0].add(clipped)
        elif len(scores) == 1:
            clipped = scores[0][2]
            self.split_reads[1 if first_breakpoint else 0].add(clipped)

    def assemble_split_reads(self):
        """
        uses the split reads and the partners of the half mapped reads to create a contig
        representing the sequence across the breakpoints
        """
        strand_specific = self.breakpoint_pair.stranded
        # gather reads for the putative assembly
        assembly_sequences = {}
        for r in itertools.chain.from_iterable(self.split_reads):
            s = r.query_sequence
            if not strand_specific or (r.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                assembly_sequences[s] = assembly_sequences.get(s, set())
                assembly_sequences[s].add(r)
            if not strand_specific or (r.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                temp = reverse_complement(s)
                assembly_sequences[temp] = assembly_sequences.get(temp, set())
                assembly_sequences[temp].add(r)
            # only collect the mates if they are unmapped
            if r.mate_is_unmapped:
                for m in self.bam_cache.get_mate(r):
                    s = m.query_sequence
                    if not strand_specific or (m.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                        assembly_sequences[s] = assembly_sequences.get(s, set())
                        assembly_sequences[s].add(m)
                    if not strand_specific or (m.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                        temp = reverse_complement(s)
                        assembly_sequences[temp] = assembly_sequences.get(temp, set())
                        assembly_sequences[temp].add(m)
        for r in itertools.chain.from_iterable(self.half_mapped):
            try:
                for m in self.bam_cache.get_mate(r):
                    s = m.query_sequence
                    if not strand_specific or (m.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                        assembly_sequences[s] = assembly_sequences.get(s, set())
                        assembly_sequences[s].add(m)
                    if not strand_specific or (m.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.CUSTOM_REALIGN)):
                        temp = reverse_complement(s)
                        assembly_sequences[temp] = assembly_sequences.get(temp, set())
                        assembly_sequences[temp].add(m)
            except KeyError:
                pass
        contigs = assemble(assembly_sequences, min_edge_weight=self.settings.assembly_min_edge_weight)
        filtered_contigs = []
        for c in sorted(contigs, key=lambda x: x.seq):  # sort so that the function is deterministic
            if c.remap_score() < self.settings.assembly_min_remap:
                continue
            filtered_contigs.append(c)
        self.contigs = filtered_contigs

    def call_events(self):
        """
        use the associated evidence and classifications and split the current evidence object
        into more specific objects.
        """
        results = []
        clss = [self.classification] if self.classification else BreakpointPair.classify(self.breakpoint_pair)
        calls = []
        for cls in clss:
            # event type = cls
            new_calls = []
            # try calling by contigs
            calls.extend(Evidence._call_by_contigs(self, cls))
        if len(calls) == 0:
            for cls in clss:
                # try calling by split reads
                calls.extend(Evidence._call_by_supporting_reads(self, cls))
        print('++++++++++++++++++++++++++++++++++++++evidence', self.breakpoint_pair, clss, 'has', len(calls), 'calls')
        for ec in calls:
            print('\n\nevent call', ec.breakpoint_pair, ec.classification, ec.call_method)
            if ec.contig:
                print('contig', ec.contig.seq, ec.contig.remap_score(), ec.contig.score)
            print('flanking counts', ec.count_flanking_support())
            print('split read support', ec.count_split_read_support())
        return calls

    @classmethod
    def _call_by_contigs(cls, ev, classification):
        # resolve the overlap if multi-read alignment
        events = []
        for ctg in ev.contigs:
            for read1, read2 in ctg.alignments:
                bpp = BreakpointPair.call_breakpoint_pair(read1, read2)
                if bpp.opposing_strands != ev.opposing_strands \
                        or (classification == SVTYPE.INS and bpp.untemplated_sequence == '') \
                        or classification not in BreakpointPair.classify(bpp):
                    continue
                new_event = EventCall(
                    ev,
                    classification,
                    bpp,
                    CALL_METHOD.CONTIG,
                    contig=ctg,
                    alignment=(read1, read2)
                )
                events.append(new_event)
        return events

    @classmethod
    def _call_by_flanking_reads(cls, ev, classification, first_breakpoint=None, second_breakpoint=None):
        # for all flanking read pairs mark the farthest possible distance to the breakpoint
        # the start/end of the read on the breakpoint side
        first_positions = set()
        second_positions = set()

        intrachromosomal = (ev.break1.chr == ev.break2.chr)

        for read in ev.flanking_reads[0]:
            s = read.reference_start
            if read.cigar[0][0] == CIGAR.S:
                s -= read.cigar[0][1]
            t = read.reference_end - 1
            if read.cigar[-1][0] == CIGAR.S:
                t += read.cigar[-1][1]
            first_positions.add(s)
            first_positions.add(t)
            second_positions.add(read.next_reference_start)

        for read in ev.flanking_reads[1]:
            s = read.reference_start
            if read.cigar[0][0] == CIGAR.S:
                s -= read.cigar[0][1]
            t = read.reference_end - 1
            if read.cigar[-1][0] == CIGAR.S:
                t += read.cigar[-1][1]
            second_positions.add(s)
            second_positions.add(t)
            first_positions.add(read.next_reference_start)

        max_insert = ev.settings.median_insert_size + ev.settings.stdev_count_abnormal * ev.settings.stdev_isize
        call_method = CALL_METHOD.MIXED

        if first_breakpoint is not None and second_breakpoint is not None:
            raise AttributeError('don\'t need to call by flanking reads')
        elif first_breakpoint is None and second_breakpoint is None:
            call_method = CALL_METHOD.FLANK

        if first_breakpoint is None:
            if len(first_positions) == 0:
                raise UserWarning('no flanking reads available')
            s = min(first_positions)
            t = max(first_positions)
            bp = sys_copy(ev.breakpoint_pair.break1)
            if ev.break1.orient == ORIENT.LEFT:
                bp.start = t
                bp.end = s + max_insert
            elif ev.break1.orient == ORIENT.RIGHT:
                bp.start = s - max_insert
                bp.end = t
            else:
                raise AttributeError('cannot call by flanking if orientation was not given')
            first_breakpoint = bp
            if intrachromosomal:
                if bp.start >= second_breakpoint.start:
                    raise UserWarning('invalid position, interval for first breakpoint must be before the second')
                elif bp.end >= second_breakpoint.start:
                    bp.end = second_breakpoint.start - 1

        if second_breakpoint is None:
            if len(second_positions) == 0:
                raise UserWarning('no flanking reads available')
            s = min(second_positions)
            t = max(second_positions)
            bp = sys_copy(ev.breakpoint_pair.break2)
            if ev.break2.orient == ORIENT.LEFT:
                bp.start = t
                bp.end = s + max_insert
            elif ev.break2.orient == ORIENT.RIGHT:
                bp.start = s - max_insert
                bp.end = t
            else:
                raise AttributeError('cannot call by flanking if orientation was not given')
            second_breakpoint = bp
            if intrachromosomal:
                if bp.end <= first_breakpoint.end:
                    raise UserWarning('invalid position, interval for first breakpoint must be before the second')
                elif bp.start <= first_breakpoint.end:
                    bp.start = first_breakpoint.end + 1

        bpp = ev.breakpoint_pair.copy()
        bpp.break1 = first_breakpoint
        bpp.break2 = second_breakpoint

        return EventCall(ev, classification, bpp, call_method)

    @classmethod
    def _call_by_supporting_reads(cls, ev, classification):
        """
        use split read evidence to resolve bp-level calls for breakpoint pairs (where possible)
        if a bp level call is not possible for one of the breakpoints then returns None
        if no breakpoints can be resolved returns the original event only with NO split read evidence
        also sets the SV type call if multiple are input
        """
        pos1 = {}
        pos2 = {}

        for i, breakpoint, d in [(0, ev.break1, pos1), (1, ev.break2, pos2)]:
            for read in ev.split_reads[i]:
                try:
                    pos = breakpoint_pos(read, breakpoint.orient)
                    if pos not in d:
                        d[pos] = []
                    d[pos].append(read)
                except AttributeError:
                    pass
            putative_positions = list(d.keys())
            for pos in putative_positions:
                if len(d[pos]) < ev.settings.min_splits_reads_resolution:
                    del d[pos]

        linked_pairings = []
        # now pair up the breakpoints with their putative partners
        for first, second in itertools.product(pos1, pos2):
            if ev.break1.chr == ev.break2.chr:
                if first >= second:
                    continue
            links = 0
            read_names = set([r.query_name for r in pos1[first]])
            for read in pos2[second]:
                if read.query_name in read_names:
                    links += 1
            if links == 0:
                continue
            bpp = ev.breakpoint_pair.copy()
            bpp.break1.start = first
            bpp.break1.end = first
            bpp.break2.start = second
            bpp.break2.end = second
            linked_pairings.append(EventCall(ev, classification, bpp, CALL_METHOD.SPLIT))

        f = [p for p in pos1 if p not in [t.breakpoint_pair.break1.start for t in linked_pairings]]
        s = [p for p in pos2 if p not in [t.breakpoint_pair.break2.start for t in linked_pairings]]

        for first, second in itertools.product(f, s):
            if ev.break1.chr == ev.break2.chr:
                if first >= second:
                    continue
            bpp = ev.breakpoint_pair.copy()
            bpp.break1.start = first
            bpp.break1.end = first
            bpp.break2.start = second
            bpp.break2.end = second
            linked_pairings.append(EventCall(ev, classification, bpp, CALL_METHOD.SPLIT))

        if len(linked_pairings) == 0:  # then call by mixed or flanking only
            assert(len(pos1.keys()) == 0 or len(pos2.keys()) == 0)
            fr = len(ev.flanking_reads[0]) + len(ev.flanking_reads[1])
            if fr > 0:
                # if can call the first breakpoint by split
                for pos in pos1:
                    bp = sys_copy(ev.break1)
                    bp.start = pos
                    bp.end = pos
                    try:
                        temp = Evidence._call_by_flanking_reads(ev, classification, first_breakpoint=bp)
                        linked_pairings.append(temp)
                    except UserWarning:
                        pass
                for pos in pos2:
                    bp = sys_copy(ev.break2)
                    bp.start = pos
                    bp.end = pos
                    try:
                        temp = Evidence._call_by_flanking_reads(ev, classification, second_breakpoint=bp)
                        linked_pairings.append(temp)
                    except UserWarning:
                        pass
        return linked_pairings

    def load_evidence(self, grab_unmapped_partners=True):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
        count = 0
        bin_gap_size = self.settings.read_length // 2

        max_dist = max(
            len(Interval.union(self.break1, self.break2)),
            len(self.untemplated_sequence if self.untemplated_sequence else '')
        )

        if max_dist < self.settings.stdev_isize * self.settings.stdev_count_abnormal:
            raise NotImplementedError('evidence gathering for small structural variants is not supported')
            # needs special consideration b/c won't have flanking reads and may have spanning reads

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
                self.half_mapped[0].add(read)
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
                self.half_mapped[1].add(read)
            else:
                try:
                    self.add_flanking_read(read)
                except UserWarning:
                    pass

    @property
    def opposing_strands(self):
        return self.breakpoint_pair.opposing_strands

if __name__ == '__main__':
    import doctest
    doctest.testmod()
