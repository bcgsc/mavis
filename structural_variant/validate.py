"""This module is primarily responsible for collecting read evidence from bams
"""

import itertools
from copy import copy as sys_copy
from .constants import *
from .error import *
from .read_tools import CigarTools, nsb_align, breakpoint_pos
from .assemble import assemble
from .interval import Interval
from .annotate.variant import overlapping_transcripts
from .breakpoint import BreakpointPair, Breakpoint
import statistics
from functools import partial
import math
from argparse import Namespace


DEFAULTS = Namespace(
    stdev_count_abnormal=3,
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
    assembly_min_edge_weight=3,
    assembly_min_remap=3,
    min_linking_split_reads=2,
    min_non_target_aligned_split_reads=1,
    min_flanking_reads_resolution=3
)
"""Namespace: holds the settings for computations with the Evidence objects

.. glossary::

    read_length
        length of reads in the bam file

    stdev_insert_size
        the standard deviation in insert sizes of paired end reads

    median_insert_size
        the median insert size of paired end reads

    stdev_count_abnormal
        the number of standard deviations away from the normal considered expected and therefore not qualifying as
        flanking reads

    call_error
        buffer zone for the evidence window

    min_splits_reads_resolution
        minimum number of split reads required to call a breakpoint by split reads

    min_mapping_quality
        the minimum mapping quality of reads to be used as evidence

    fetch_reads_limit
        maximum number of reads, cap, to loop over for any given evidence window

    fetch_reads_bins
        number of bins to split an evidence window into to ensure more even sampling of high coverage regions

    filter_secondary_alignments
        filter secondary alignments when gathering read evidence

    assembly_min_edge_weight
        when building the initial deBruijn graph edge weights are determined by the frequency of the kmer they represent
        in all the input sequences. The parameter here discards edges to the kmer they represent in all the input
        sequences. The parameter here discards edges to simply the graph if they have a weight less than specified

    assembly_min_remap
        The minimum input sequences that must remap for an assembly to be used

    min_non_target_aligned_split_reads
        The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local alignment
        to the target region to call a breakpoint by split read evidence

    min_linking_split_reads
        The minimum number of split reads which aligned to both breakpoints

    min_flanking_reads_resolution
        the minimum number of flanking reads required to call a breakpoint by flanking evidence
"""


class EventCall(BreakpointPair):
    """
    class for holding evidence and the related calls since we can't freeze the evidence object
    directly without a lot of copying. Instead we use call objects which are basically
    just a reference to the evidence object and decisions on class, exact breakpoints, etc
    """
    def __init__(
        self,
        b1, b2,
        ev,
        classification,
        call_method=None,
        break2_call_method=None,
        contig=None,
        alignment=None,
        stranded=None,
        opposing_strands=None,
        untemplated_sequence=None,
        data=None
    ):
        """
        Args:
            evidence (Evidence): the evidence object we are calling based on
            classification (SVTYPE): the type of structural variant
            breakpoint_pair (BreakpointPair): the breakpoint pair representing the exact breakpoints
            call_method (CALL_METHOD): the way the breakpoints were called
            contig (Contig): the contig used to call the breakpoints (if applicable)
        """
        BreakpointPair.__init__(
            self, b1, b2,
            stranded=stranded if stranded is not None else ev.breakpoint_pair.stranded,
            opposing_strands=opposing_strands if opposing_strands is not None else ev.opposing_strands,
            untemplated_sequence=untemplated_sequence if untemplated_sequence is not None else ev.untemplated_sequence,
            data=data if data is not None else ev.breakpoint_pair.data
        )
        self.evidence = ev
        self.classification = classification
        self.contig = contig
        self.call_method = None
        if contig:
            self.call_method = (CALL_METHOD.CONTIG, CALL_METHOD.CONTIG)
            if call_method or break2_call_method:
                raise AttributeError('contig overrides call method arguments')
        elif break2_call_method:
            self.call_method = (CALL_METHOD.enforce(call_method), CALL_METHOD.enforce(break2_call_method))
        else:
            self.call_method = (CALL_METHOD.enforce(call_method), CALL_METHOD.enforce(call_method))
        self.alignment = alignment

    def count_flanking_support(self):
        """
        counts the flanking read-pair support for the event called

        Returns:
            tuple[int, int, int]:
            * (*int*) - the number of flanking read pairs
            * (*int*) - the median insert size
            * (*int*) - the standard deviation (from the median) of the insert size
        """
        support = set()
        exp_isize_range = self.evidence.expected_insert_size_range()

        insert_sizes = []
        for read in itertools.chain.from_iterable(self.evidence.flanking_reads):
            isize = abs(read.template_length)
            if (self.classification == SVTYPE.INS and isize < exp_isize_range.start) \
                    or (self.classification == SVTYPE.DEL and isize > exp_isize_range.end) \
                    or self.classification not in [SVTYPE.DEL, SVTYPE.INS]:
                support.add(read.query_name)
                insert_sizes.append(isize)

        if len(support) > 0:
            median = statistics.median(insert_sizes)
            err = 0
            for insert in insert_sizes:
                err += math.pow(insert - median, 2)
            err /= len(insert_sizes)
            stdev = math.sqrt(err)
            return len(support), median, stdev
        else:
            return 0, 0, 0

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
        realigns1 = set()
        support2 = set()
        realigns2 = set()

        for read in self.evidence.split_reads[0]:
            try:
                bpos = breakpoint_pos(read, self.break1.orient)
                if Interval.overlaps((bpos, bpos), self.break1):
                    support1.add(read.query_name)
                    if read.has_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN) and \
                            read.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN):
                        realigns1.add(read.query_name)
            except AttributeError:
                pass

        for read in self.evidence.split_reads[1]:
            try:
                bpos = breakpoint_pos(read, self.break2.orient)
                if Interval.overlaps((bpos, bpos), self.break2):
                    support2.add(read.query_name)
                    if read.has_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN) and \
                            read.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN):
                        realigns2.add(read.query_name)
            except AttributeError:
                pass

        return len(support1), len(realigns1), len(support2), len(realigns2), len(support1 & support2)

    def __hash__(self):
        raise NotImplementedError('this object type does not support hashing')

    def __eq__(self, other):
        object.__eq__(self, other)


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
        """(str) the sequence that falls between the two breakpoints and does not map to the reference template"""
        return self.breakpoint_pair.untemplated_sequence

    @classmethod
    def generate_window(
        cls, breakpoint, read_length, median_insert_size, call_error, stdev_insert_size, stdev_count_abnormal
    ):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            read_length (int): the read length
            median_insert_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            stdev_insert_size (int):
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        fragment_size = read_length * 2 + median_insert_size + stdev_insert_size * stdev_count_abnormal
        start = breakpoint.start - fragment_size - call_error
        end = breakpoint.end + fragment_size + call_error

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + read_length
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - read_length
        return Interval(start, end)

    @classmethod
    def generate_transcriptome_window(cls, breakpoint, annotations, read_length, median_insert_size, call_error,
                                      stdev_insert_size, stdev_count_abnormal):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            annotations (dict of str and list of Gene): the set of reference annotations: genes, transcripts, etc
            read_length (int): the read length
            median_insert_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            stdev_insert_size:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        print('generate_transcriptome_window')
        transcripts = overlapping_transcripts(annotations, breakpoint)
        window = cls.generate_window(
            breakpoint, read_length, median_insert_size, call_error, stdev_insert_size, stdev_count_abnormal)

        tgt_left = breakpoint.start - window.start + 1  # amount to expand to the left
        tgt_right = window.end - breakpoint.end + 1  # amount to expand to the right
        print(tgt_left, tgt_right)
        if len(transcripts) == 0:  # case 1. no overlapping transcripts
            print('no overlapping transcripts')
            return window

        intervals = [breakpoint]

        for ust in transcripts:
            print('ust', ust)
            mapping = {}
            cdna_length = sum([len(e) for e in ust.exons])
            s = 1
            for ex in ust.exons:
                mapping[Interval(ex.start, ex.end)] = Interval(s, s + len(ex) - 1)
                s += len(ex)
            reverse_mapping = {}
            for k, v in mapping.items():
                reverse_mapping[v] = k
            print('t', ust.name)
            curr = Interval(breakpoint.start, breakpoint.end)

            if breakpoint.start < ust.start:
                print('before the start')
                curr = curr | Interval(breakpoint.start - tgt_left, breakpoint.start)
            elif breakpoint.start > ust.end:
                print('after the end')
                tgt = tgt_left - (breakpoint.start - ust.end)
                c = Interval.convert_pos(ust.end)
                g = ust.start - (tgt - c)
                if c >= tgt:
                    g = Interval.convert_pos(reverse_mapping, c - tgt + 1)
                curr = curr | Interval(g, breakpoint.start)
            else:
                for ex1, ex2 in zip(ust.exons, ust.exons[1:]):
                    if (breakpoint.start >= ex1.start and breakpoint.start <= ex1.end) or \
                            (breakpoint.start >= ex2.start and breakpoint.start <= ex2.end):
                        print('exonic', ex1, ex2, breakpoint.start, 'left')
                        # in an exon
                        c = Interval.convert_pos(mapping, breakpoint.start)
                        g = ust.start - (tgt_left - c + 1)
                        if c >= tgt_left:
                            g = Interval.convert_pos(reverse_mapping, c - tgt_left + 1)
                        curr = curr | Interval(g, breakpoint.start)
                        break
                    elif breakpoint.start > ex1.end and breakpoint.start < ex2.start:
                        print('intronic', ex1, ex2, breakpoint.start, 'left')
                        isize = breakpoint.start - ex1.end
                        if isize >= tgt_left:
                            curr = curr | Interval(breakpoint.start - tgt_left + 1, breakpoint.start)
                        else:
                            # in an intron
                            tgt = tgt_left - isize
                            print('adjusted tgt', tgt_left, tgt)
                            c = Interval.convert_pos(mapping, ex1.end)
                            print('c', c)
                            g = ust.start - (tgt - c + 1)
                            if c >= tgt:
                                print('dont need the whole thing', c - tgt)
                                g = Interval.convert_pos(reverse_mapping, c - tgt + 1)
                            curr = curr | Interval(g, breakpoint.start)
                            break
            if breakpoint.end > ust.end:
                print('after the end')
                curr = curr | Interval(breakpoint.end, breakpoint.end + tgt_right)
            elif breakpoint.end < ust.start:
                print('before the start')
                tgt = tgt_right - (ust.start - breakpoint.end)
                c = Interval.convert_pos(mapping, ust.end)
                cr = cdna_length - c
                g = ust.end + (tgt - cr)
                if cr > tgt:
                    cr -= tgt
                    g = Interval.convert_pos(reverse_mapping, cdna_length - cr + tgt)
                curr = curr | Interval(g, breakpoint.end)
            else:
                for ex1, ex2 in zip(ust.exons, ust.exons[1:]):
                    if (breakpoint.end >= ex1.start and breakpoint.end <= ex1.end) or \
                            (breakpoint.end >= ex2.start and breakpoint.end <= ex2.end):
                        # in an exon
                        print('exonic', ex1, ex2, breakpoint.end, 'right')
                        c = Interval.convert_pos(mapping, breakpoint.end)
                        rem_c = cdna_length - c
                        g = ust.end + (tgt_right - rem_c)
                        if rem_c > tgt_right:
                            g = Interval.convert_pos(reverse_mapping, c + tgt_right - 1)
                        curr = curr | Interval(breakpoint.end, g)
                        break
                    elif breakpoint.end > ex1.end and breakpoint.end < ex2.start:
                        # in an intron
                        print('intronic', ex1, ex2, breakpoint.end, 'right')
                        isize = ex2.start - breakpoint.end
                        if isize > tgt_right:
                            curr = curr | Interval(breakpoint.end, breakpoint.end + tgt_right - 1)
                        else:
                            tgt = tgt_right - isize
                            c = Interval.convert_pos(mapping, ex2.start)
                            rem_c = cdna_length - c
                            g = ust.end + (tgt - rem_c)
                            if rem_c >= tgt:
                                print('dont need the whole thing')
                                g = Interval.convert_pos(reverse_mapping, c + tgt - 1)
                            curr = curr | Interval(breakpoint.end, g)
                        break
            intervals.append(curr)
        return Interval.union(*intervals)

    def __init__(
            self,
            breakpoint_pair,
            bam_cache,
            REFERENCE_GENOME,
            data={},
            classification=None,
            protocol=PROTOCOL.GENOME,
            annotations={},
            **kwargs):
        """
        Args:
            breakpoint_pair (BreakpointPair): the breakpoint pair to collect evidence for
            bam_cache (BamCache): the bam cache (and assc file) to collect evidence from
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by template/chr name
            data (dict): a dictionary of data to associate with the evidence object
            classification (SVTYPE): the event type
            protocol (PROTOCOL): genome or transcriptome
            **kwargs: named arguments to be passed to EvidenceSettings
        """
        d = dict()
        d.update(DEFAULTS.__dict__)
        d.update(kwargs)
        REQ = ['read_length', 'stdev_insert_size', 'median_insert_size']
        if any([x not in d for x in REQ]):
            raise KeyError('required argument missing', [x for x in REQ if x not in d])
        self.settings = Namespace(**d)
        self.bam_cache = bam_cache
        self.data = data
        self.classification = classification
        self.protocol = PROTOCOL.enforce(protocol)
        self.REFERENCE_GENOME = REFERENCE_GENOME
        self.annotations = annotations

        if self.classification is not None and self.classification not in BreakpointPair.classify(breakpoint_pair):
            raise AttributeError('breakpoint pair improper classification',
                                 BreakpointPair.classify(breakpoint_pair), self.classification)

        if not self.annotations and self.protocol == PROTOCOL.TRANS:
            raise AttributeError('must specify the reference annotations for transcriptome evidence objects')

        self.breakpoint_pair = breakpoint_pair

        if self.break1.orient == ORIENT.NS or self.break2.orient == ORIENT.NS:
            raise NotSpecifiedError(
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
        self.counts = [0, 0]
        self.contigs = []

        if self.protocol == PROTOCOL.GENOME:
            w1 = Evidence.generate_window(
                self.break1,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_insert_size=self.settings.stdev_insert_size,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            w2 = Evidence.generate_window(
                self.break2,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_insert_size=self.settings.stdev_insert_size,
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
                stdev_insert_size=self.settings.stdev_insert_size,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            w2 = Evidence.generate_transcriptome_window(
                self.break2,
                self.annotations,
                read_length=self.settings.read_length,
                median_insert_size=self.settings.median_insert_size,
                call_error=self.settings.call_error,
                stdev_insert_size=self.settings.stdev_insert_size,
                stdev_count_abnormal=self.settings.stdev_count_abnormal
            )
            self.windows = (w1, w2)
        else:
            raise AttributeError('invalid protocol', self.protocol)
        self.half_mapped = (set(), set())

    def expected_insert_size_range(self):
        v = self.settings.stdev_count_abnormal * self.settings.stdev_insert_size
        return Interval(self.settings.median_insert_size - v, self.settings.median_insert_size + v)

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

        .. todo::
            add support for indels
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

        expected_isize = self.expected_insert_size_range()

        if classifications == sorted([SVTYPE.DEL, SVTYPE.INS]):
            if insert_size <= expected_isize.end and insert_size >= expected_isize.start:
                raise UserWarning('insert size is not abnormal. does not support del/ins', insert_size)
        elif classifications == [SVTYPE.DEL]:
            if insert_size <= expected_isize.end:
                raise UserWarning('insert size is smaller than expected for a deletion type event', insert_size)
        elif classifications == [SVTYPE.INS]:
            if insert_size >= expected_isize.start:
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
                and (
                    not self.breakpoint_pair.stranded or not read.is_read1 or
                    read.is_reverse == (self.break1.strand == STRAND.NEG)) \
                and (self.breakpoint_pair.interchromosomal or read.reference_end < read.next_reference_start):
            # current read falls in the first breakpoint window, mate in the
            # second
            self.flanking_reads[0].add(read)
            added = True
        if read.reference_start >= w2[0] and read.reference_end <= w2[1] \
                and self.bam_cache.reference_id(self.break2.chr) == read.reference_id \
                and read.next_reference_start >= w1[0] and read.next_reference_start <= w1[1] \
                and self.bam_cache.reference_id(self.break1.chr) == read.next_reference_id \
                and (
                    not self.breakpoint_pair.stranded or not read.is_read2 or
                    read.is_reverse != (self.break2.strand == STRAND.NEG)) \
                and (self.breakpoint_pair.interchromosomal or read.reference_start > read.next_reference_start):
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
            first_breakpoint (bool): add to the first breakpoint (or second if false)
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
            raise NotSpecifiedError('cannot assign split reads to a breakpoint where the orientation has not been '
                                 'specified')
        if len(primary) < self.settings.min_anchor_exact or len(clipped) < self.settings.min_anchor_exact:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria', primary, clipped)
        elif len(read.query_sequence) - (read.query_alignment_end + 2) < self.settings.min_anchor_exact \
                and (read.query_alignment_start + 1) < self.settings.min_anchor_exact:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria')
        elif len(primary) < self.settings.min_anchor_exact or len(clipped) < self.settings.min_anchor_exact:
            raise UserWarning(
                'split read does not meet the minimum anchor criteria')

        if not read.has_tag('ca') or read.get_tag('ca') != 1:
            read.set_tag('ca', 1, value_type='i')
            # recalculate the read cigar string to ensure M is replaced with = or X
            c = CigarTools.recompute_cigar_mismatch(
                read,
                self.REFERENCE_GENOME[self.bam_cache.chr(read)].seq
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
            raise UserWarning('alignment of too poor quality')
        if CigarTools.longest_exact_match(read.cigar) < self.settings.min_anchor_exact \
                and CigarTools.longest_fuzzy_match(read.cigar, 1) < self.settings.min_anchor_fuzzy:
            raise UserWarning('alignment of too poor quality')
        else:
            self.split_reads[0 if first_breakpoint else 1].add(read)

        # try mapping the soft-clipped portion to the other breakpoint
        w = (opposite_window[0], opposite_window[1])
        opposite_breakpoint_ref = self.REFERENCE_GENOME[
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
            a.set_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN, 1, value_type='i')
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

        if it is not strand specific then sequences are sorted alphanumerically and only the
        first of a pair is kept (paired by sequence)
        """
        strand_specific = self.breakpoint_pair.stranded
        # gather reads for the putative assembly
        assembly_sequences = {}
        for r in itertools.chain.from_iterable(self.split_reads):
            s = r.query_sequence
            if not strand_specific or (r.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                assembly_sequences[s] = assembly_sequences.get(s, set())
                assembly_sequences[s].add(r)
            if not strand_specific or (r.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                temp = reverse_complement(s)
                assembly_sequences[temp] = assembly_sequences.get(temp, set())
                assembly_sequences[temp].add(r)
            # only collect the mates if they are unmapped
            if r.mate_is_unmapped:
                for m in self.bam_cache.get_mate(r):
                    s = m.query_sequence
                    if not strand_specific or (m.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                        assembly_sequences[s] = assembly_sequences.get(s, set())
                        assembly_sequences[s].add(m)
                    if not strand_specific or (m.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                        temp = reverse_complement(s)
                        assembly_sequences[temp] = assembly_sequences.get(temp, set())
                        assembly_sequences[temp].add(m)
        for r in itertools.chain.from_iterable(self.half_mapped):
            try:
                for m in self.bam_cache.get_mate(r):
                    s = m.query_sequence
                    if not strand_specific or (m.is_read1 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                        assembly_sequences[s] = assembly_sequences.get(s, set())
                        assembly_sequences[s].add(m)
                    if not strand_specific or (m.is_read2 and not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN)):
                        temp = reverse_complement(s)
                        assembly_sequences[temp] = assembly_sequences.get(temp, set())
                        assembly_sequences[temp].add(m)
            except KeyError:
                pass
        contigs = assemble(assembly_sequences, min_edge_weight=self.settings.assembly_min_edge_weight)
        filtered_contigs = {}
        for c in sorted(contigs, key=lambda x: x.seq):  # sort so that the function is deterministic
            if c.remap_score() < self.settings.assembly_min_remap:
                continue
            rseq = reverse_complement(c.seq)
            if c.seq not in filtered_contigs and (strand_specific or rseq not in filtered_contigs):
                filtered_contigs[c.seq] = c
        self.contigs = list(filtered_contigs.values())

    def call_events(self):
        """
        use the associated evidence and classifications and split the current evidence object
        into more specific objects.
        """
        results = []
        clss = [self.classification] if self.classification else BreakpointPair.classify(self.breakpoint_pair)
        calls = []
        errors = set()
        for cls in clss:
            # event type = cls
            new_calls = []
            # try calling by contigs
            calls.extend(Evidence._call_by_contigs(self, cls))
        if len(calls) == 0:
            for cls in clss:
                # try calling by split reads
                try:
                    calls.extend(Evidence._call_by_supporting_reads(self, cls))
                except UserWarning as err:
                    errors.add(str(err))
        if len(calls) == 0:
            raise UserWarning(';'.join(sorted(list(errors))))
        return calls

    def clean_flanking_reads(self):
        """ TODO
        cleans the current set of flanking read support
        1. remove any read that appear as evidence for flanking both breakpoints
        2. use the stdev insert size to greedy remove read pairs until the error is within acceptable limits
        """
        for read1, read2 in itertools.product(self.flanking_reads[0], self.flanking_reads[1]):
            if read1 == read2:
                self.flanking_reads[0].remove(read1)
                self.flanking_reads[1].remove(read1)

        # calculate the insert size stdev

    @classmethod
    def _call_by_contigs(cls, ev, classification):
        # resolve the overlap if multi-read alignment
        events = []
        for ctg in ev.contigs:
            for read1, read2 in ctg.alignments:
                try:
                    bpp = BreakpointPair.call_breakpoint_pair(read1, read2)
                except UserWarning as err:
                    continue
                if bpp.opposing_strands != ev.opposing_strands \
                        or (classification == SVTYPE.INS and bpp.untemplated_sequence == '') \
                        or classification not in BreakpointPair.classify(bpp):
                    continue
                new_event = EventCall(
                    bpp.break1,
                    bpp.break2,
                    ev,
                    classification,
                    contig=ctg,
                    alignment=(read1, read2),
                    opposing_strands=bpp.opposing_strands,
                    stranded=bpp.stranded,
                    untemplated_sequence=bpp.untemplated_sequence,
                    data=ev.data
                )
                events.append(new_event)
        return events

    @classmethod
    def _call_by_flanking_reads(cls, ev, classification, first_breakpoint=None, second_breakpoint=None):
        # for all flanking read pairs mark the farthest possible distance to the breakpoint
        # the start/end of the read on the breakpoint side
        first_positions = set()
        second_positions = set()

        expected_isize = ev.expected_insert_size_range()

        for read in ev.flanking_reads[0]:
            if classification == SVTYPE.DEL:
                if abs(read.template_length) <= expected_isize.end:
                    continue
            elif classification == SVTYPE.INS:
                if abs(read.template_length) >= expected_isize.start:
                    continue
            first_positions.add(read.reference_start)
            first_positions.add(read.reference_end)
            second_positions.add(read.next_reference_start)

        for read in ev.flanking_reads[1]:
            if classification == SVTYPE.DEL:
                if abs(read.template_length) <= expected_isize.end:
                    continue
            elif classification == SVTYPE.INS:
                if abs(read.template_length) >= expected_isize.start:
                    continue
            second_positions.add(read.reference_start)
            second_positions.add(read.reference_end - 1)
            first_positions.add(read.next_reference_start)

        break1_call_method = CALL_METHOD.SPLIT if first_breakpoint else CALL_METHOD.FLANK
        break2_call_method = CALL_METHOD.SPLIT if second_breakpoint else CALL_METHOD.FLANK

        if first_breakpoint is not None and second_breakpoint is not None:
            raise AttributeError('don\'t need to call by flanking reads')

        cover1 = None
        cover2 = None
        if len(first_positions) >= ev.settings.min_flanking_reads_resolution:
            cover1 = Interval(min(first_positions), max(first_positions))
        if len(second_positions) >= ev.settings.min_flanking_reads_resolution:
            cover2 = Interval(min(second_positions), max(second_positions))

        if cover1 and cover2:
            if ev.breakpoint_pair.interchromosomal:
                pass
            elif Interval.overlaps(cover1, cover2):
                raise AssertionError(
                    'Cannot resolve {} by flanking reads, flanking read coverage overlaps at the breakpoint: '
                    '{} and {}'.format(classification, cover1, cover2))
            elif cover1.start > cover2.start:
                raise AssertionError(
                    'Cannot resolve {} by flanking reads. Region of coverage for breakpoint1 falls ahead of the region '
                    'of coverage for breakpoint2'.format(classification))
            # have coverage for flanking evidence for both breakpoints
        else:
            raise UserWarning(
                'Unable to call {} by flanking reads, insufficient flanking reads available'.format(classification))

        if first_breakpoint is None:
            shift = max([0, len(cover1) - ev.settings.read_length])
            if shift > expected_isize.end:
                raise AssertionError(
                    'Cannot resolve {} by flanking reads. Flanking coverage region is larger than expected. It is'
                    'likely that these flanking reads are not all supporting the same event'.format(classification))

            if ev.break1.orient == ORIENT.LEFT:
                """
                                                                  *  *
                coverage        ----------------=======---------------======----------------
                other break     --------------------------------|=|-------------------------
                raw window      ----------------------|==========================|----------
                refined window  ----------------------|===========|-------------------------
                """
                right_side_bound = cover1.end + expected_isize.end - shift
                if not ev.breakpoint_pair.interchromosomal:
                    if cover2:
                        right_side_bound = min([right_side_bound, cover2.start - 1])
                    if second_breakpoint:
                        if second_breakpoint.end < cover1.end:
                            raise AssertionError(
                                'Cannot call by flanking reads. Coverage region for the first breakpoint extends past '
                                'the call region for the second breakpoint')
                        right_side_bound = min([right_side_bound, second_breakpoint.end])

                first_breakpoint = Breakpoint(
                    ev.break1.chr,
                    cover1.end,
                    right_side_bound,
                    orient=ev.break1.orient,
                    strand=ev.break1.strand
                )
            elif ev.break1.orient == ORIENT.RIGHT:
                first_breakpoint = Breakpoint(
                    ev.break1.chr,
                    max([cover1.start - expected_isize.end + shift, 0]),
                    max([cover1.start, 0]),
                    orient=ev.break1.orient,
                    strand=ev.break1.strand
                )
            else:
                raise NotSpecifiedError('Cannot call by flanking if orientation was not given')

        if second_breakpoint is None:
            shift = max([0, len(cover2) - ev.settings.read_length])
            if shift > expected_isize.end:
                raise AssertionError(
                    'Cannot resolve {} by flanking reads. Flanking coverage region is larger than expected. It is '
                    'likely that these flanking reads are not all supporting the same event'.format(classification))
            if ev.break2.orient == ORIENT.LEFT:
                second_breakpoint = Breakpoint(
                    ev.break2.chr,
                    cover2.end,
                    cover2.end + expected_isize.end - shift,
                    orient=ev.break2.orient,
                    strand=ev.break2.strand
                )
            elif ev.break2.orient == ORIENT.RIGHT:
                """
                                                     * *
                coverage        ----------------=======---------------======----------------
                other break     ---------------------|====|---------------------------------
                raw window      ------------------|===================|---------------------
                refined window  -----------------------|==============|---------------------
                """
                left_side_bound = cover2.start - expected_isize.end + shift
                if not ev.breakpoint_pair.interchromosomal:
                    if cover1:
                        left_side_bound = max([left_side_bound, cover1.end + 1])
                    if first_breakpoint:
                        if first_breakpoint.start > cover2.start:
                            raise AssertionError(
                                'Cannot call by flanking reads. Coverage region for the first breakpoint extends past '
                                'the call region for the second breakpoint')
                        left_side_bound = max([left_side_bound, first_breakpoint.start])

                second_breakpoint = Breakpoint(
                    ev.break2.chr,
                    left_side_bound,
                    cover2.start,
                    orient=ev.break2.orient,
                    strand=ev.break2.strand
                )
            else:
                raise NotSpecifiedError('Cannot call {} by flanking if orientation was not given'.format(classification))
        return first_breakpoint, second_breakpoint

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
                else:
                    count = 0
                    for r in d[pos]:
                        if not r.has_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN) or \
                                not r.get_tag(PYSAM_READ_FLAGS.FORCED_TARGET_ALIGN):
                            count += 1
                    if count < ev.settings.min_non_target_aligned_split_reads:
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
            if links < ev.settings.min_linking_split_reads:
                continue
            first_breakpoint = Breakpoint(ev.break1.chr, first, strand=ev.break1.strand, orient=ev.break1.orient)
            second_breakpoint = Breakpoint(ev.break2.chr, second, strand=ev.break2.strand, orient=ev.break2.orient)
            call = EventCall(
                first_breakpoint, second_breakpoint, ev, classification,
                call_method=CALL_METHOD.SPLIT,
                data=ev.data)
            linked_pairings.append(call)

        f = [p for p in pos1 if p not in [t.break1.start for t in linked_pairings]]
        s = [p for p in pos2 if p not in [t.break2.start for t in linked_pairings]]

        for first, second in itertools.product(f, s):
            if ev.break1.chr == ev.break2.chr:
                if first >= second:
                    continue
            first_breakpoint = Breakpoint(ev.break1.chr, first, strand=ev.break1.strand, orient=ev.break1.orient)
            second_breakpoint = Breakpoint(ev.break2.chr, second, strand=ev.break2.strand, orient=ev.break2.orient)
            call = EventCall(
                first_breakpoint, second_breakpoint, ev, classification,
                call_method=CALL_METHOD.SPLIT,
                data=ev.data
            )
            linked_pairings.append(call)

        if len(linked_pairings) == 0:  # then call by mixed or flanking only
            fr = len(ev.flanking_reads[0]) + len(ev.flanking_reads[1])
            error_messages = set()
            if fr > 0:
                # if can call the first breakpoint by split
                for pos in pos1:
                    bp = sys_copy(ev.break1)
                    bp.start = pos
                    bp.end = pos
                    try:
                        f, s = Evidence._call_by_flanking_reads(ev, classification, first_breakpoint=bp)
                        call = EventCall(
                            f, s, ev, classification,
                            call_method=CALL_METHOD.SPLIT,
                            break2_call_method=CALL_METHOD.FLANK,
                            data=ev.data
                        )
                        linked_pairings.append(call)
                    except AssertionError as err:
                        error_messages.add(str(err))
                    except UserWarning as err:
                        error_messages.add(str(err))

                for pos in pos2:
                    bp = sys_copy(ev.break2)
                    bp.start = pos
                    bp.end = pos
                    try:
                        f, s = Evidence._call_by_flanking_reads(ev, classification, second_breakpoint=bp)
                        call = EventCall(
                            f, s, ev, classification,
                            call_method=CALL_METHOD.FLANK,
                            break2_call_method=CALL_METHOD.SPLIT,
                            data=ev.data
                        )
                        linked_pairings.append(call)
                    except AssertionError as err:
                        error_messages.add(str(err))
                    except UserWarning as err:
                        error_messages.add(str(err))

                if len(linked_pairings) == 0:  # call by flanking only
                    try:
                        f, s = Evidence._call_by_flanking_reads(ev, classification)
                        call = EventCall(
                            f, s, ev, classification,
                            call_method=CALL_METHOD.FLANK,
                            data=ev.data
                        )
                        linked_pairings.append(call)
                    except AssertionError as err:
                        error_messages.add(str(err))
                    except UserWarning as err:
                        error_messages.add(str(err))
        if len(linked_pairings) == 0:
            raise UserWarning(';'.join(list(error_messages)))
        return linked_pairings

    def load_evidence(self, grab_unmapped_partners=True):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering

        .. todo::
            support gathering evidence for small structural variants
        """
        bin_gap_size = self.settings.read_length // 2

        max_dist = max(
            len(Interval.union(self.break1, self.break2)),
            len(self.untemplated_sequence if self.untemplated_sequence else '')
        )

        if max_dist < self.settings.stdev_insert_size * self.settings.stdev_count_abnormal:
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
            self.counts[0] += 1
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
            self.counts[1] += 1
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
