import itertools
import logging
from .constants import DEFAULTS
from ..assemble import assemble
from ..bam import cigar as _cigar
from ..bam import read as _read
from ..bam.cache import BamCache
from ..breakpoint import BreakpointPair
from ..constants import CIGAR, COLUMNS, NA_MAPPING_QUALITY, ORIENT, PROTOCOL, PYSAM_READ_FLAGS, reverse_complement, STRAND, SVTYPE
from ..error import NotSpecifiedError
from ..interval import Interval
from ..util import DEVNULL


class Evidence(BreakpointPair):

    @property
    def min_expected_fragment_size(self):
        # cannot be negative
        return int(round(max([self.median_fragment_size - self.stdev_fragment_size * self.stdev_count_abnormal, 0]), 0))

    @property
    def max_expected_fragment_size(self):
        return int(round(self.median_fragment_size + self.stdev_fragment_size * self.stdev_count_abnormal, 0))

    def __init__(
            self,
            break1, break2,
            bam_cache,
            reference_genome,
            read_length,
            stdev_fragment_size,
            median_fragment_size,
            stranded=False,
            opposing_strands=None,
            untemplated_seq=None,
            data={},
            classification=None,
            **kwargs):
        """
        Args:
            breakpoint_pair (BreakpointPair): the breakpoint pair to collect evidence for
            bam_cache (BamCache): the bam cache (and assc file) to collect evidence from
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`):
              dict of reference sequence by template/chr name
            data (dict): a dictionary of data to associate with the evidence object
            classification (SVTYPE): the event type
            protocol (PROTOCOL): genome or transcriptome
        """
        # initialize the breakpoint pair
        self.bam_cache = bam_cache
        self.stranded = stranded and bam_cache.stranded
        BreakpointPair.__init__(
            self, break1, break2,
            stranded=stranded,
            opposing_strands=opposing_strands,
            untemplated_seq=untemplated_seq,
            **data
        )
        # check that the breakpoints are within the reference length
        if reference_genome:
            if self.break1.start < 1 or self.break1.end > len(reference_genome[self.break1.chr].seq):
                raise ValueError('Breakpoint {}-{} is outside the range of the reference sequence {} (1-{})'.format(
                    self.break1.start, self.break1.end, self.break1.chr, len(reference_genome[self.break1.chr].seq)))
            if self.break2.start < 1 or self.break2.end > len(reference_genome[self.break2.chr].seq):
                raise ValueError('Breakpoint {}-{} is outside the range of the reference sequence {} (1-{})'.format(
                    self.break2.start, self.break2.end, self.break2.chr, len(reference_genome[self.break2.chr].seq)))
        defaults = dict()
        for arg in kwargs:
            if arg not in DEFAULTS:
                raise AttributeError('unrecognized attribute', arg)
        defaults.update(DEFAULTS.items())
        kwargs.setdefault('assembly_max_kmer_size', int(read_length * 0.7))
        defaults.update(kwargs)  # input arguments should override the defaults
        for arg, val in defaults.items():
            setattr(self, arg, val)

        self.bam_cache = bam_cache
        self.classification = classification
        self.reference_genome = reference_genome
        self.read_length = read_length
        self.stdev_fragment_size = stdev_fragment_size
        self.median_fragment_size = median_fragment_size
        self.compatible_window1 = None
        self.compatible_window2 = None

        if self.classification is not None and self.classification not in BreakpointPair.classify(self):
            raise AttributeError(
                'breakpoint pair improper classification', BreakpointPair.classify(self), self.classification)

        if self.break1.orient == ORIENT.NS or self.break2.orient == ORIENT.NS:
            raise NotSpecifiedError(
                'input breakpoint pair must specify strand and orientation. Cannot be \'not specified'
                '\' for evidence gathering')

        self.split_reads = (set(), set())
        self.flanking_pairs = set()
        self.compatible_flanking_pairs = set()
        self.spanning_reads = set()
        # for each breakpoint stores the number of reads that were read from the associated
        # bamfile for the window surrounding the breakpoint
        self.counts = [0, 0]  # has to be a list to assign
        self.contigs = []

        self.half_mapped = (set(), set())

        try:
            self.compute_fragment_size(None, None)
        except NotImplementedError:
            raise NotImplementedError('abstract class cannot be initialized')
        except BaseException:
            pass

    @staticmethod
    def distance(start, end):
        return Interval(abs(end - start))

    @staticmethod
    def traverse(start, distance, direction):
        if direction == ORIENT.LEFT:
            return Interval(start - distance)
        return Interval(start + distance)

    def collect_from_outer_window(self):
        """
        determines if evidence should be collected from the outer window (looking for flanking evidence)
        or should be limited to the inner window (split/spanning/contig only)

        Returns:
            bool: True or False
        """
        if self.interchromosomal:
            return True
        elif len(self.break1 | self.break2) >= self.outer_window_min_event_size:
            return True
        return False

    def standardize_read(self, read):
        # recomputing to standardize b/c split reads can be used to call breakpoints exactly
        read.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
        # recalculate the read cigar string to ensure M is replaced with = or X
        cigar = _cigar.recompute_cigar_mismatch(
            read,
            self.reference_genome[self.bam_cache.get_read_reference_name(read)].seq
        )
        prefix = 0
        try:
            cigar, prefix = _cigar.extend_softclipping(cigar, self.min_anchor_exact)
        except AttributeError:
            pass
        read.cigar = _cigar.join(cigar)
        read.cigar = _cigar.merge_internal_events(
            read.cigar,
            inner_anchor=self.contig_aln_merge_inner_anchor,
            outer_anchor=self.contig_aln_merge_outer_anchor
        )
        read.reference_start = read.reference_start + prefix

        # makes sure all indels are called as far 'right' as possible
        read.cigar = _cigar.hgvs_standardize_cigar(
            read, self.reference_genome[self.bam_cache.get_read_reference_name(read)].seq)
        return read

    def putative_event_types(self):
        """
        Returns:
            list of :class:`~mavis.constants.SVTYPE`: list of the possible classifications
        """
        if self.classification:
            return {self.classification}
        return BreakpointPair.classify(self)

    @property
    def compatible_type(self):
        if SVTYPE.INS in self.putative_event_types():
            return SVTYPE.DUP
        elif SVTYPE.DUP in self.putative_event_types():
            return SVTYPE.INS
        return None

    def compute_fragment_size(self, read, mate):
        """
        Args:
            read (pysam.AlignedSegment):
            mate (pysam.AlignedSegment):
        Returns:
            Interval: interval representing the range of possible fragment sizes for this read pair
        """
        raise NotImplementedError('abstract method must be overridden')

    def supporting_reads(self):
        """
        convenience method to return all flanking, split and spanning reads associated with an evidence object
        """
        result = set()
        for read_set in self.flanking_pairs:
            result.update(read_set)
        for read_set in self.split_reads:
            result.update(read_set)
        result.update(self.spanning_reads)
        return result

    def collect_spanning_read(self, read):
        """
        spanning read: a read covering BOTH breakpoints

        This is only applicable to small events. Do not need to look for soft clipped reads
        here since they will be collected already

        Args:
            read (pysam.AlignedSegment): the putative spanning read

        Returns:
            bool:
                - True: the read was collected and stored in the current evidence object
                - False: the read was not collected
        """
        if self.interchromosomal:
            return False
        elif not Interval.overlaps(self.inner_window1, self.inner_window2):
            # too far apart to call spanning reads
            return False

        if self.stranded:
            strand = _read.sequenced_strand(read, self.strand_determining_read)
            if strand != self.break1.strand and strand != self.break2.strand:
                return False

        combined = self.inner_window1 & self.inner_window2
        read_interval = Interval(read.reference_start + 1, read.reference_end)

        if Interval.overlaps(combined, read_interval):

            if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):

                read = self.standardize_read(read)
            # in the correct position, now determine if it can support the event types
            for event_type in self.putative_event_types():
                if event_type in [SVTYPE.DUP, SVTYPE.INS]:
                    if CIGAR.I in [c[0] for c in read.cigar]:
                        self.spanning_reads.add(read)
                        return True
                elif event_type == SVTYPE.DEL:
                    if CIGAR.D in [c[0] for c in read.cigar]:
                        self.spanning_reads.add(read)
                        return True
                elif event_type == SVTYPE.INV:
                    if CIGAR.X in [c[0] for c in read.cigar]:
                        return True
        return False

    def collect_compatible_flanking_pair(self, read, mate, compatible_type):
        """
        checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
        this event

        Args:
            read (pysam.AlignedSegment): the read to add
            mate (pysam.AlignedSegment): the mate
            compatible_type (SVTYPE): the type we are collecting for

        Returns:
            bool:
                - True: the pair was collected and stored in the current evidence object
                - False: the pair was not collected
        Raises:
            ValueError: if the input reads are not a valid pair

        see :ref:`theory - types of flanking evidence <theory-compatible-flanking-pairs>`
        """
        if read.is_unmapped or mate.is_unmapped or read.query_name != mate.query_name or read.is_read1 == mate.is_read1:
            raise ValueError('input reads must be a mapped and mated pair')
        if not self.compatible_window1:
            raise ValueError('compatible windows were not given')
        if self.interchromosomal:
            raise NotImplementedError('interchromosomal events do not have compatible flanking pairs')
        # check that the references are right for the pair
        if read.reference_id != read.next_reference_id:
            return False
        elif read.mapping_quality < self.min_mapping_quality or mate.mapping_quality < self.min_mapping_quality:
            return False
        if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            read = self.standardize_read(read)
        if not mate.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not mate.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            mate = self.standardize_read(mate)
        # order the read pairs so that they are in the same order that we expect for the breakpoints
        if read.reference_start > mate.reference_start:
            read, mate = mate, read

        if self.bam_cache.get_read_reference_name(read) != self.break1.chr or \
                self.bam_cache.get_read_reference_name(mate) != self.break2.chr:
            return False

        # check if this read falls in the first breakpoint window
        iread = Interval(read.reference_start + 1, read.reference_end)
        imate = Interval(mate.reference_start + 1, mate.reference_end)

        if self.stranded:
            strand1 = _read.sequenced_strand(read, self.strand_determining_read)
            strand2 = _read.sequenced_strand(mate, self.strand_determining_read)
            if strand1 != self.break1.strand or strand2 != self.break2.strand:
                return False

        # check that the pair orientation is correct
        if not _read.orientation_supports_type(read, compatible_type):
            return False

        # check that the fragment size is reasonable
        fragment_size = self.compute_fragment_size(read, mate)

        if compatible_type == SVTYPE.DEL:
            if fragment_size.end <= self.max_expected_fragment_size:
                return False
        elif compatible_type == SVTYPE.INS:
            if fragment_size.start >= self.min_expected_fragment_size:
                return False

        # check that the positions of the reads and the strands make sense
        if Interval.overlaps(iread, self.compatible_window1) and \
                Interval.overlaps(imate, self.compatible_window2):
            self.compatible_flanking_pairs.add((read, mate))
            return True

        return False

    def collect_flanking_pair(self, read, mate):
        """
        checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
        this event

        Args:
            read (pysam.AlignedSegment): the read to add
            mate (pysam.AlignedSegment): the mate

        Returns:
            bool:
                - True: the pair was collected and stored in the current evidence object
                - False: the pair was not collected
        Raises:
            ValueError: if the input reads are not a valid pair

        see :ref:`theory - types of flanking evidence <theory-types-of-flanking-evidence>`
        """
        if read.is_unmapped or mate.is_unmapped:
            raise ValueError('input reads must be a mapped and mated pair. One or both of the reads is unmapped')
        elif read.query_name != mate.query_name:
            raise ValueError('input reads must be a mapped and mated pair. The query names do not match')
        elif abs(read.template_length) != abs(mate.template_length):
            raise ValueError(
                'input reads must be a mapped and mated pair. The template lengths (abs value) do not match',
                abs(read.template_length), abs(mate.template_length))
        elif read.mapping_quality < self.min_mapping_quality or mate.mapping_quality < self.min_mapping_quality:
            return False  # do not meet the minimum mapping quality requirements
        # check that the references are right for the pair
        if self.interchromosomal:
            if read.reference_id == read.next_reference_id:
                return False
        elif read.reference_id != read.next_reference_id:
            return False

        if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            read = self.standardize_read(read)
        if not mate.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not mate.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            mate = self.standardize_read(mate)
        # order the read pairs so that they are in the same order that we expect for the breakpoints
        if read.reference_id != mate.reference_id:
            if self.bam_cache.get_read_reference_name(read) > self.bam_cache.get_read_reference_name(mate):
                read, mate = mate, read
        elif read.reference_start > mate.reference_start:
            read, mate = mate, read

        if self.bam_cache.get_read_reference_name(read) != self.break1.chr or \
                self.bam_cache.get_read_reference_name(mate) != self.break2.chr:
            return False
        if any([
            self.break1.orient == ORIENT.LEFT and read.is_reverse,
            self.break1.orient == ORIENT.RIGHT and not read.is_reverse,
            self.break2.orient == ORIENT.LEFT and mate.is_reverse,
            self.break2.orient == ORIENT.RIGHT and not mate.is_reverse
        ]):
            return False
        # check if this read falls in the first breakpoint window
        iread = Interval(read.reference_start + 1, read.reference_end)
        imate = Interval(mate.reference_start + 1, mate.reference_end)

        if self.stranded:
            strand1 = _read.sequenced_strand(read, self.strand_determining_read)
            strand2 = _read.sequenced_strand(mate, self.strand_determining_read)

            if strand1 != self.break1.strand or strand2 != self.break2.strand:
                return False

        for event_type in self.putative_event_types():

            # check that the pair orientation is correct
            if not _read.orientation_supports_type(read, event_type):
                continue

            # check that the fragment size is reasonable
            fragment_size = self.compute_fragment_size(read, mate)

            if event_type == SVTYPE.DEL:
                if fragment_size.end <= self.max_expected_fragment_size:
                    continue
            elif event_type == SVTYPE.INS:
                if fragment_size.start >= self.min_expected_fragment_size:
                    continue

            # check that the positions of the reads and the strands make sense
            if Interval.overlaps(iread, self.outer_window1) and Interval.overlaps(imate, self.outer_window2):
                self.flanking_pairs.add((read, mate))
                return True

        return False

    def collect_half_mapped(self, read, mate):
        """
        Args:
            read (pysam.AlignedSegment): the read to add
            mate (pysam.AlignedSegment): the unmapped mate

        Returns:
            bool:
                - True: the read was collected and stored in the current evidence object
                - False: the read was not collected
        Raises:
            AssertionError: if the mate is not unmapped
        """
        if not mate.is_unmapped:
            raise AssertionError('expected the mate to be unmapped')
        read_itvl = Interval(read.reference_start + 1, read.reference_end)
        added = False
        if self.break1.chr == read.reference_name and Interval.overlaps(self.outer_window1, read_itvl):
            self.half_mapped[0].add(mate)
            added = True
        if self.break2.chr == read.reference_name and Interval.overlaps(self.outer_window2, read_itvl):
            self.half_mapped[1].add(mate)
            added = True
        return added

    def collect_split_read(self, read, first_breakpoint):
        """
        adds a split read if it passes the criteria filters and raises a warning if it does not

        Args:
            read (pysam.AlignedSegment): the read to add
            first_breakpoint (bool): add to the first breakpoint (or second if false)
        Returns:
            bool:
                - True: the read was collected and stored in the current evidence object
                - False: the read was not collected
        Raises:
            NotSpecifiedError: if the breakpoint orientation is not specified
        """
        breakpoint = self.break1 if first_breakpoint else self.break2
        window = self.inner_window1 if first_breakpoint else self.inner_window2
        opposite_breakpoint = self.break2 if first_breakpoint else self.break1
        opposite_window = self.inner_window2 if first_breakpoint else self.inner_window1

        if read.cigar[0][0] != CIGAR.S and read.cigar[-1][0] != CIGAR.S:
            return False  # split read is not softclipped
        elif breakpoint.orient == ORIENT.LEFT and read.cigar[-1][0] != CIGAR.S:
            return False  # split read is not softclipped
        elif breakpoint.orient == ORIENT.RIGHT and read.cigar[0][0] != CIGAR.S:
            return False  # split read is not softclipped

        # the first breakpoint of a BreakpointPair is always the lower breakpoint
        # if this is being added to the second breakpoint then we'll need to check if the
        # read soft-clipping needs to be adjusted

        # check if the read falls within the evidence collection window
        # recall that pysam is 0-indexed and the window is 1-indexed
        if read.reference_start > window.end - 1 or read.reference_end < window.start \
                or self.bam_cache.get_read_reference_name(read) != breakpoint.chr:
            return False  # read not in breakpoint evidence window
        # can only enforce strand if both the breakpoint and the bam are stranded
        if self.stranded and self.bam_cache.stranded:
            strand = _read.sequenced_strand(read, strand_determining_read=self.strand_determining_read)
            if strand != breakpoint.strand:
                return False  # split read not on the appropriate strand
        unused = ''
        primary = ''
        clipped = ''
        if breakpoint.orient == ORIENT.LEFT:
            unused = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            # end is exclusive in pysam
            clipped = read.query_sequence[read.query_alignment_end:]
        elif breakpoint.orient == ORIENT.RIGHT:
            clipped = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            unused = read.query_sequence[read.query_alignment_end:]
        else:
            raise NotSpecifiedError(
                'cannot assign split reads to a breakpoint where the orientation has not been specified')
        if len(primary) + len(clipped) + len(unused) != len(read.query_sequence):
            raise AssertionError(
                'unused, primary, and clipped sequences should make up the original sequence',
                unused, primary, clipped, read.query_sequence, len(read.query_sequence))

        if len(primary) < self.min_anchor_exact or len(clipped) < self.min_softclipping:
            # split read does not meet the minimum anchor criteria
            return False
        if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            read = self.standardize_read(read)
        # data quality filters
        if _cigar.alignment_matches(read.cigar) >= self.min_sample_size_to_apply_percentage \
                and _cigar.match_percent(read.cigar) < self.min_anchor_match:
            return False  # too poor quality of an alignment
        if _cigar.longest_exact_match(read.cigar) < self.min_anchor_exact \
                and _cigar.longest_fuzzy_match(read.cigar, self.fuzzy_mismatch_number) < self.min_anchor_fuzzy:
            return False  # too poor quality of an alignment
        else:
            self.split_reads[0 if first_breakpoint else 1].add(read)

        # try mapping the soft-clipped portion to the other breakpoint
        w = (opposite_window[0], opposite_window[1])
        opposite_breakpoint_ref = self.reference_genome[opposite_breakpoint.chr].seq[w[0] - 1: w[1]]

        putative_alignments = None
        # figure out how much of the read must match when remaped
        min_match_tgt = read.cigar[-1][1] if breakpoint.orient == ORIENT.LEFT else read.cigar[0][1]
        min_match_tgt = min(min_match_tgt * self.min_anchor_match, min_match_tgt - 1) / len(read.query_sequence)
        if not self.opposing_strands:  # same strand
            sc_align = _read.nsb_align(
                opposite_breakpoint_ref, read.query_sequence, min_consecutive_match=self.min_anchor_exact,
                min_match=min_match_tgt, min_overlap_percent=min_match_tgt)  # split half to this side

            for alignment in sc_align:
                alignment.flag = read.flag
            putative_alignments = sc_align
        else:
            # should align opposite the current read
            revcomp_sc_align = reverse_complement(read.query_sequence)
            revcomp_sc_align = _read.nsb_align(
                opposite_breakpoint_ref, revcomp_sc_align, min_consecutive_match=self.min_anchor_exact,
                min_match=min_match_tgt, min_overlap_percent=min_match_tgt)

            for alignment in revcomp_sc_align:
                alignment.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
            putative_alignments = revcomp_sc_align

        scores = []
        for alignment in putative_alignments:  # loop over the alignments
            alignment.flag = alignment.flag | PYSAM_READ_FLAGS.SUPPLEMENTARY
            # set this flag so we don't recompute the cigar multiple
            alignment.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
            # add information from the original read
            alignment.reference_start = w[0] - 1 + alignment.reference_start
            alignment._reference_name = opposite_breakpoint.chr  # must be set since not associated with an alignment file
            alignment.reference_id = self.bam_cache.reference_id(opposite_breakpoint.chr)
            alignment.query_name = read.query_name
            alignment.set_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1, value_type='i')
            alignment.next_reference_start = read.next_reference_start
            alignment.next_reference_id = read.next_reference_id
            alignment.mapping_quality = NA_MAPPING_QUALITY
            try:
                cigar, offset = _cigar.extend_softclipping(alignment.cigar, self.min_anchor_exact)
                alignment.cigar = cigar
                alignment.reference_start = alignment.reference_start + offset
            except AttributeError:
                # if the matches section is too small you can't extend the
                # softclipping
                pass
            s = _cigar.score(alignment.cigar)
            alignment.cigar = _cigar.join(alignment.cigar)
            if alignment.reference_id == alignment.next_reference_id:
                # https://samtools.github.io/hts-specs/SAMv1.pdf
                # unsigned observed template length equals the number of bases from the leftmost
                # mapped base to the rightmost mapped base
                tlen = abs(alignment.reference_start - alignment.next_reference_start) + 1
                if alignment.reference_start < alignment.next_reference_start:
                    alignment.template_length = tlen
                else:
                    alignment.template_length = -1 * tlen
            else:
                alignment.template_length = 0
            if _cigar.alignment_matches(alignment.cigar) >= self.min_sample_size_to_apply_percentage \
                    and _cigar.match_percent(alignment.cigar) < self.min_anchor_match:
                continue
            if _cigar.longest_exact_match(alignment.cigar) < self.min_anchor_exact \
                    and _cigar.longest_fuzzy_match(alignment.cigar, self.fuzzy_mismatch_number) < self.min_anchor_fuzzy:
                continue
            if self.max_sc_preceeding_anchor is not None:
                if opposite_breakpoint.orient == ORIENT.LEFT:
                    if alignment.cigar[0][0] == CIGAR.S and alignment.cigar[0][1] > self.max_sc_preceeding_anchor:
                        continue
                elif opposite_breakpoint.orient == ORIENT.RIGHT:
                    if alignment.cigar[-1][0] == CIGAR.S and alignment.cigar[-1][1] > self.max_sc_preceeding_anchor:
                        continue
            alignment.set_key()  # set the hash key before we add the read as evidence
            scores.append((s, _cigar.match_percent(alignment.cigar), alignment))

        scores = sorted(scores, key=lambda x: (x[0], x[1]), reverse=True) if scores else []

        if len(scores) > 1:
            if scores[0][0] != scores[1][0] and scores[0][1] != scores[1][1]:
                # not multimap, pick highest scoring alignment
                clipped = scores[0][2]
                self.split_reads[1 if first_breakpoint else 0].add(clipped)  # add to the opposite breakpoint
        elif len(scores) == 1:
            clipped = scores[0][2]
            self.split_reads[1 if first_breakpoint else 0].add(clipped)  # add to the opposite breakpoint
        return True

    def decide_sequenced_strand(self, reads):
        """
        given a set of reads, determines the sequenced strand (if possible) and then returns the majority
        strand found

        Args:
            reads (set of :class:`pysam.AlignedSegment`): set of reads

        Returns:
            STRAND: the sequenced strand

        Raises:
            ValueError: input was an empty set or the ratio was not sufficient to decide on a strand
        """
        if not reads:
            raise ValueError('cannot determine the strand of a set of reads if the set is empty')

        strand_calls = {STRAND.POS: 0, STRAND.NEG: 0}
        for read in reads:
            try:
                strand = _read.sequenced_strand(read, self.strand_determining_read)
                strand_calls[strand] = strand_calls.get(strand, 0) + 1
            except ValueError:
                pass
        if sum(strand_calls.values()) == 0:
            raise ValueError('Could not determine strand. Insufficient mapped reads')
        if strand_calls[STRAND.POS] == 0:
            return STRAND.NEG
        elif strand_calls[STRAND.NEG] == 0:
            return STRAND.POS
        else:
            ratio = strand_calls[STRAND.POS] / (strand_calls[STRAND.NEG] + strand_calls[STRAND.POS])
            neg_ratio = 1 - ratio
            if ratio >= self.assembly_strand_concordance:
                return STRAND.POS
            elif neg_ratio >= self.assembly_strand_concordance:
                return STRAND.NEG
            raise ValueError('Could not determine the strand. Equivocal POS/(NEG + POS) ratio', ratio, strand_calls)

    def assemble_contig(self, log=DEVNULL):
        """
        uses the split reads and the partners of the half mapped reads to create a contig
        representing the sequence across the breakpoints

        if it is not strand specific then sequences are sorted alphanumerically and only the
        first of a pair is kept (paired by sequence)
        """
        # gather reads for the putative assembly
        assembly_sequences = {}
        # add split reads
        for read in list(itertools.chain.from_iterable(self.split_reads)) + list(self.spanning_reads):
            # ignore targeted realignments
            if read.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and read.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                continue
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

        # add half-mapped reads
        for read in itertools.chain.from_iterable(self.half_mapped):
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

        # add flanking reads
        for read, mate in self.flanking_pairs:
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

            assembly_sequences.setdefault(mate.query_sequence, set()).add(mate)
            rqs_comp = reverse_complement(mate.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(mate)

        log('assembly size of {} sequences'.format(len(assembly_sequences) // 2))

        kmer_size = self.read_length * self.assembly_kmer_size
        remap_min_overlap = max(self.read_length - self.assembly_min_exact_match_to_remap, kmer_size)

        contigs = assemble(
            assembly_sequences, kmer_size,
            min_edge_trim_weight=self.assembly_min_edge_trim_weight,
            assembly_max_paths=self.assembly_max_paths,
            min_contig_length=self.read_length,
            log=log,
            remap_min_overlap=remap_min_overlap,
            remap_min_exact_match=self.assembly_min_exact_match_to_remap,
            assembly_min_uniq=self.assembly_min_uniq,
            min_complexity=self.min_call_complexity
        )

        # add the input reads
        # drop any contigs without reads from both breakpoints
        filtered_contigs = []
        for ctg in contigs:
            for read_seq in ctg.remapped_sequences:
                ctg.input_reads.update(assembly_sequences[read_seq.query_sequence])
            break1_reads = {r.query_sequence for r in self.split_reads[0] | self.half_mapped[0] | self.spanning_reads}
            break2_reads = {r.query_sequence for r in self.split_reads[1] | self.half_mapped[1] | self.spanning_reads}
            for read, mate in self.flanking_pairs | self.compatible_flanking_pairs:
                break1_reads.add(read.query_sequence)
                break2_reads.add(mate.query_sequence)

            ctg_reads = {r.query_sequence for r in ctg.input_reads}
            ctg_reads.update({reverse_complement(r) for r in ctg_reads})
            if (ctg_reads & break1_reads and ctg_reads & break2_reads) or \
                    (not self.interchromosomal and len(self.break1 | self.break2) < self.read_length):
                filtered_contigs.append(ctg)
        log('filtered contigs from {} to {} based on remapped reads from both breakpoints'.format(len(contigs), len(filtered_contigs)), time_stamp=False)
        contigs = filtered_contigs

        # now determine the strand from the remapped reads if possible
        if self.stranded and self.bam_cache.stranded:  # strand specific
            for contig in contigs:
                build_strand = {STRAND.POS: 0, STRAND.NEG: 0}  # if neg will have to flip
                for read_seq in contig.remapped_sequences:
                    for read in assembly_sequences[read_seq.query_sequence]:
                        if read.is_unmapped:
                            continue
                        flip = False
                        if read.query_sequence != read_seq.query_sequence:
                            flip = not flip
                        try:
                            seq_strand = _read.sequenced_strand(read, self.strand_determining_read)
                            if seq_strand == STRAND.NEG:
                                flip = not flip
                            build_strand[STRAND.NEG if flip else STRAND.POS] += 1
                        except ValueError:
                            pass
                if sum(build_strand.values()) == 0:
                    continue
                elif build_strand[STRAND.POS] == 0:
                    flipped_build = True
                elif build_strand[STRAND.NEG] == 0:
                    flipped_build = False
                else:
                    ratio = build_strand[STRAND.POS] / (build_strand[STRAND.NEG] + build_strand[STRAND.POS])
                    neg_ratio = 1 - ratio
                    if ratio >= self.assembly_strand_concordance:
                        flipped_build = False
                    elif neg_ratio >= self.assembly_strand_concordance:
                        flipped_build = True
                    else:
                        continue
                if flipped_build:
                    contig.seq = reverse_complement(contig.seq)
                contig.strand_specific = True

        filtered_contigs = {}
        # sort so that the function is deterministic
        for contig in sorted(contigs, key=lambda x: (x.remap_score() * -1, x.seq)):
            # filter on evidence level
            if contig.remap_score() < self.assembly_min_remapped_seq or \
                    contig.remap_coverage() < self.assembly_min_remap_coverage:
                continue
            if self.stranded and self.bam_cache.stranded:
                filtered_contigs.setdefault(contig.seq, contig)
            else:
                rseq = reverse_complement(contig.seq)
                if contig.seq not in filtered_contigs and rseq not in filtered_contigs:
                    filtered_contigs[contig.seq] = contig
        self.contigs = list(filtered_contigs.values())

    def load_evidence(self, log=DEVNULL):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
        def cache_if_true(read):
            if read.is_unmapped or read.mate_is_unmapped:
                return True
            elif any([
                self.filter_secondary_alignments and read.is_secondary,
                read.mapping_quality < self.min_mapping_quality
            ]):
                return False
            elif set([x[0] for x in read.cigar]) & {CIGAR.S, CIGAR.H}:
                return True
            elif not read.is_proper_pair:
                if any([_read.orientation_supports_type(read, e) for e in self.putative_event_types()]):
                    return True
                elif self.compatible_type and _read.orientation_supports_type(read, self.compatible_type):
                    return True
            elif not self.interchromosomal and not self.opposing_strands:
                min_frag_est = abs(read.reference_start - read.next_reference_start) - self.read_length
                max_frag_est = min_frag_est + 3 * self.read_length
                if min_frag_est < self.min_expected_fragment_size or max_frag_est > self.max_expected_fragment_size:
                    return True

            return False

        def filter_if_true(read):
            if not cache_if_true(read):
                if any([
                    self.filter_secondary_alignments and read.is_secondary,
                    read.mapping_quality < self.min_mapping_quality
                ]):
                    return True
                elif not self.interchromosomal and set([x[0] for x in read.cigar]) & {CIGAR.I, CIGAR.D}:
                    return False
                return True
            return False

        flanking_pairs = set()  # collect putative pairs
        half_mapped_partners1 = set()
        half_mapped_partners2 = set()

        for read in self.bam_cache.fetch_from_bins(
                '{0}'.format(self.break1.chr),
                self.outer_window1[0],
                self.outer_window1[1],
                read_limit=self.fetch_reads_limit,
                sample_bins=self.fetch_reads_bins,
                min_bin_size=self.fetch_min_bin_size,
                cache=True,
                cache_if=cache_if_true,
                filter_if=filter_if_true):
            if read.mapping_quality < self.min_mapping_quality:
                continue
            self.counts[0] += 1
            if read.is_unmapped:
                continue
            if not self.collect_split_read(read, True):
                self.collect_spanning_read(read)
            if read.mate_is_unmapped:
                half_mapped_partners1.add(read)
            elif any([_read.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
                    (read.reference_id != read.next_reference_id) == self.interchromosomal:
                flanking_pairs.add(read)

        for read in self.bam_cache.fetch_from_bins(
                '{0}'.format(self.break2.chr),
                self.outer_window2[0],
                self.outer_window2[1],
                read_limit=self.fetch_reads_limit,
                sample_bins=self.fetch_reads_bins,
                min_bin_size=self.fetch_min_bin_size,
                cache=True,
                cache_if=cache_if_true,
                filter_if=filter_if_true):
            if read.mapping_quality < self.min_mapping_quality:
                continue

            self.counts[1] += 1

            if read.is_unmapped:
                continue
            if not self.collect_split_read(read, False):
                self.collect_spanning_read(read)
            if read.mate_is_unmapped:
                half_mapped_partners2.add(read)
            elif any([_read.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
                    (read.reference_id != read.next_reference_id) == self.interchromosomal:
                flanking_pairs.add(read)
        for flanking_read in sorted(flanking_pairs, key=lambda x: (x.query_name, x.reference_start)):
            # try and get the mate from the cache
            try:
                mates = self.bam_cache.get_mate(flanking_read, allow_file_access=False)
                for mate in mates:
                    if mate.is_unmapped:
                        log('ignoring unmapped mate', mate.query_name, level=logging.DEBUG)
                        continue
                    self.collect_flanking_pair(flanking_read, mate)
            except KeyError:
                pass

        if self.compatible_window1:
            compatible_type = SVTYPE.DUP
            if SVTYPE.DUP in self.putative_event_types():
                compatible_type = SVTYPE.INS

            compt_flanking = set()
            for read in self.bam_cache.fetch_from_bins(
                    '{0}'.format(self.break1.chr),
                    self.compatible_window1[0],
                    self.compatible_window1[1],
                    read_limit=self.fetch_reads_limit,
                    sample_bins=self.fetch_reads_bins,
                    min_bin_size=self.fetch_min_bin_size,
                    cache=True,
                    cache_if=cache_if_true,
                    filter_if=filter_if_true):
                if _read.orientation_supports_type(read, compatible_type):
                    compt_flanking.add(read)

            for read in self.bam_cache.fetch_from_bins(
                    '{0}'.format(self.break2.chr),
                    self.compatible_window2[0],
                    self.compatible_window2[1],
                    read_limit=self.fetch_reads_limit,
                    sample_bins=self.fetch_reads_bins,
                    min_bin_size=self.fetch_min_bin_size,
                    cache=True,
                    cache_if=cache_if_true,
                    filter_if=filter_if_true):
                if _read.orientation_supports_type(read, compatible_type):
                    compt_flanking.add(read)

            for flanking_read in compt_flanking:
                # try and get the mate from the cache
                try:
                    mates = self.bam_cache.get_mate(flanking_read, allow_file_access=False)
                    for mate in mates:
                        if mate.is_unmapped:
                            log('ignoring unmapped mate', mate.query_name, level=logging.DEBUG)
                            continue
                        try:
                            self.collect_compatible_flanking_pair(flanking_read, mate, compatible_type)
                        except ValueError:
                            pass
                except KeyError:
                    pass

        # now collect the half mapped reads
        log(
            'collected', len(half_mapped_partners1 | half_mapped_partners2),
            'putative half mapped reads', time_stamp=False)
        mates_found = 0
        for read in half_mapped_partners1 | half_mapped_partners2:
            # try and get the mate from the cache
            try:
                mates = self.bam_cache.get_mate(read, allow_file_access=False)
                mates_found += 1
                for mate in mates:
                    self.collect_half_mapped(read, mate)
            except KeyError:
                pass
        log(mates_found, 'half-mapped mates found')

    def copy(self):
        raise NotImplementedError('not appropriate for copy of evidence')

    def flatten(self):
        row = BreakpointPair.flatten(self)
        row.update({
            COLUMNS.raw_flanking_pairs: len(self.flanking_pairs),
            COLUMNS.raw_spanning_reads: len(self.spanning_reads),
            COLUMNS.raw_break1_split_reads: len(self.split_reads[0]),
            COLUMNS.raw_break2_split_reads: len(self.split_reads[1]),
            COLUMNS.raw_break1_half_mapped_reads: len(self.half_mapped[0]),
            COLUMNS.raw_break2_half_mapped_reads: len(self.half_mapped[1]),
            COLUMNS.protocol: self.protocol,
            COLUMNS.event_type: ';'.join(sorted(self.putative_event_types())),
            COLUMNS.contigs_assembled: len(self.contigs),
            COLUMNS.break1_ewindow: '{}-{}'.format(*self.outer_window1),
            COLUMNS.break2_ewindow: '{}-{}'.format(*self.outer_window2),
            COLUMNS.break1_ewindow_count: self.counts[0],
            COLUMNS.break2_ewindow_count: self.counts[1],
            COLUMNS.contigs_assembled: len(self.contigs)
        })
        return row

    def get_bed_repesentation(self):
        bed = []
        name = self.data.get(COLUMNS.cluster_id, None)
        bed.append((self.break1.chr, self.outer_window1[0] - 1, self.outer_window1[1], name))
        bed.append((self.break1.chr, self.inner_window1[0] - 1, self.inner_window1[1], name))
        bed.append((self.break2.chr, self.outer_window2[0] - 1, self.outer_window2[1], name))
        bed.append((self.break2.chr, self.inner_window2[0] - 1, self.inner_window2[1], name))
        return bed
