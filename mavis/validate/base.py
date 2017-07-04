import itertools
from ..constants import *
from ..error import *
from .constants import DEFAULTS
from ..assemble import assemble
from ..interval import Interval
from ..breakpoint import BreakpointPair
from ..bam import read as read_tools
from ..bam import cigar as cigar_tools
from ..bam.cache import BamCache
from ..util import devnull, log


class Evidence(BreakpointPair):
    @property
    def outer_window1(self):
        """:class:`~mavis.interval.Interval`: the window where evidence will be gathered for the first
        breakpoint

        see :ref:`theory - calculating the evidence window <theory-calculating-the-evidence-window>`
        """
        if self.collect_from_outer_window():
            try:
                return self.outer_windows[0]
            except AttributeError:
                raise NotImplementedError('abstract property must be overridden')
        else:
            return self.inner_window1

    @property
    def outer_window2(self):
        """:class:`~mavis.interval.Interval`: the window where evidence will be gathered for the second
        breakpoint

        see :ref:`theory - calculating the evidence window <theory-calculating-the-evidence-window>`
        """
        if self.collect_from_outer_window():
            try:
                return self.outer_windows[1]
            except AttributeError:
                raise NotImplementedError('abstract property must be overridden')
        else:
            return self.inner_window2

    @property
    def compatible_window1(self):
        """:class:`~mavis.interval.Interval`: the window/region where it is expected to 
        find reads in a compatible flanking pair (mate must be in compatible_window2)

        see :ref:`theory - calculating the evidence window <theory-calculating-the-evidence-window>`
        """
        return self.compatible_windows[0]

    @property
    def compatible_window2(self):
        """:class:`~mavis.interval.Interval`: the window/region where it is expected to 
        find reads in a compatible flanking pair (mate must be in compatible_window1)

        see :ref:`theory - calculating the evidence window <theory-calculating-the-evidence-window>`
        """
        return self.compatible_windows[1]

    @property
    def inner_window1(self):
        """:class:`~mavis.interval.Interval`: the window where evidence will be gathered for the first
        breakpoint
        """
        try:
            return self.inner_windows[0]
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

    @property
    def inner_window2(self):
        """:class:`~mavis.interval.Interval`: the window where evidence will be gathered for the second
        breakpoint
        """
        try:
            return self.inner_windows[1]
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

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
            REFERENCE_GENOME,
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
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`):
              dict of reference sequence by template/chr name
            data (dict): a dictionary of data to associate with the evidence object
            classification (SVTYPE): the event type
            protocol (PROTOCOL): genome or transcriptome
        """
        cls = self.__class__
        # initialize the breakpoint pair
        self.bam_cache = bam_cache
        if stranded and bam_cache.stranded:
            self.stranded = True
        else:
            self.stranded = False
        BreakpointPair.__init__(
            self, break1, break2,
            stranded=stranded,
            opposing_strands=opposing_strands,
            untemplated_seq=untemplated_seq,
            **data
        )
        d = dict()
        for arg in kwargs:
            if arg not in DEFAULTS.__dict__:
                raise AttributeError('unrecognized attribute', arg)
        d.update(DEFAULTS.__dict__)
        kwargs.setdefault('assembly_min_contig_length', int(read_length * 1.25))
        kwargs.setdefault('assembly_max_kmer_size', int(read_length * 0.7))
        d.update(kwargs)  # input arguments should override the defaults
        for arg, val in d.items():
            setattr(self, arg, val)

        self.bam_cache = bam_cache
        self.classification = classification
        self.REFERENCE_GENOME = REFERENCE_GENOME
        self.read_length = read_length
        self.stdev_fragment_size = stdev_fragment_size
        self.median_fragment_size = median_fragment_size
        self.compatible_windows = None

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
            self.compute_fragment_size(None)
        except NotImplementedError:
            raise NotImplementedError('abstract class cannot be initialized')
        except BaseException:
            pass

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
        else:
            return False

    def standardize_read(self, read):
        # recomputing to standardize b/c split reads can be used to call breakpoints exactly
        read.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
        # recalculate the read cigar string to ensure M is replaced with = or X
        c = cigar_tools.recompute_cigar_mismatch(
            read,
            self.REFERENCE_GENOME[self.bam_cache.get_read_reference_name(read)].seq
        )
        prefix = 0
        try:
            c, prefix = cigar_tools.extend_softclipping(c, self.sc_extension_stop)
        except AttributeError as err:
            pass
        read.cigar = cigar_tools.join(c)

        # makes sure all insertions are called as far 'right' as possible
        read.cigar = cigar_tools.hgvs_standardize_cigar(
            read, self.REFERENCE_GENOME[self.bam_cache.get_read_reference_name(read)].seq)
        read.reference_start = read.reference_start + prefix
        return read

    def putative_event_types(self):
        """
        Returns:
            list of :class:`~mavis.constants.SVTYPE`: list of the possible classifications
        """
        if self.classification:
            return [self.classification]
        else:
            return BreakpointPair.classify(self)

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
        for s in self.flanking_pairs:
            result.update(s)
        for s in self.split_reads:
            result.update(s)
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
            strand = read_tools.sequenced_strand(read, self.strand_determining_read)
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
        if not self.compatible_windows:
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
            strand1 = read_tools.sequenced_strand(read, self.strand_determining_read)
            strand2 = read_tools.sequenced_strand(mate, self.strand_determining_read)
            if strand1 != self.break1.strand or strand2 != self.break2.strand:
                return False

        # check that the pair orientation is correct
        if not read_tools.orientation_supports_type(read, compatible_type):
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
        if Interval.overlaps(iread, self.compatible_windows[0]) and \
                Interval.overlaps(imate, self.compatible_windows[1]):
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
            strand1 = read_tools.sequenced_strand(read, self.strand_determining_read)
            strand2 = read_tools.sequenced_strand(mate, self.strand_determining_read)

            if strand1 != self.break1.strand or strand2 != self.break2.strand:
                return False

        for event_type in self.putative_event_types():

            # check that the pair orientation is correct
            if not read_tools.orientation_supports_type(read, event_type):
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
            strand = read_tools.sequenced_strand(read, strand_determining_read=self.strand_determining_read)
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
        if cigar_tools.alignment_matches(read.cigar) >= self.min_sample_size_to_apply_percentage \
                and cigar_tools.match_percent(read.cigar) < self.min_anchor_match:
            return False  # too poor quality of an alignment
        if cigar_tools.longest_exact_match(read.cigar) < self.min_anchor_exact \
                and cigar_tools.longest_fuzzy_match(read.cigar, self.fuzzy_mismatch_number) < self.min_anchor_fuzzy:
            return False  # too poor quality of an alignment
        else:
            self.split_reads[0 if first_breakpoint else 1].add(read)

        # try mapping the soft-clipped portion to the other breakpoint
        w = (opposite_window[0], opposite_window[1])
        opposite_breakpoint_ref = self.REFERENCE_GENOME[opposite_breakpoint.chr].seq[w[0] - 1: w[1]]

        putative_alignments = None

        if not self.opposing_strands:  # same strand
            sc_align = read_tools.nsb_align(
                opposite_breakpoint_ref, read.query_sequence, min_consecutive_match=self.min_anchor_exact)

            for a in sc_align:
                a.flag = read.flag
            putative_alignments = sc_align
        else:
            # should align opposite the current read
            revcomp_sc_align = reverse_complement(read.query_sequence)
            revcomp_sc_align = read_tools.nsb_align(
                opposite_breakpoint_ref, revcomp_sc_align, min_consecutive_match=self.min_anchor_exact)

            for a in revcomp_sc_align:
                a.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
            putative_alignments = revcomp_sc_align

        scores = []

        for a in putative_alignments:  # loop over the alignments
            a.flag = a.flag | PYSAM_READ_FLAGS.SUPPLEMENTARY
            # set this flag so we don't recompute the cigar multiple
            a.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
            # add information from the original read
            a.reference_start = w[0] - 1 + a.reference_start
            a.reference_id = self.bam_cache.reference_id(opposite_breakpoint.chr)
            a.query_name = read.query_name
            a.set_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1, value_type='i')
            a.next_reference_start = read.next_reference_start
            a.next_reference_id = read.next_reference_id
            a.mapping_quality = NA_MAPPING_QUALITY
            try:
                cigar, offset = cigar_tools.extend_softclipping(
                    a.cigar, self.sc_extension_stop)
                a.cigar = cigar
                a.reference_start = a.reference_start + offset
            except AttributeError:
                # if the matches section is too small you can't extend the
                # softclipping
                pass
            s = cigar_tools.score(a.cigar)
            a.cigar = cigar_tools.join(a.cigar)
            if a.reference_id == a.next_reference_id:
                # https://samtools.github.io/hts-specs/SAMv1.pdf
                # unsigned observed template length equals the number of bases from the leftmost
                # mapped base to the rightmost mapped base
                tlen = abs(a.reference_start - a.next_reference_start) + 1
                if a.reference_start < a.next_reference_start:
                    a.template_length = tlen
                else:
                    a.template_length = -1 * tlen
            else:
                a.template_length = 0

            if cigar_tools.alignment_matches(a.cigar) >= self.min_sample_size_to_apply_percentage \
                    and cigar_tools.match_percent(a.cigar) < self.min_anchor_match:
                continue
            if cigar_tools.longest_exact_match(a.cigar) < self.min_anchor_exact \
                    and cigar_tools.longest_fuzzy_match(a.cigar, self.fuzzy_mismatch_number) < self.min_anchor_fuzzy:
                continue
            if self.max_sc_preceeding_anchor is not None:
                if opposite_breakpoint.orient == ORIENT.LEFT:
                    if a.cigar[0][0] == CIGAR.S and a.cigar[0][1] > self.max_sc_preceeding_anchor:
                        continue
                elif opposite_breakpoint.orient == ORIENT.RIGHT:
                    if a.cigar[-1][0] == CIGAR.S and a.cigar[-1][1] > self.max_sc_preceeding_anchor:
                        continue
            scores.append((s, cigar_tools.match_percent(a.cigar), a))

        scores = sorted(scores, key=lambda x: (x[0], x[1]), reverse=True) if len(scores) > 0 else []

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
        if len(reads) == 0:
            raise ValueError('cannot determine the strand of a set of reads if the set is empty')

        strand_calls = {STRAND.POS: 0, STRAND.NEG: 0}
        for read in reads:
            try:
                strand = read_tools.sequenced_strand(read, self.strand_determining_read)
                strand_calls[strand] = strand_calls.get(strand, 0) + 1
            except ValueError as err:
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

    def assemble_contig(self, log=devnull):
        """
        uses the split reads and the partners of the half mapped reads to create a contig
        representing the sequence across the breakpoints

        if it is not strand specific then sequences are sorted alphanumerically and only the
        first of a pair is kept (paired by sequence)
        """
        strand_specific = self.stranded and self.bam_cache.stranded
        # gather reads for the putative assembly
        assembly_sequences = {}
        targeted = 0
        # add split reads
        for r in list(itertools.chain.from_iterable(self.split_reads)) + list(self.spanning_reads):
            # ignore targeted realignments
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                targeted += 1
                continue
            assembly_sequences.setdefault(r.query_sequence, set()).add(r)
            rqs_comp = reverse_complement(r.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(r)

        # add half-mapped reads
        # exclude half-mapped reads if there is 'n' split reads that target align
        if targeted < self.assembly_min_tgt_to_exclude_half_map and \
                self.assembly_include_half_mapped_reads:
            for r in itertools.chain.from_iterable(self.half_mapped):
                assembly_sequences.setdefault(r.query_sequence, set()).add(r)
                rqs_comp = reverse_complement(r.query_sequence)
                assembly_sequences.setdefault(rqs_comp, set()).add(r)

        # add flanking reads
        if self.assembly_include_flanking_pairs:
            for read, mate in self.flanking_pairs:
                assembly_sequences.setdefault(read.query_sequence, set()).add(read)
                rqs_comp = reverse_complement(read.query_sequence)
                assembly_sequences.setdefault(rqs_comp, set()).add(read)

                assembly_sequences.setdefault(mate.query_sequence, set()).add(mate)
                rqs_comp = reverse_complement(mate.query_sequence)
                assembly_sequences.setdefault(rqs_comp, set()).add(mate)

        log('assembly size of {} sequences'.format(len(assembly_sequences) // 2))

        contigs = assemble(
            assembly_sequences,
            assembly_min_edge_weight=self.assembly_min_edge_weight,
            assembly_min_nc_edge_weight=self.assembly_min_nc_edge_weight,
            assembly_max_paths=self.assembly_max_paths,
            log=log,
            assembly_min_exact_match_to_remap=self.assembly_min_exact_match_to_remap,
            assembly_min_contig_length=self.assembly_min_contig_length,
            assembly_max_kmer_size=self.assembly_max_kmer_size,
            assembly_max_kmer_strict=self.assembly_max_kmer_strict
        )

        # add the input reads
        for ctg in contigs:
            for read_seq in ctg.remapped_sequences:
                ctg.input_reads.update(assembly_sequences[read_seq.query_sequence])

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
                            seq_strand = read_tools.sequenced_strand(read, self.strand_determining_read)
                            if seq_strand == STRAND.NEG:
                                flip = not flip
                            build_strand[STRAND.NEG if flip else STRAND.POS] += 1
                        except ValueError as err:
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
        for c in sorted(contigs, key=lambda x: (x.remap_score() * -1, x.seq)):
            if c.remap_score() < self.assembly_min_remapped_seq:  # filter on evidence level
                continue
            if self.stranded and self.bam_cache.stranded:
                filtered_contigs.setdefault(c.seq, c)
            else:
                rseq = reverse_complement(c.seq)
                if c.seq not in filtered_contigs and rseq not in filtered_contigs:
                    filtered_contigs[c.seq] = c
        self.contigs = list(filtered_contigs.values())

    @classmethod
    def load_multiple(cls, evidence, log=devnull):
        """
        loads evidence from the bam file for multiple evidence objects at once
        
        Args:
            evidence (list of :class:`Evidence`): list of evidence objects to collect evidence for

        Warning:
            this is not exactly equivalent to multiple calls of load_evidence because it
            does not keep a running total of reads between bins and cannot adjust dynamically.
            Effectively this means load_evidence may potentially collect more evidence
        """
        log('gathering evidence for {} events'.format(len(evidence)))
        cache = evidence[0].bam_cache
        filter_secondary_alignments = evidence[0].filter_secondary_alignments
        min_mapping_quality = evidence[0].min_mapping_quality
        max_expected_fragment_size = evidence[0].max_expected_fragment_size
        min_expected_fragment_size = evidence[0].min_expected_fragment_size
        read_length = evidence[0].read_length
        fetch_min_bin_size = evidence[0].fetch_min_bin_size
        protocol = evidence[0].protocol

        # check the inputs make sense
        for ev in evidence:
            if cache is not ev.bam_cache:
                raise UserWarning('cannot load_multiple when bam_cache is not the same object')
            if filter_secondary_alignments != ev.filter_secondary_alignments:
                raise UserWarning('cannot load_multiple when filter_secondary_alignments differs')
            if min_mapping_quality != ev.min_mapping_quality:
                raise UserWarning('cannot load_multiple when the min_mapping_quality differs')
            if read_length != ev.read_length:
                raise UserWarning('cannot load_multiple when the read_length differs')
            max_expected_fragment_size = min([ev.max_expected_fragment_size, max_expected_fragment_size])
            min_expected_fragment_size = max([ev.min_expected_fragment_size, min_expected_fragment_size])
            fetch_min_bin_size = max([ev.fetch_min_bin_size, fetch_min_bin_size])
            if ev.protocol == PROTOCOL.TRANS:
                protocol = PROTOCOL.TRANS

        def cache_if_true(read):
            if read.is_unmapped or read.mate_is_unmapped:
                return True
            elif any([
                filter_secondary_alignments and read.is_secondary,
                read.mapping_quality < min_mapping_quality
            ]):
                return False
            elif set([x[0] for x in read.cigar]) & {CIGAR.D, CIGAR.N, CIGAR.S, CIGAR.I}:
                return True
            elif read.is_proper_pair and protocol != PROTOCOL.TRANS:
                min_frag_est = abs(read.reference_start - read.next_reference_start) - read_length
                max_frag_est = min_frag_est + 3 * read_length
                if min_frag_est >= min_expected_fragment_size and max_frag_est <= max_expected_fragment_size:
                    return False
            return True

        def filter_if_true(read):
            if not cache_if_true(read):
                return True
            elif read.is_unmapped:
                return True
            return False

        # create the intervals to investigate
        read_limits = {}  # chr => interval => cap
        for ev in evidence:
            # first breakpoint window
            bins = BamCache._generate_fetch_bins(
                ev.outer_window1[0], ev.outer_window1[1], ev.fetch_reads_bins, fetch_min_bin_size)
            read_cap = ev.fetch_reads_limit // len(bins)
            for b in bins:
                read_limits.setdefault(ev.break1.chr, {})
                read_limits[ev.break1.chr][b] = max([read_limits[ev.break1.chr].get(b, 0), read_cap])
            # second breakpoint window
            bins = BamCache._generate_fetch_bins(
                ev.outer_window2[0], ev.outer_window2[1], ev.fetch_reads_bins, fetch_min_bin_size)
            read_cap = ev.fetch_reads_limit // len(bins)
            for b in bins:
                read_limits.setdefault(ev.break2.chr, {})
                read_limits[ev.break2.chr][b] = max([read_limits[ev.break2.chr].get(b, 0), read_cap])
            # compatible_windows
            if ev.compatible_windows and ev.collect_from_outer_window():
                # first compatible breakpoint window
                bins = BamCache._generate_fetch_bins(
                    ev.compatible_window1[0], ev.compatible_window1[1], ev.fetch_reads_bins, fetch_min_bin_size)
                read_cap = ev.fetch_reads_limit // len(bins)
                for b in bins:
                    read_limits.setdefault(ev.break1.chr, {})
                    read_limits[ev.break1.chr][b] = max([read_limits[ev.break1.chr].get(b, 0), read_cap])
                # second compatible breakpoint window
                bins = BamCache._generate_fetch_bins(
                    ev.compatible_window2[0], ev.compatible_window2[1], ev.fetch_reads_bins, fetch_min_bin_size)
                read_cap = ev.fetch_reads_limit // len(bins)
                for b in bins:
                    read_limits.setdefault(ev.break2.chr, {})
                    read_limits[ev.break2.chr][b] = max([read_limits[ev.break2.chr].get(b, 0), read_cap])
        # get the reads for any given interval, map over each applicable evidence
        fetch_regions = []

        for chr in read_limits:
            weighted_intervals = Interval.split_overlap(*read_limits[chr], weight_mapping=read_limits[chr])
            # Merge any intervals that are < min size?
            temp = {}
            intervals = sorted(weighted_intervals)
            i = 0
            while i < len(intervals):
                curr = intervals[i]
                if len(curr) < fetch_min_bin_size:
                    # merge with the following interval
                    merge_itvl = curr
                    merge_weight = weighted_intervals[curr]
                    i += 1
                    while len(merge_itvl) < fetch_min_bin_size and i < len(intervals):
                        merge_itvl = intervals[i] | merge_itvl
                        merge_weight += weighted_intervals[intervals[i]]
                        i += 1
                    temp[merge_itvl] = merge_weight
                else:
                    temp[curr] = weighted_intervals[curr]
                    i += 1
            weighted_intervals = temp
            for bin, limit in weighted_intervals.items():
                fetch_regions.append((chr, bin.start, bin.end, limit))

        putative_half_maps = set()
        putative_flanking = set()
        fetch_regions = sorted(fetch_regions, reverse=True)

        for i in range(0, len(fetch_regions)):
            chr, start, end, limit = fetch_regions[i]
            log('({} of {}) loading the bin {}:{}-{} (limit {}; size {})'.format(
                i + 1, len(fetch_regions), chr, start, end, limit, end - start + 1))
            for read in cache.fetch(
                chr, start, end,
                limit=limit,
                cache_if=cache_if_true,
                filter_if=filter_if_true,
                stop_on_cached_read=True
            ):
                read_itvl = Interval(read.reference_start + 1, read.reference_end)
                if read.mate_is_unmapped:
                    putative_half_maps.add(read)
                    continue
                # breakpoint region 1
                for ev in [e for e in evidence if e.break1.chr == chr]:
                    compatible_type = SVTYPE.INS if SVTYPE.DUP in ev.putative_event_types() else SVTYPE.DUP
                    if not Interval.overlaps(read_itvl, ev.outer_window1):
                        continue
                    ev.counts[0] += 1
                    if not ev.collect_split_read(read, first_breakpoint=True):
                        ev.collect_spanning_read(read)
                    try:
                        mates = cache.get_mate(read, allow_file_access=False)
                        for mate in mates:
                            ev.collect_flanking_pair(read, mate)
                            if ev.compatible_windows:
                                ev.collect_compatible_flanking_pair(read, mate, compatible_type)
                    except KeyError:
                        putative_flanking.add(read)
                # breakpoint region 2
                for ev in [e for e in evidence if e.break2.chr == chr]:
                    compatible_type = SVTYPE.INS if SVTYPE.DUP in ev.putative_event_types() else SVTYPE.DUP
                    if not Interval.overlaps(read_itvl, ev.outer_window2):
                        continue
                    ev.counts[1] += 1
                    ev.collect_split_read(read, first_breakpoint=False)
                    try:
                        mates = cache.get_mate(read, allow_file_access=False)
                        for mate in mates:
                            ev.collect_flanking_pair(read, mate)
                            if ev.compatible_windows:
                                ev.collect_compatible_flanking_pair(read, mate, compatible_type)
                    except KeyError:
                        putative_flanking.add(read)
        # go over the flanking /half-mapped that have uncached mates previously
        # if we didn't go over them yet, then they are outside the target region and should be ignored
        log('gathering cached flanking pairs', len(putative_flanking))
        for read in putative_flanking:
            try:
                mates = cache.get_mate(read, allow_file_access=False)
                for mate in mates:
                    for ev in evidence:
                        ev.collect_flanking_pair(read, mate)
                        if ev.compatible_windows:
                            compatible_type = SVTYPE.INS if SVTYPE.DUP in ev.putative_event_types() else SVTYPE.DUP
                            ev.collect_compatible_flanking_pair(read, mate, compatible_type)
            except KeyError:
                pass
        log('gathering cached half-mapped reads', len(putative_half_maps))
        mates_found = 0
        for read in putative_half_maps:
            try:
                mates = cache.get_mate(read, allow_file_access=False)
                mates_found += 1
                for mate in mates:
                    for ev in evidence:
                        ev.collect_half_mapped(read, mate)
            except KeyError:
                pass
        log(mates_found, 'half-mapped mates found')

    def load_evidence(self, log=devnull):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering
        """
        max_dist = max(
            len(Interval.union(self.break1, self.break2)),
            len(self.untemplated_seq if self.untemplated_seq else '')
        )
        
        def cache_if_true(read):
            if read.is_unmapped or read.mate_is_unmapped:
                return True
            elif any([
                self.filter_secondary_alignments and read.is_secondary,
                read.mapping_quality < self.min_mapping_quality
            ]):
                return False
            elif set([x[0] for x in read.cigar]) & {CIGAR.D, CIGAR.N, CIGAR.S, CIGAR.I}:
                return True
            elif read.is_proper_pair and self.protocol != PROTOCOL.TRANS:
                min_frag_est = abs(read.reference_start - read.next_reference_start) - self.read_length
                max_frag_est = min_frag_est + 3 * self.read_length
                if min_frag_est >= self.min_expected_fragment_size and max_frag_est <= self.max_expected_fragment_size:
                    return False
            return True

        def filter_if_true(read):
            if not cache_if_true(read):
                return True
            elif read.is_unmapped:
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
            elif any([read_tools.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
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
            elif any([read_tools.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
                    (read.reference_id != read.next_reference_id) == self.interchromosomal:
                flanking_pairs.add(read)
        for fl in sorted(flanking_pairs, key=lambda x: (x.query_name, x.reference_start)):
            # try and get the mate from the cache
            try:
                mates = self.bam_cache.get_mate(fl, allow_file_access=False)
                for mate in mates:
                    self.collect_flanking_pair(fl, mate)
            except KeyError:
                pass

        if self.compatible_windows:
            compatible_type = SVTYPE.DUP
            if SVTYPE.DUP in self.putative_event_types():
                compatible_type = SVTYPE.INS

            compt_flanking = set()
            for read in self.bam_cache.fetch_from_bins(
                    '{0}'.format(self.break1.chr),
                    self.compatible_windows[0][0],
                    self.compatible_windows[0][1],
                    read_limit=self.fetch_reads_limit,
                    sample_bins=self.fetch_reads_bins,
                    min_bin_size=self.fetch_min_bin_size,
                    cache=True,
                    cache_if=cache_if_true,
                    filter_if=filter_if_true):
                if read_tools.orientation_supports_type(read, compatible_type):
                    compt_flanking.add(read)

            for read in self.bam_cache.fetch_from_bins(
                    '{0}'.format(self.break2.chr),
                    self.compatible_windows[1][0],
                    self.compatible_windows[1][1],
                    read_limit=self.fetch_reads_limit,
                    sample_bins=self.fetch_reads_bins,
                    min_bin_size=self.fetch_min_bin_size,
                    cache=True,
                    cache_if=cache_if_true,
                    filter_if=filter_if_true):
                if read_tools.orientation_supports_type(read, compatible_type):
                    compt_flanking.add(read)

            for fl in compt_flanking:
                # try and get the mate from the cache
                try:
                    mates = self.bam_cache.get_mate(fl, allow_file_access=False)
                    for mate in mates:
                        try:
                            self.collect_compatible_flanking_pair(fl, mate, compatible_type)
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
            COLUMNS.contigs_aligned: sum([len(c.alignments) for c in self.contigs]),
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
