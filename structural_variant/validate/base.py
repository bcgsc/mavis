import itertools
from ..constants import *
from ..error import *
from ..assemble import assemble
from ..interval import Interval
from ..breakpoint import BreakpointPair
from ..bam import read as read_tools
from ..bam import cigar as cigar_tools


class Evidence(BreakpointPair):
    @property
    def window1(self):
        """(:class:`~structural_variant.interval.Interval`): the window where evidence will be gathered for the first
        breakpoint
        """
        try:
            return self.windows[0]
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

    @property
    def window2(self):
        """(:class:`~structural_variant.interval.Interval`): the window where evidence will be gathered for the second
        breakpoint
        """
        try:
            return self.windows[1]
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

    @property
    def min_expected_fragment_size(self):
        try:
            return self._min_expected_fragment_size
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

    @property
    def max_expected_fragment_size(self):
        try:
            return self._max_expected_fragment_size
        except AttributeError:
            raise NotImplementedError('abstract property must be overridden')

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
            untemplated_sequence=None,
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
        # initialize the breakpoint pair
        BreakpointPair.__init__(
            self, break1, break2,
            stranded=stranded,
            opposing_strands=opposing_strands,
            untemplated_sequence=untemplated_sequence,
            data=data
        )
        d = dict()
        d.update(VALIDATION_DEFAULTS.__dict__)
        d.update(kwargs)  # input arguments should override the defaults
        for arg, val in d.items():
            setattr(self, arg, val)

        self.bam_cache = bam_cache
        self.classification = classification
        self.REFERENCE_GENOME = REFERENCE_GENOME
        self.read_length = read_length
        self.stdev_fragment_size = stdev_fragment_size
        self.median_fragment_size = median_fragment_size

        if self.classification is not None and self.classification not in BreakpointPair.classify(self):
            raise AttributeError(
                'breakpoint pair improper classification', BreakpointPair.classify(self), self.classification)

        if self.break1.orient == ORIENT.NS or self.break2.orient == ORIENT.NS:
            raise NotSpecifiedError(
                'input breakpoint pair must specify strand and orientation. Cannot be \'not specified'
                '\' for evidence gathering')

        self.split_reads = (set(), set())
        self.flanking_pairs = set()
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
        except:
            pass
    
    def putative_event_types(self):
        if self.classification:
            return [self.classification]
        else:
            return BreakpointPair.classify(self)

    def compute_fragment_size(self, read, mate):
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

    def add_spanning_read(self, read):  # TODO
        """
        spanning read: a read covering BOTH breakpoints

        this is only applicable to small events

        .. todo::
            add support for indels
        """
        # check that the read fully covers BOTH breakpoints
        read_start = read.reference_start + 1 - \
            self.call_error  # adjust b/c pysam is 0-indexed
        # don't adjust b/c pysam is 0-indexed but end coord are one past
        read_end = read.reference_end + self.call_error
        if self.break1.start >= read_start and self.break1.end <= read_end  \
                and self.break2.start >= read_start and self.break2.end <= read_end:
            pass  # now check if this supports the putative event types
            raise NotImplementedError('have not added support for indels yet')
        else:
            raise UserWarning(
                'this does not cover/span both breakpoints and cannot be added as spanning evidence')

    def add_flanking_pair(self, read, mate):
        """
        checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
        this event

        Args:
            read (pysam.AlignedSegment): the read to add
        Raises:
            UserWarning: the read does not support this event or does not pass quality filters
        """
        if read.is_unmapped or mate.is_unmapped or read.query_name != mate.query_name or read.is_read1 == mate.is_read1:
            raise ValueError('input reads must be a mapped and mated pair')
        # check that the references are right for the pair
        if self.interchromosomal:
            if read.reference_id == read.next_reference_id:
                return False
        elif read.reference_id != read.next_reference_id:
            return False

        # order the read pairs so that they are in the same order that we expect for the breakpoints
        if read.reference_id != mate.reference_id:
            if self.bam_cache.chr(read) > self.bam_cache.chr(mate):
                read, mate = mate, read
        elif read.reference_start > mate.reference_start:
            read, mate = mate, read

        if self.bam_cache.chr(read) != self.break1.chr or self.bam_cache.chr(mate) != self.break2.chr:
            return False

        # check if this read falls in the first breakpoint window
        is_stranded = self.bam_cache.stranded and self.stranded
        iread = Interval(read.reference_start + 1, read.reference_end)
        imate = Interval(mate.reference_start + 1, mate.reference_end)
        added = False
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
            if Interval.overlaps(iread, self.window1) and Interval.overlaps(imate, self.window2):
                if not is_stranded or self.read_pair_strand(read) == self.break1.strand:
                    self.flanking_pairs.add((read, mate))
                    added = True
        return added

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

        # check if the read falls within the evidence collection window
        # recall that pysam is 0-indexed and the window is 1-indexed
        if read.reference_start > window.end - 1 or read.reference_end < window.start \
                or self.bam_cache.chr(read) != breakpoint.chr:
            return False  # read not in breakpoint evidence window
        # can only enforce strand if both the breakpoint and the bam are stranded
        if self.stranded and self.bam_cache.stranded:
            if read_tools.read_pair_strand(read) != (breakpoint.strand == STRAND.NEG):
                return False  # split read not on the appropriate strand
        
        primary = ''
        clipped = ''
        if breakpoint.orient == ORIENT.LEFT:
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
            # end is exclusive in pysam
            clipped = read.query_sequence[read.query_alignment_end:]
        elif breakpoint.orient == ORIENT.RIGHT:
            clipped = read.query_sequence[:read.query_alignment_start]
            primary = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
        else:
            raise NotSpecifiedError(
                'cannot assign split reads to a breakpoint where the orientation has not been specified')
        assert(len(primary) + len(clipped) == len(read.query_sequence))
        if len(primary) < self.min_anchor_exact or len(clipped) < self.min_anchor_exact:
            # split read does not meet the minimum anchor criteria
            return False
        elif len(read.query_sequence) - (read.query_alignment_end + 2) < self.min_anchor_exact \
                and (read.query_alignment_start + 1) < self.min_anchor_exact:
            # split read does not meet the minimum anchor criteria
            return False
        elif len(primary) < self.min_anchor_exact or len(clipped) < self.min_anchor_exact:
            # split read does not meet the minimum anchor criteria
            return False

        if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR):
            # recomputing to standardize b/c split reads can be used to call breakpoints exactly
            read.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
            # recalculate the read cigar string to ensure M is replaced with = or X
            c = cigar_tools.recompute_cigar_mismatch(
                read,
                self.REFERENCE_GENOME[self.bam_cache.chr(read)].seq
            )
            prefix = 0
            try:
                c, prefix = cigar_tools.extend_softclipping(c, self.sc_extension_stop)
            except AttributeError:
                pass
            read.cigar = c
            read.reference_start = read.reference_start + prefix
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
            sc_align = read_tools.nsb_align(opposite_breakpoint_ref, read.query_sequence)

            for a in sc_align:
                a.flag = read.flag
            putative_alignments = sc_align
        else:
            # should align opposite the current read
            revcomp_sc_align = reverse_complement(read.query_sequence)
            revcomp_sc_align = read_tools.nsb_align(opposite_breakpoint_ref, revcomp_sc_align)

            for a in revcomp_sc_align:
                a.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
            putative_alignments = revcomp_sc_align

        scores = []

        for a in putative_alignments:  # loop over the alignments
            a.flag = a.flag | PYSAM_READ_FLAGS.SECONDARY  # this is a secondary alignment
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

    def assemble_split_reads(self, log=lambda *x: None):
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
        for r in itertools.chain.from_iterable(self.split_reads):
            # ignore targeted realignments
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                targeted += 1
                continue
            assembly_sequences.setdefault(r.query_sequence, set()).add(r)
            if self.opposing_strands:
                rqs_comp = reverse_complement(r.query_sequence)
                assembly_sequences.setdefault(rqs_comp, set()).add(r)

        # add half-mapped reads
        # exclude half-mapped reads if there is 'n' split reads that target align
        if targeted < self.assembly_min_tgt_to_exclude_half_map and \
                self.assembly_include_half_mapped_reads:
            for r in itertools.chain.from_iterable(self.half_mapped):
                assembly_sequences.setdefault(r.query_sequence, set()).add(r)
                if self.opposing_strands:
                    rqs_comp = reverse_complement(r.query_sequence)
                    assembly_sequences.setdefault(rqs_comp, set()).add(r)
                try:
                    for m in self.bam_cache.get_mate(r):
                        if not m.is_unmapped:
                            continue
                        assembly_sequences.setdefault(m.query_sequence, set()).add(m)
                        if self.opposing_strands:
                            rqs_comp = reverse_complement(m.query_sequence)
                            assembly_sequences.setdefault(rqs_comp, set()).add(m)
                except KeyError:
                    pass

        # add flanking reads
        if self.assembly_include_flanking_pairs:
            for r in itertools.chain.from_iterable(self.flanking_pairs):
                if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                    continue
                assembly_sequences.setdefault(r.query_sequence, set()).add(r)
                if self.opposing_strands:
                    rqs_comp = reverse_complement(r.query_sequence)
                    assembly_sequences.setdefault(rqs_comp, set()).add(r)

        log('assembly size of {} sequences'.format(len(assembly_sequences) // 2))
        contigs = assemble(
            assembly_sequences,
            assembly_min_edge_weight=self.assembly_min_edge_weight,
            assembly_max_paths=self.assembly_max_paths,
            log=log,
            assembly_min_consec_match_remap=self.min_anchor_exact,
            assembly_min_contig_length=self.read_length
        )

        # now determine the strand from the remapped reads if possible
        if self.stranded and self.bam_cache.stranded:  # strand specific
            for contig in contigs:
                if len(contig.remapped_reads.keys()) == 0:
                    continue
                strand_calls = {STRAND.POS: 0, STRAND.NEG: 0}
                for seq in contig.remapped_reads:
                    for read in assembly_sequences[seq.query_sequence]:
                        if read.is_unmapped:
                            continue
                        if read.is_read1:
                            if read.is_reverse:
                                strand_calls[STRAND.NEG] += 1
                            else:
                                strand_calls[STRAND.POS] += 1
                        else:
                            if read.is_reverse:
                                strand_calls[STRAND.POS] += 1
                            else:
                                strand_calls[STRAND.NEG] += 1

                if strand_calls[STRAND.NEG] == 0 and strand_calls[STRAND.POS] == 0:
                    raise AssertionError('could not determine strand', strand_calls)
                elif strand_calls[STRAND.POS] > strand_calls[STRAND.NEG]:
                    ratio = strand_calls[STRAND.POS] / (strand_calls[STRAND.NEG] + strand_calls[STRAND.POS])
                    if ratio >= self.assembly_strand_concordance:
                        continue
                elif strand_calls[STRAND.NEG] > strand_calls[STRAND.POS]:
                    ratio = strand_calls[STRAND.NEG] / (strand_calls[STRAND.NEG] + strand_calls[STRAND.POS])
                    if ratio >= self.assembly_strand_concordance:
                        # negative strand. should reverse
                        contig.sequence = reverse_complement(contig.sequence)
                        continue
                reverse = True if strand_calls[STRAND.NEG] > strand_calls[STRAND.POS] else False
                for seq in contig.remapped_reads:
                    for read in assembly_sequences[seq.query_sequence]:
                        if read.is_reverse != reverse:
                            print('conflicting read', read.cigar, read.query_sequence)
                print(contig.sequence)
                raise AssertionError('could not determine contig strand', strand_calls)

        filtered_contigs = {}
        # sort so that the function is deterministic
        for c in sorted(contigs, key=lambda x: (x.remap_score() * -1, x.sequence)):
            if c.remap_score() < self.assembly_min_remap:  # filter on evidence level
                continue
            if not self.stranded or not self.bam_cache.stranded:  # not strand specific
                rseq = reverse_complement(c.sequence)
                if c.sequence not in filtered_contigs and (strand_specific or rseq not in filtered_contigs):
                    filtered_contigs[c.sequence] = c
            else:
                if c.sequence not in filtered_contigs:
                    filtered_contigs[c.sequence] = c
        self.contigs = list(filtered_contigs.values())

    def load_evidence(self, grab_unmapped_partners=True):
        """
        open the associated bam file and read and store the evidence
        does some preliminary read-quality filtering

        .. todo::
            support gathering evidence for small structural variants
        """
        bin_gap_size = self.read_length // 2

        max_dist = max(
            len(Interval.union(self.break1, self.break2)),
            len(self.untemplated_sequence if self.untemplated_sequence else '')
        )

        if not self.interchromosomal and max_dist < self.stdev_fragment_size * self.stdev_count_abnormal:
            raise NotImplementedError('evidence gathering for small structural variants is not supported')
            # needs special consideration b/c won't have flanking reads and may have spanning reads
        
        def filter_if_true(read):
            if self.filter_secondary_alignments and read.is_secondary:
                return True
            elif read.mapping_quality < self.min_mapping_quality:
                return True
            return False
        
        flanking_pairs = []  # collect putative pairs

        for read in self.bam_cache.fetch(
                '{0}'.format(self.break1.chr),
                self.window1[0],
                self.window1[1],
                read_limit=self.fetch_reads_limit,
                sample_bins=self.fetch_reads_bins,
                bin_gap_size=bin_gap_size,
                cache=True,
                filter_if=filter_if_true):
            self.counts[0] += 1
            if read.is_unmapped:
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
            elif any([read_tools.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
                    (read.reference_id != read.next_reference_id) == self.interchromosomal:
                flanking_pairs.append(read)
        print('flanking_pairs', len(flanking_pairs))
        for read in self.bam_cache.fetch(
                '{0}'.format(self.break2.chr),
                self.window2[0],
                self.window2[1],
                read_limit=self.fetch_reads_limit,
                sample_bins=self.fetch_reads_bins,
                bin_gap_size=bin_gap_size,
                cache=True,
                filter_if=filter_if_true):
            self.counts[1] += 1

            if read.is_unmapped:
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
            elif any([read_tools.orientation_supports_type(read, et) for et in self.putative_event_types()]) and \
                    (read.reference_id != read.next_reference_id) == self.interchromosomal:
                flanking_pairs.append(read)
        
        added = set()
        for fl in flanking_pairs:
            # try and get the mate from the cache
            try:
                mates = self.bam_cache.get_mate(fl, allow_file_access=False)
                for mate in mates:
                    self.add_flanking_pair(fl, mate)
            except KeyError:
                pass
        print('counts', self.counts)

    def copy(self):
        raise NotImplementedError('not appropriate for copy of evidence')

    def flatten(self):
        raise NotImplementedError('TODO')
