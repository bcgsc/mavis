from ..breakpoint import BreakpointPair, Breakpoint
from ..constants import CALL_METHOD, SVTYPE, PYSAM_READ_FLAGS, ORIENT, PROTOCOL, COLUMNS, STRAND
from ..bam import read as read_tools
from ..interval import Interval
from .evidence import TranscriptomeEvidence
import itertools
import statistics
import math


class EventCall(BreakpointPair):
    """
    class for holding evidence and the related calls since we can't freeze the evidence object
    directly without a lot of copying. Instead we use call objects which are basically
    just a reference to the evidence object and decisions on class, exact breakpoints, etc
    """
    
    @property
    def has_compatible(self):
        """bool: True if compatible flanking pairs are appropriate to collect"""
        try:
            self.compatible_type
            return True
        except AttributeError:
            return False

    def __init__(
        self,
        b1, b2,
        source_evidence,
        event_type,
        call_method,
        break2_call_method=None,
        contig=None,
        contig_alignment=None,
        untemplated_seq=None
    ):
        """
        Args:
            evidence (Evidence): the evidence object we are calling based on
            event_type (SVTYPE): the type of structural variant
            breakpoint_pair (BreakpointPair): the breakpoint pair representing the exact breakpoints
            call_method (CALL_METHOD): the way the breakpoints were called
            contig (Contig): the contig used to call the breakpoints (if applicable)
        """
        if untemplated_seq is None:
            untemplated_seq = source_evidence.untemplated_seq
        if break2_call_method is None:
            break2_call_method = call_method

        BreakpointPair.__init__(
            self, b1, b2,
            stranded=source_evidence.stranded and source_evidence.bam_cache.stranded,
            opposing_strands=source_evidence.opposing_strands,
            untemplated_seq=untemplated_seq
        )
        self.data.update(source_evidence.data)
        if not source_evidence.bam_cache.stranded:
            self.break1.strand = STRAND.NS
            self.break2.strand = STRAND.NS
        self.source_evidence = source_evidence
        self.spanning_reads = set()
        self.flanking_pairs = set()
        self.break1_split_reads = set()
        self.break2_split_reads = set()
        self.compatible_flanking_pairs = set()
        # check that the event type is compatible
        self.event_type = SVTYPE.enforce(event_type)
        if event_type == SVTYPE.DUP:
            self.compatible_type = SVTYPE.INS
        elif event_type == SVTYPE.INS:
            self.compatible_type = SVTYPE.DUP
        if event_type not in BreakpointPair.classify(self):
            raise ValueError(
                'event_type is not compatible with the breakpoint call', event_type, BreakpointPair.classify(self))
        self.contig = contig
        self.call_method = (CALL_METHOD.enforce(call_method), CALL_METHOD.enforce(break2_call_method))
        if contig and self.call_method != (CALL_METHOD.CONTIG, CALL_METHOD.CONTIG):
            raise ValueError('if a contig is given the call method must be by contig')
        self.contig_alignment = contig_alignment

    def support(self):
        support = set()
        support.update(self.spanning_reads)
        for read, mate in self.flanking_pairs | self.compatible_flanking_pairs:
            support.add(read)
            support.add(mate)
        support.update(self.break1_split_reads)
        support.update(self.break2_split_reads)
        if self.contig:
            support.update(self.contig.input_reads)
        return support

    def add_flanking_support(self, flanking_pairs, is_compatible=False):
        """
        counts the flanking read-pair support for the event called. The original source evidence may
        have contained evidence for multiple events and uses a larger range so flanking pairs here
        are checked specifically against the current breakpoint call

        Returns:
            tuple:
            * :class:`set` of :class:`str` - set of the read query_names
            * :class:`int` - the median insert size
            * :class:`int` - the standard deviation (from the median) of the insert size

        see :ref:`theory - determining flanking support <theory-determining-flanking-support>`
        """
        support = set()

        fragment_sizes = []

        min_frag = max([
            self.source_evidence.min_expected_fragment_size + Interval.dist(self.break1, self.break2),
            self.source_evidence.max_expected_fragment_size])
        max_frag = len(self.break1 | self.break2) + self.source_evidence.max_expected_fragment_size

        encompass = len(self.break1 | self.break2)

        for read, mate in flanking_pairs:
            # check that the fragment size is reasonable
            fragment_size = self.source_evidence.compute_fragment_size(read, mate)
            if self.event_type == SVTYPE.DEL:
                if fragment_size.end < min_frag or fragment_size.start > max_frag:
                    continue
            elif self.event_type == SVTYPE.INS:
                if fragment_size.start >= self.source_evidence.min_expected_fragment_size:
                    continue
            if self.interchromosomal != (read.reference_id != mate.reference_id):
                continue
            # check that the flanking reads work with the current call
            if not read_tools.orientation_supports_type(
                    read, self.event_type if not is_compatible else self.compatible_type):
                continue
            # check that the positions make sense
            LEFT = ORIENT.LEFT if not is_compatible else ORIENT.RIGHT
            if self.break1.orient == LEFT:
                if self.break2.orient == LEFT:  # L L
                    if not all([
                        read.reference_start + 1 <= self.break1.end,
                        mate.reference_start + 1 <= self.break2.end,
                        mate.reference_end > self.break1.start
                    ]):
                        continue
                else:  # L R
                    if not all([
                        read.reference_start + 1 <= self.break1.end,
                        mate.reference_end >= self.break2.start
                    ]):
                        continue
            else:
                if self.break2.orient == LEFT:  # R L
                    if not all([
                        read.reference_end >= self.break1.start,
                        mate.reference_start + 1 <= self.break2.end
                    ]):
                        continue
                else:  # R R
                    if not all([
                        read.reference_end >= self.break1.start,
                        mate.reference_end >= self.break2.start,
                        read.reference_end < self.break2.end
                    ]):
                        continue
            if is_compatible:
                self.compatible_flanking_pairs.add((read, mate))
            else:
                self.flanking_pairs.add((read, mate))

    def add_break1_split_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): putative split read supporting the first breakpoint
        """
        try:
            p = read_tools.breakpoint_pos(read, self.break1.orient) + 1
            if Interval.overlaps((p, p), self.break1):
                self.break1_split_reads.add(read)
        except AttributeError:
            pass

    def add_break2_split_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): putative split read supporting the second breakpoint
        """
        try:
            p = read_tools.breakpoint_pos(read, self.break2.orient) + 1
            if Interval.overlaps((p, p), self.break2):
                self.break2_split_reads.add(read)
        except AttributeError:
            pass

    def add_spanning_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): putative spanning read
        """
        bpp, event_types = _call_by_reads(self.source_evidence, read)
        if self.event_type in event_types:
            if bpp == self:
                self.spanning_reads.add(read)

    def __hash__(self):
        raise NotImplementedError('this object type does not support hashing')

    def flanking_metrics(self):
        """
        computes the median and standard deviation of the flanking pairs. Note that standard
        deviation is calculated wrt the median and not the average. Also that the fragment size
        is calculated as a range so the start and end of the range are used in computing these
        metrics

        Returns:
            tuple:
                - ``float`` - the median fragment size
                - ``float`` - the fragment size standard deviation wrt the median
        """
        fragment_sizes = []
        for read, mate in self.flanking_pairs:
            # check that the fragment size is reasonable
            f = self.source_evidence.compute_fragment_size(read, mate)
            fragment_sizes.append(f.start)
            fragment_sizes.append(f.end)
        median = 0
        stdev = 0
        if len(fragment_sizes) > 0:
            median = statistics.median(fragment_sizes)
            err = 0
            for insert in fragment_sizes:
                err += math.pow(insert - median, 2)
            err /= len(fragment_sizes)
            stdev = math.sqrt(err)
        return median, stdev

    def break1_tgt_align_split_read_names(self):
        reads = set()
        for r in self.break1_split_reads:
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                reads.add(r.query_name)
        return reads

    def break2_tgt_align_split_read_names(self):
        reads = set()
        for r in self.break2_split_reads:
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                reads.add(r.query_name)
        return reads

    def linking_split_read_names(self):
        reads1 = set()
        for r in self.break1_split_reads:
            reads1.add(r.query_name)

        reads2 = set()
        for r in self.break2_split_reads:
            reads2.add(r.query_name)

        return reads1 & reads2

    def flatten(self):
        row = self.source_evidence.flatten()
        row.update(BreakpointPair.flatten(self))  # this will overwrite the evidence breakpoint which is what we want
        row.update({
            COLUMNS.break1_call_method: self.call_method[0],
            COLUMNS.break2_call_method: self.call_method[1],
            COLUMNS.event_type: self.event_type
        })
        median, stdev = self.flanking_metrics()
        flank = set()
        for f, m in self.flanking_pairs:
            flank.update({f.query_name, m.query_name})
        row.update({
            COLUMNS.flanking_pairs: len(self.flanking_pairs),
            COLUMNS.flanking_median_fragment_size: median,
            COLUMNS.flanking_stdev_fragment_size: stdev,
            COLUMNS.flanking_pairs_read_names: ';'.join(sorted(list(flank)))
        })

        b1 = set()
        b1_tgt = set()
        b2 = set()
        b2_tgt = set()

        for r in self.break1_split_reads:
            name = r.query_name
            b1.add(name)
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                b1_tgt.add(name)
        for r in self.break2_split_reads:
            name = r.query_name
            b2.add(name)
            if r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                b2_tgt.add(name)

        linking = b1 & b2

        row.update({
            COLUMNS.break1_split_reads: len(b1),
            COLUMNS.break1_split_reads_forced: len(b1_tgt),
            COLUMNS.break1_split_read_names: ';'.join(sorted(b1)),
            COLUMNS.break2_split_reads: len(b2),
            COLUMNS.break2_split_reads_forced: len(b2_tgt),
            COLUMNS.break2_split_read_names: ';'.join(sorted(b2)),
            COLUMNS.linking_split_reads: len(linking),
            COLUMNS.linking_split_read_names: ';'.join(sorted(linking)),
            COLUMNS.spanning_reads: len(self.spanning_reads),
            COLUMNS.spanning_read_names: ';'.join(sorted([r.query_name for r in self.spanning_reads]))
        })

        if self.contig:
            r1, r2 = self.contig_alignment
            ascore = r1.get_tag('br')
            if r2:
                ascore = int(round((r1.get_tag('br') + r2.get_tag('br')) / 2, 0))
            cseq = self.contig_alignment[0].query_sequence
            qc1 = r1.query_coverage_interval()
            qc2 = qc1
            if r2:
                qc2 = r2.query_coverage_interval()
                if r2.is_reverse != r1.is_reverse:
                    qc2 = Interval(len(self.contig.seq) - qc2.end, len(self.contig.seq) - qc2.start)
            caqc = len(qc1 | qc2) if not Interval.overlaps(qc1, qc2) else len(qc1) + len(qc2)
            row.update({
                COLUMNS.contig_seq: cseq,  # don't output sequence directly from contig b/c must always be wrt to the positive strand
                COLUMNS.contig_remap_score: self.contig.remap_score(),
                COLUMNS.contig_alignment_score: ascore,
                COLUMNS.contig_remapped_reads: len(self.contig.input_reads),
                COLUMNS.contig_remapped_read_names:
                    ';'.join(sorted(set([r.query_name for r in self.contig.input_reads]))),
                COLUMNS.contig_strand_specific: self.contig.strand_specific,
                COLUMNS.contig_alignment_query_coverage: caqc,
                COLUMNS.contig_build_score: self.contig.score,
                COLUMNS.contig_alignment_query_name: self.contig_alignment[0].query_name
            })
        return row


def _call_by_reads(source_evidence, read1, read2=None):
    """
    for any read or given set of reads calls a breakpoint pair
    also ensures that the call is compatible with the source_evidence object
    putative event types
    """
    try:
        bpp = BreakpointPair.call_breakpoint_pair(read1, read2)
        if bpp.opposing_strands != source_evidence.opposing_strands:
            return None, []
        putative_event_types = set(source_evidence.putative_event_types())
        if set([SVTYPE.DUP, SVTYPE.INS]) & putative_event_types:
            putative_event_types.update([SVTYPE.DUP, SVTYPE.INS])

        if len(set(BreakpointPair.classify(bpp)) & putative_event_types) == 0:
            return None, []
        if source_evidence.stranded:  # strand specific
            if any([
                bpp.break1.strand != source_evidence.break1.strand,
                bpp.break2.strand != source_evidence.break2.strand
            ]):
                return None, []
        else:
            bpp.stranded = False
            bpp.break1.strand = STRAND.NS
            bpp.break2.strand = STRAND.NS

        calls = []
        for event_type in putative_event_types:
            if event_type == SVTYPE.INS:
                if len(bpp.untemplated_seq) == 0 or \
                        len(bpp.untemplated_seq) <= abs(Interval.dist(bpp.break1, bpp.break2)):
                    continue
            elif event_type == SVTYPE.DEL:
                if len(bpp.untemplated_seq) > abs(Interval.dist(bpp.break1, bpp.break2)):
                    continue
            if event_type not in BreakpointPair.classify(bpp):
                continue
            calls.append(event_type)
        return bpp, calls
    except UserWarning as err:
        return None, []


def _call_by_contigs(source_evidence):
    # try calling by contigs
    contig_calls = []
    for ctg in source_evidence.contigs:
        for read1, read2 in ctg.alignments:
            bpp, event_types = _call_by_reads(source_evidence, read1, read2)
            for event_type in event_types:
                new_event = EventCall(
                    bpp.break1,
                    bpp.break2,
                    source_evidence,
                    event_type,
                    contig=ctg,
                    contig_alignment=(read1, read2),
                    untemplated_seq=bpp.untemplated_seq,
                    call_method=CALL_METHOD.CONTIG
                )
                # add the flanking support
                new_event.add_flanking_support(source_evidence.flanking_pairs)
                if new_event.has_compatible:
                    new_event.add_flanking_support(source_evidence.compatible_flanking_pairs, is_compatible=True)

                # add any spanning reads that call the same event
                for read in source_evidence.spanning_reads:
                    new_event.add_spanning_read(read)
                # add any split read support (this will be consumed for non-contig calls)
                for read in source_evidence.split_reads[0]:
                    new_event.add_break1_split_read(read)
                for read in source_evidence.split_reads[1]:
                    new_event.add_break2_split_read(read)

                contig_calls.append(new_event)
    return contig_calls


def filter_consumed_pairs(pairs, consumed_reads):
    """
    given a set of read tuples, returns all tuples where neither read in the tuple is in the consumed set

    Args:
        pairs (set of tuples of :class:`pysam.AlignedSegment` and :class:`pysam.AlignedSegment`): pairs to be filtered
        consumed_reads: (set of :class:`pysam.AlignedSegment`): set of reads that have been used/consumed

    Returns:
        set of tuples of :class:`pysam.AlignedSegment` and :class:`pysam.AlignedSegment`: set of filtered tuples
    
    Note:
        this will work with any hash-able object

    Example:
        >>> pairs = {(1, 2), (3, 4), (5, 6)}
        >>> consumed_reads = {1, 2, 4}
        >>> filter_consumed_pairs(pairs, consumed_reads)
        {(5, 6)}
    """
    temp = set()
    for read, mate in pairs:
        if read not in consumed_reads and mate not in consumed_reads:
            temp.add((read, mate))
    return temp


def _call_by_spanning_reads(source_evidence, consumed_evidence):
    spanning_calls = {}
    available_flanking_pairs = filter_consumed_pairs(source_evidence.flanking_pairs, consumed_evidence)

    for read in source_evidence.spanning_reads - consumed_evidence:
        bpp, event_types = _call_by_reads(source_evidence, read)
        for event_type in event_types:
            spanning_calls.setdefault((bpp, event_type), set()).add(read)
    result = []
    for k, reads in spanning_calls.items():
        if len(reads) < source_evidence.min_spanning_reads_resolution:
            continue
        bpp, event_type = k
        bpp.break1.seq = None  # unless we are collecting a consensus we shouldn't assign sequences to the breaks
        bpp.break2.seq = None
        new_event = EventCall(
            bpp.break1, bpp.break2,
            source_evidence,
            event_type,
            CALL_METHOD.SPAN,
            break2_call_method=CALL_METHOD.SPAN,
            untemplated_seq=bpp.untemplated_seq
        )
        new_event.spanning_reads.update(reads)
        # add any supporting split reads
        # add the flanking support
        new_event.add_flanking_support(available_flanking_pairs)
        if new_event.has_compatible:
            new_event.add_flanking_support(available_flanking_pairs, is_compatible=True)
        # add any split read support (this will be consumed for non-contig calls)
        for read in source_evidence.split_reads[0] - consumed_evidence:
            new_event.add_break1_split_read(read)
        for read in source_evidence.split_reads[1] - consumed_evidence:
            new_event.add_break2_split_read(read)
        result.append(new_event)
    return result


def call_events(source_evidence):
    """
    generates a set of event calls based on the evidence associated with the source_evidence object
    will also narrow down the event type

    Args:
        source_evidence (Evidence): the input evidence
        event_type (SVTYPE): the type of event we are collecting evidence for

    Returns:
        :class:`list` of :class:`EventCall`: list of calls
    """
    consumed_evidence = set()  # keep track to minimize evidence re-use
    calls = []
    errors = set()

    contig_calls = _call_by_contigs(source_evidence)
    calls.extend(contig_calls)
    for call in contig_calls:
        consumed_evidence.update(call.support())

    spanning_calls = _call_by_spanning_reads(source_evidence, consumed_evidence)
    for call in spanning_calls:
        consumed_evidence.update(call.support())
    calls.extend(spanning_calls)

    for event_type in sorted(source_evidence.putative_event_types()):
        # try calling by split/flanking reads
        try:
            contig_consumed_evidence = set()
            contig_consumed_evidence.update(consumed_evidence)
            calls.extend(_call_by_supporting_reads(source_evidence, event_type, contig_consumed_evidence))
        except UserWarning as err:
            errors.add(str(err))

    if len(calls) == 0 and len(errors) > 0:
        raise UserWarning(';'.join(sorted(list(errors))))
    elif len(calls) == 0:
        raise UserWarning('insufficient evidence to call events')
    return calls


def distance(start, end):
    return Interval(abs(end - start) + 1)


def traverse(start, distance, orientation):
    if orientation == ORIENT.LEFT:
        return Interval(start - distance)
    elif orientation == ORIENT.RIGHT:
        return Interval(start + distance)
    else:
        raise ValueError('invalid value for orientation', orientation)


def _call_interval_by_flanking_coverage(
    coverage, orientation, max_expected_fragment_size, read_length,
    distance=distance, traverse=traverse
):
    if max_expected_fragment_size <= 0 or read_length <= 0:
        raise ValueError(
            'max_expected_fragment_size and read_length must be positive integers',
            max_expected_fragment_size, read_length)

    coverage_d = distance(coverage.start, coverage.end).start  # minimum distance of the coverage
    max_interval = max_expected_fragment_size - read_length
    assert(coverage_d <= max_interval)
    if orientation == ORIENT.LEFT:
        s = coverage.end
        t = traverse(coverage.end, max_interval - coverage_d, ORIENT.RIGHT).end
        return Interval(s, t)
    elif orientation == ORIENT.RIGHT:
        t = coverage.start
        s = max([1, traverse(coverage.start, max_interval - coverage_d, ORIENT.LEFT).start])
        return Interval(s, t)
    else:
        raise ValueError('orientation must be specific', orientation)


def _call_by_flanking_pairs(
        ev, event_type, first_breakpoint_called=None, second_breakpoint_called=None, consumed_evidence=None):
    """
    Given a set of flanking reads, computes the coverage interval (the area that is covered by flanking read alignments)
    this area gives the starting position for computing the breakpoint interval.

    .. todo::

        pre-split pairs into clusters by position and fragment size. This will enable calling multiple
        events in close proximity by flanking reads only. It will also aid in stopping FP reads from
        interfering with resolving events by flanking pairs.
    """
    if consumed_evidence is None:
        consumed_evidence = set()
    # for all flanking read pairs mark the farthest possible distance to the breakpoint
    # the start/end of the read on the breakpoint side
    first_positions = []
    second_positions = []

    flanking_count = 0
    cover1_reads = []
    cover2_reads = []
    available_flanking_pairs = filter_consumed_pairs(ev.flanking_pairs, consumed_evidence)

    for read, mate in available_flanking_pairs:
        # check that the fragment size is reasonable
        fragment_size = ev.compute_fragment_size(read, mate)
        if event_type == SVTYPE.DEL:
            if fragment_size.end <= ev.max_expected_fragment_size:
                continue
        elif event_type == SVTYPE.INS:
            if fragment_size.start >= ev.min_expected_fragment_size:
                continue
        flanking_count += 1
        cover1_reads.append(read)
        cover2_reads.append(mate)
        first_positions.extend([read.reference_start + 1, read.reference_end])
        second_positions.extend([mate.reference_start + 1, mate.reference_end])

    if flanking_count < ev.min_flanking_pairs_resolution:
        raise AssertionError('insufficient coverage to call by flanking reads')

    cover1 = Interval(min(first_positions), max(first_positions))
    cover2 = Interval(min(second_positions), max(second_positions))
    if not ev.interchromosomal:
        if Interval.overlaps(cover1, cover2) and event_type != SVTYPE.DUP:
            raise AssertionError('flanking read coverage overlaps. cannot call by flanking reads', cover1, cover2)
        elif event_type == SVTYPE.DUP and (cover1.start > cover2.start or cover2.end < cover1.end):
            raise AssertionError('flanking coverage for duplications must have some distinct positions', cover1, cover2)

    break1_strand = STRAND.NS
    break2_strand = STRAND.NS
    if ev.stranded:
        break1_strand = ev.decide_sequenced_strand(cover1_reads)
        break2_strand = ev.decide_sequenced_strand(cover2_reads)

    distance_func = distance
    traverse_func = traverse

    if first_breakpoint_called is None and second_breakpoint_called is None:
        # call the first breakpoint
        if ev.protocol == PROTOCOL.TRANS:
            def distance_func(s, t):
                return TranscriptomeEvidence.compute_exonic_distance(
                    s, t, ev.overlapping_transcripts[0])

            def traverse_func(s, d, o):
                return TranscriptomeEvidence.traverse_exonic_distance(
                    s, d, o, ev.overlapping_transcripts[0])

        window1 = _call_interval_by_flanking_coverage(
            cover1, ev.break1.orient, ev.max_expected_fragment_size, ev.read_length,
            distance=distance_func, traverse=traverse_func
        )

        # call the second breakpoint
        if ev.protocol == PROTOCOL.TRANS:
            def distance_func(s, t):
                return TranscriptomeEvidence.compute_exonic_distance(
                    s, t, ev.overlapping_transcripts[1])

            def traverse_func(s, d, o):
                return TranscriptomeEvidence.traverse_exonic_distance(
                    s, d, o, ev.overlapping_transcripts[1])

        window2 = _call_interval_by_flanking_coverage(
            cover2, ev.break2.orient, ev.max_expected_fragment_size, ev.read_length,
            distance=distance_func, traverse=traverse_func
        )

        if not ev.interchromosomal:
            if window1.start > window2.end:
                raise AssertionError('flanking window regions are incompatible', window1, window2)
            window1.end = min([window1.end, window2.end, cover2.start - (0 if event_type == SVTYPE.DUP else 1)])
            window2.start = max([window1.start, window2.start, cover1.end + (0 if event_type == SVTYPE.DUP else 1)])
        first_breakpoint_called = Breakpoint(
            ev.break1.chr, window1.start, window1.end,
            orient=ev.break1.orient,
            strand=break1_strand
        )

        second_breakpoint_called = Breakpoint(
            ev.break2.chr, window2.start, window2.end,
            orient=ev.break2.orient,
            strand=break2_strand
        )
        return first_breakpoint_called, second_breakpoint_called

    elif second_breakpoint_called is None:
        # does the input breakpoint make sense with the coverage?
        if any([
            first_breakpoint_called.orient == ORIENT.LEFT and cover1.end > first_breakpoint_called.end,
            first_breakpoint_called.orient == ORIENT.RIGHT and cover1.start < first_breakpoint_called.start
        ]):
            raise AssertionError(
                'input breakpoint is incompatible with flanking coverage', cover1, first_breakpoint_called)

        # call the second breakpoint
        if ev.protocol == PROTOCOL.TRANS:
            def distance_func(s, t):
                return TranscriptomeEvidence.compute_exonic_distance(
                    s, t, ev.overlapping_transcripts[1])

            def traverse_func(s, d, o):
                return TranscriptomeEvidence.traverse_exonic_distance(
                    s, d, o, ev.overlapping_transcripts[1])

        window = _call_interval_by_flanking_coverage(
            cover2, ev.break2.orient, ev.max_expected_fragment_size, ev.read_length,
            distance=distance_func, traverse=traverse_func
        )
        # trim the putative window by the input breakpoint location for intrachromosomal events
        if not ev.interchromosomal:
            window.start = max([
                window.start, first_breakpoint_called.start + (0 if event_type == SVTYPE.DUP else 1), cover1.end + 1])
            if window.start > window.end or window.end < first_breakpoint_called.start:
                raise AssertionError('input breakpoint incompatible with call', window, first_breakpoint_called)
        second_breakpoint_called = Breakpoint(
            ev.break2.chr, window.start, window.end,
            orient=ev.break2.orient,
            strand=break2_strand
        )
        return first_breakpoint_called, second_breakpoint_called

    elif first_breakpoint_called is None:
        # does the input breakpoint make sense with the coverage?
        if any([
            second_breakpoint_called.orient == ORIENT.LEFT and cover2.end > second_breakpoint_called.end,
            second_breakpoint_called.orient == ORIENT.RIGHT and cover2.start < second_breakpoint_called.start
        ]):
            raise AssertionError(
                'input breakpoint is incompatible with flanking coverage', cover2, second_breakpoint_called)

        # call the first breakpoint
        if ev.protocol == PROTOCOL.TRANS:
            def distance_func(s, t):
                return TranscriptomeEvidence.compute_exonic_distance(
                    s, t, ev.overlapping_transcripts[0])

            def traverse_func(s, d, o):
                return TranscriptomeEvidence.traverse_exonic_distance(
                    s, d, o, ev.overlapping_transcripts[0])

        window = _call_interval_by_flanking_coverage(
            cover1, ev.break1.orient, ev.max_expected_fragment_size, ev.read_length,
            distance=distance_func, traverse=traverse_func
        )

        # trim the putative window by the input breakpoint location for intrachromosomal events
        if not ev.interchromosomal:
            window.end = min([
                window.end, second_breakpoint_called.end - (0 if event_type == SVTYPE.DUP else 1), cover2.start - 1])
            if window.end < window.start or window.start > second_breakpoint_called.end:
                raise AssertionError('input breakpoint incompatible with call', window, second_breakpoint_called)
        first_breakpoint_called = Breakpoint(
            ev.break1.chr, window.start, window.end,
            orient=ev.break1.orient,
            strand=break1_strand
        )
        return first_breakpoint_called, second_breakpoint_called
    else:
        raise ValueError('cannot input both breakpoints')


def _call_by_supporting_reads(ev, event_type, consumed_evidence=None):
    """
    use split read evidence to resolve bp-level calls for breakpoint pairs (where possible)
    if a bp level call is not possible for one of the breakpoints then returns None
    if no breakpoints can be resolved returns the original event only with NO split read evidence
    also sets the SV type call if multiple are input
    """
    if consumed_evidence is None:
        consumed_evidence = set()
    pos1 = {}
    pos2 = {}

    available_flanking_pairs = filter_consumed_pairs(ev.flanking_pairs, consumed_evidence)

    for i, breakpoint, d in [(0, ev.break1, pos1), (1, ev.break2, pos2)]:
        for read in ev.split_reads[i] - consumed_evidence:
            try:
                pos = read_tools.breakpoint_pos(read, breakpoint.orient) + 1
                if pos not in d:
                    d[pos] = set()
                d[pos].add(read)
            except AttributeError:
                pass
        putative_positions = list(d.keys())
        for pos in putative_positions:
            if len(d[pos]) < ev.min_splits_reads_resolution:
                del d[pos]
            else:
                count = 0
                for r in d[pos]:
                    if not r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) or \
                            not r.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                        count += 1
                if count < ev.min_non_target_aligned_split_reads:
                    del d[pos]

    linked_pairings = []
    # now pair up the breakpoints with their putative partners
    for first, second in itertools.product(pos1, pos2):
        if ev.break1.chr == ev.break2.chr:
            if first >= second:
                continue
        links = 0
        read_names = set([r.query_name for r in pos1[first]])
        reads = set([(r.query_name, r.query_sequence) for r in pos1[first]])
        tgt_align = 0
        for read in pos2[second]:
            if read.query_name in read_names:
                links += 1
            if (read.query_name, read.query_sequence) in reads:
                tgt_align += 1
        if links < ev.min_linking_split_reads:
            continue
        deletion_size = second - first - 1
        if tgt_align >= ev.min_double_aligned_to_estimate_insertion_size:
            # we can estimate the fragment size
            max_insert = ev.read_length - 2 * ev.min_softclipping
            if event_type == SVTYPE.INS and max_insert < deletion_size:
                continue
            elif event_type == SVTYPE.DEL and deletion_size < max_insert:
                continue
        elif links >= ev.min_double_aligned_to_estimate_insertion_size:
            if deletion_size > ev.max_expected_fragment_size and event_type == SVTYPE.INS:
                continue

        first_breakpoint = Breakpoint(ev.break1.chr, first, strand=ev.break1.strand, orient=ev.break1.orient)
        second_breakpoint = Breakpoint(ev.break2.chr, second, strand=ev.break2.strand, orient=ev.break2.orient)
        call = EventCall(
            first_breakpoint, second_breakpoint, ev, event_type,
            call_method=CALL_METHOD.SPLIT
        )
        call.add_flanking_support(available_flanking_pairs)
        if call.has_compatible:
            call.add_flanking_support(available_flanking_pairs, is_compatible=True)
        call.break1_split_reads.update(pos1[first])
        call.break2_split_reads.update(pos2[second])
        linked_pairings.append(call)

    for call in linked_pairings:
        consumed_evidence.update(call.support())

    error_messages = set()
    available_flanking_pairs = filter_consumed_pairs(available_flanking_pairs, consumed_evidence)

    try:
        f, s = _call_by_flanking_pairs(ev, event_type, consumed_evidence=consumed_evidence)
        call = EventCall(
            f, s, ev, event_type,
            call_method=CALL_METHOD.FLANK
        )
        call.add_flanking_support(available_flanking_pairs)
        if call.has_compatible:
            call.add_flanking_support(available_flanking_pairs, is_compatible=True)
        linked_pairings.append(call)
    except (AssertionError, UserWarning) as err:
        error_messages.add(str(err))

    if len(linked_pairings) == 0:
        raise UserWarning(';'.join(list(error_messages)))
    return linked_pairings
