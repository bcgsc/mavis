from ..breakpoint import BreakpointPair, Breakpoint
from ..constants import CALL_METHOD, SVTYPE, PYSAM_READ_FLAGS, ORIENT, PROTOCOL, COLUMNS
from ..bam import read as read_tools
from ..interval import Interval
from ..error import NotSpecifiedError
import itertools
import statistics
import math
from copy import copy as sys_copy


class EventCall(BreakpointPair):
    """
    class for holding evidence and the related calls since we can't freeze the evidence object
    directly without a lot of copying. Instead we use call objects which are basically
    just a reference to the evidence object and decisions on class, exact breakpoints, etc
    """
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
            untemplated_seq=untemplated_seq,
            data=source_evidence.data
        )
        self.source_evidence = source_evidence
        self.flanking_pairs = set()
        self.break1_split_reads = set()
        self.break2_split_reads = set()
        self.compatible_flanking_pairs = set()
        # check that the event type is compatible
        self.event_type = SVTYPE.enforce(event_type)
        if event_type not in BreakpointPair.classify(self):
            raise ValueError(
                'event_type is not compatible with the breakpoint call', event_type, BreakpointPair.classify(self))
        self.contig = contig
        self.call_method = (CALL_METHOD.enforce(call_method), CALL_METHOD.enforce(break2_call_method))
        if contig and self.call_method != (CALL_METHOD.CONTIG, CALL_METHOD.CONTIG):
            raise ValueError('if a contig is given the call method must be by contig')
        self.contig_alignment = contig_alignment

    def supporting_reads(self):
        support = set()
        for read, mate in self.flanking_pairs:
            support.add(read)
            support.add(mate)
        support.update(self.break1_split_reads)
        support.update(self.break2_split_reads)
        if self.contig:
            support.update(self.contig.input_reads)
        return support

    def pull_flanking_support(self, flanking_pairs):
        """
        counts the flanking read-pair support for the event called. The original source evidence may
        have contained evidence for multiple events and uses a larger range so flanking pairs here
        are checked specifically against the current breakpoint call

        Returns:
            tuple:
            * :class:`set` of :class:`str` - set of the read query_names
            * :class:`int` - the median insert size
            * :class:`int` - the standard deviation (from the median) of the insert size
        """
        support = set()
        print('pull_flanking_support')
        
        fragment_sizes = []

        min_frag = max([
            self.source_evidence.min_expected_fragment_size + Interval.dist(self.break1, self.break2),
            self.source_evidence.max_expected_fragment_size])
        max_frag = len(self.break1 | self.break2) + self.source_evidence.max_expected_fragment_size

        encompass = len(self.break1 | self.break2)
        if self.event_type == SVTYPE.DEL and encompass < self.source_evidence.stdev_fragment_size:
            print('deletion', len(self.break1 | self.break2), self.untemplated_seq)
        elif self.event_type == SVTYPE.INS:
            print('insertion', len(self.break1 | self.break2), self.untemplated_seq)

        for read, mate in flanking_pairs:
            print(read.reference_start, read.reference_end, read.cigar, read.is_reverse)
            print(mate.reference_start, mate.reference_end, mate.cigar, mate.is_reverse)
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
            if not read_tools.orientation_supports_type(read, self.event_type):
                continue
            # check that the positions make sense
            if self.break1.orient == ORIENT.LEFT:
                if self.break2.orient == ORIENT.LEFT:  # L L
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
                if self.break2.orient == ORIENT.LEFT:  # R L
                    if not all([
                        read.reference_end >= self.break1.start,
                        mate.reference_start + 1 <= self.break2.end
                    ]):
                        print('position failed')
                        continue
                else:  # R R
                    if not all([
                        read.reference_end >= self.break1.start,
                        mate.reference_end >= self.break2.start,
                        read.reference_end < self.break2.end
                    ]):
                        continue
            self.flanking_pairs.add((read, mate))

    def __hash__(self):
        raise NotImplementedError('this object type does not support hashing')

    def __eq__(self, other):
        object.__eq__(self, other)

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
            COLUMNS.linking_split_read_names: ';'.join(sorted(linking))
        })

        if self.contig:
            r1, r2 = self.contig_alignment
            ascore = r1.get_tag('br')
            if r2:
                ascore = int(round((r1.get_tag('br') + r2.get_tag('br')) / 2, 0))
            row.update({
                COLUMNS.contig_seq: self.contig.seq,
                COLUMNS.contig_remap_score: self.contig.remap_score(),
                COLUMNS.contig_alignment_score: ascore,
                COLUMNS.contig_remapped_reads: len(self.contig.input_reads),
                COLUMNS.contig_remapped_read_names:
                    ';'.join(sorted(set([r.query_name for r in self.contig.input_reads])))
            })
        return row


def call_events(source_evidence):
    """
    generates a set of event calls based on the evidence associated with the source_evidence object
    will also narrow down the event type?

    Args:
        source_evidence (Evidence): the input evidence
        event_type (SVTYPE): the type of event we are collecting evidence for

    Returns:
        :class:`list` of :class:`EventCall`: list of calls

    .. figure:: ../_static/call_breakpoint_by_flanking_reads.svg

        model used in calculating the uncertainty interval for breakpoints called by flanking read pair evidence
    """
    consumed_evidence = set()  # keep track to minimize evidence re-use
    calls = []
    errors = set()
    contig_calls = []

    # try calling by contigs
    for ctg in source_evidence.contigs:
        for read1, read2 in ctg.alignments:
            curr = []
            try:
                bpp = BreakpointPair.call_breakpoint_pair(read1, read2)
            except UserWarning as err:
                continue
            if bpp.opposing_strands != source_evidence.opposing_strands:
                continue
            putative_event_types = set(source_evidence.putative_event_types())
            if set([SVTYPE.DUP, SVTYPE.INS]) & putative_event_types:
                putative_event_types.update([SVTYPE.DUP, SVTYPE.INS])

            if len(set(BreakpointPair.classify(bpp)) & putative_event_types) == 0:
                continue

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
                new_event.pull_flanking_support(source_evidence.flanking_pairs)
                # add any split read support (this will be consumed for non-contig calls)
                for read in source_evidence.split_reads[0]:
                    try:
                        p = read_tools.breakpoint_pos(read, source_evidence.break1.orient) + 1
                        if Interval.overlaps((p, p), source_evidence.break1):
                            new_event.break1_split_reads.add(read)
                    except AttributeError:
                        pass
                for read in source_evidence.split_reads[1]:
                    try:
                        p = read_tools.breakpoint_pos(read, source_evidence.break2.orient) + 1
                        if Interval.overlaps((p, p), source_evidence.break2):
                            new_event.break2_split_reads.add(read)
                    except AttributeError:
                        pass

                contig_calls.append(new_event)
    calls.extend(contig_calls)

    for call in contig_calls:
        consumed_evidence.update(call.contig.input_reads)
        consumed_evidence.update(call.break1_split_reads)
        consumed_evidence.update(call.break2_split_reads)
        consumed_evidence.update(call.flanking_pairs)

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


def _call_by_flanking_pairs(
        ev, event_type, first_breakpoint_called=None, second_breakpoint_called=None, consumed_evidence=None):
    """
    Given a set of flanking reads, computes the coverage interval (the area that is covered by flanking read alignments)
    this area gives the starting position for computing the breakpoint interval.

    .. figure:: _static/call_breakpoint_by_flanking_reads.svg

        model used in calculating the uncertainty interval for breakpoints called by flanking read pair evidence

    .. todo::

        pre-split pairs into clusters by position and fragment size. This will enable calling mutliple
        events in close proximity by flanking reads only. It will also aid in stopping FP reads from
        interfering with resolving events by flanking pairs.
    """
    if consumed_evidence is None:
        consumed_evidence = set()
    # for all flanking read pairs mark the farthest possible distance to the breakpoint
    # the start/end of the read on the breakpoint side
    first_positions = []
    second_positions = []
    if first_breakpoint_called and second_breakpoint_called:
        raise ValueError('do not bother calling when both breakpoints have already been called')

    flanking_count = 0
    for read, mate in ev.flanking_pairs:
        if (read, mate) in consumed_evidence:
            continue
        # check that the fragment size is reasonable
        fragment_size = ev.compute_fragment_size(read, mate)
        if event_type == SVTYPE.DEL:
            if fragment_size.end <= ev.max_expected_fragment_size:
                continue
        elif event_type == SVTYPE.INS:
            if fragment_size.start >= ev.min_expected_fragment_size:
                continue
        flanking_count += 1
        first_positions.extend([read.reference_start + 1, read.reference_end, mate.next_reference_start + 1])
        second_positions.extend([mate.reference_start + 1, mate.reference_end, read.next_reference_start + 1])

    if flanking_count < ev.min_flanking_pairs_resolution:
        raise AssertionError('insufficient coverage to call by flanking reads')

    cover1 = Interval(min(first_positions), max(first_positions))
    cover2 = Interval(min(second_positions), max(second_positions))

    if not ev.interchromosomal and Interval.overlaps(cover1, cover2) and event_type != SVTYPE.DUP:
        raise AssertionError('flanking read coverage overlaps. cannot call by flanking reads', cover1, cover2)
    if len(cover1) + ev.read_length * 2 > ev.max_expected_fragment_size or \
            len(cover2) + ev.read_length * 2 > ev.max_expected_fragment_size:
        raise AssertionError(
            'Cannot resolve by flanking reads. Coverage interval of flanking reads is larger than '
            'expected for normal variation. It is likely there are flanking reads for multiple events',
            cover1, cover2, ev.max_expected_fragment_size
        )

    cover1_length = len(cover1)
    cover2_length = len(cover2)

    if ev.protocol == PROTOCOL.TRANS:
        cover1_length = ev.compute_fragment_size

    if first_breakpoint_called is None:
        max_breakpoint_width = ev.max_expected_fragment_size - len(cover1) - ev.read_length * 2

        if ev.break1.orient == ORIENT.LEFT:
            end = cover1.end + max_breakpoint_width
            if not ev.interchromosomal:
                end = min([end, cover2.start - 1])
                if second_breakpoint_called:
                    end = min([end, second_breakpoint_called.end - 1])
            try:
                first_breakpoint_called = Breakpoint(
                    ev.break1.chr,
                    cover1.end, end,
                    orient=ev.break1.orient,
                    strand=ev.break1.strand
                )
            except AttributeError:
                raise AssertionError(
                    'input breakpoint is incompatible with flanking coverage region', cover1, second_breakpoint_called)
        elif ev.break1.orient == ORIENT.RIGHT:
            first_breakpoint_called = Breakpoint(
                ev.break1.chr,
                max([cover1.start - max_breakpoint_width, 1]),
                max([cover1.start, 1]),
                orient=ev.break1.orient,
                strand=ev.break1.strand
            )
        else:
            raise NotSpecifiedError('Cannot call by flanking if orientation was not given')

    if second_breakpoint_called is None:
        max_breakpoint_width = ev.max_expected_fragment_size - len(cover2) - ev.read_length * 2

        if ev.break2.orient == ORIENT.LEFT:
            second_breakpoint_called = Breakpoint(
                ev.break2.chr,
                cover2.end,
                cover2.end + max_breakpoint_width,
                orient=ev.break2.orient,
                strand=ev.break2.strand
            )
        elif ev.break2.orient == ORIENT.RIGHT:
            start = max([cover2.start - max_breakpoint_width, 1])
            if not ev.interchromosomal:
                start = max([start, cover1.end + 1])
                if first_breakpoint_called:
                    start = max([start, first_breakpoint_called.start + 1])
            try:
                second_breakpoint_called = Breakpoint(
                    ev.break2.chr,
                    start,
                    cover2.start,
                    orient=ev.break2.orient,
                    strand=ev.break2.strand
                )
            except AttributeError:
                raise AssertionError(
                    'input breakpoint is incompatible with flanking coverage region', cover2, first_breakpoint_called)
        else:
            raise NotSpecifiedError('Cannot call by flanking if orientation was not given')
    return first_breakpoint_called, second_breakpoint_called


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

    available_flanking_pairs = set()
    for pair in ev.flanking_pairs:
        if pair in consumed_evidence:
            continue
        available_flanking_pairs.add(pair)

    for i, breakpoint, d in [(0, ev.break1, pos1), (1, ev.break2, pos2)]:
        for read in ev.split_reads[i]:
            if read not in consumed_evidence:
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
        call.pull_flanking_support(available_flanking_pairs)
        call.break1_split_reads.update(pos1[first])
        call.break2_split_reads.update(pos2[second])
        linked_pairings.append(call)

    for call in linked_pairings:
        if call.break1.start in pos1:
            del pos1[call.break1.start]
        if call.break2.start in pos2:
            del pos2[call.break2.start]

    for first, second in itertools.product(pos1, pos2):
        if ev.break1.chr == ev.break2.chr:
            if first >= second:  # illegal combination, first breakpoint has to be before the second if intrachromosomal
                continue
        first_breakpoint = Breakpoint(ev.break1.chr, first, strand=ev.break1.strand, orient=ev.break1.orient)
        second_breakpoint = Breakpoint(ev.break2.chr, second, strand=ev.break2.strand, orient=ev.break2.orient)
        call = EventCall(
            first_breakpoint, second_breakpoint, ev, event_type,
            call_method=CALL_METHOD.SPLIT
        )
        call.pull_flanking_support(available_flanking_pairs)
        call.break1_split_reads.update(pos1[first])
        call.break2_split_reads.update(pos2[second])
        linked_pairings.append(call)

    for call in linked_pairings:
        consumed_evidence.update(call.flanking_pairs)

    available_flanking_pairs = available_flanking_pairs - consumed_evidence

    error_messages = set()
    # if can call the first breakpoint by split
    for pos in pos1:
        bp = sys_copy(ev.break1)
        bp.start = pos
        bp.end = pos
        try:
            f, s = _call_by_flanking_pairs(
                ev, event_type, first_breakpoint_called=bp, consumed_evidence=consumed_evidence)
            call = EventCall(
                f, s, ev, event_type,
                call_method=CALL_METHOD.SPLIT,
                break2_call_method=CALL_METHOD.FLANK
            )
            call.break1_split_reads.update(pos1[pos])
            call.pull_flanking_support(available_flanking_pairs)
            linked_pairings.append(call)

        except (AssertionError, UserWarning) as err:
            error_messages.add(str(err))

    for pos in pos2:
        bp = sys_copy(ev.break2)
        bp.start = pos
        bp.end = pos
        try:
            f, s = _call_by_flanking_pairs(
                ev, event_type, second_breakpoint_called=bp, consumed_evidence=consumed_evidence)
            call = EventCall(
                f, s, ev, event_type,
                call_method=CALL_METHOD.FLANK,
                break2_call_method=CALL_METHOD.SPLIT
            )
            call.break2_split_reads.update(pos2[pos])
            call.pull_flanking_support(available_flanking_pairs)
            linked_pairings.append(call)
        except (AssertionError, UserWarning) as err:
            error_messages.add(str(err))

    if len(linked_pairings) == 0:  # call by flanking only
        try:
            f, s = _call_by_flanking_pairs(ev, event_type, consumed_evidence=consumed_evidence)
            call = EventCall(
                f, s, ev, event_type,
                call_method=CALL_METHOD.FLANK
            )
            call.pull_flanking_support(available_flanking_pairs)
            linked_pairings.append(call)
        except (AssertionError, UserWarning) as err:
            error_messages.add(str(err))

    if len(linked_pairings) == 0:
        raise UserWarning(';'.join(list(error_messages)))
    return linked_pairings
