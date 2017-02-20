from .base import Evidence
from ..interval import Interval
from ..constants import ORIENT
from ..annotate.variant import overlapping_transcripts
import itertools


class GenomeEvidence(Evidence):

    def __init__(self, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)

        self.outer_windows = (
            GenomeEvidence._generate_window(
                self.break1,
                max_expected_fragment_size=self.max_expected_fragment_size,
                call_error=self.call_error,
                read_length=self.read_length
            ),
            GenomeEvidence._generate_window(
                self.break2,
                max_expected_fragment_size=self.max_expected_fragment_size,
                call_error=self.call_error,
                read_length=self.read_length
            )
        )
        self.inner_windows = (
            Interval(
                max([self.break1.start - self.call_error - self.read_length + 1, 1]),
                self.break1.end + self.call_error + self.read_length - 1
            ),
            Interval(
                max([self.break2.start - self.call_error - self.read_length + 1, 1]),
                self.break2.end + self.call_error + self.read_length - 1
            )
        )

    @staticmethod
    def _generate_window(breakpoint, max_expected_fragment_size, call_error, read_length):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            read_length (int): the read length
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        start = breakpoint.start - max_expected_fragment_size - call_error + 1
        end = breakpoint.end + max_expected_fragment_size + call_error - 1

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - read_length + 1
        return Interval(max([1, start]), end)

    def compute_fragment_size(self, read, mate=None):
        return Interval(abs(read.template_length))


class TranscriptomeEvidence(Evidence):
    def __init__(self, ANNOTATIONS, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)

        # get the list of overlapping transcripts
        self.overlapping_transcripts = (
            overlapping_transcripts(ANNOTATIONS, self.break1),
            overlapping_transcripts(ANNOTATIONS, self.break2)
        )

        self.outer_windows = (
            TranscriptomeEvidence._generate_window(
                self.break1,
                transcripts=self.overlapping_transcripts[0],
                read_length=self.read_length,
                call_error=self.call_error,
                max_expected_fragment_size=self.max_expected_fragment_size
            ),
            TranscriptomeEvidence._generate_window(
                self.break2,
                transcripts=self.overlapping_transcripts[1],
                read_length=self.read_length,
                call_error=self.call_error,
                max_expected_fragment_size=self.max_expected_fragment_size
            )
        )
        w1 = TranscriptomeEvidence._expand_breakpoint_interval_by_exonic(
            self.break1, self.overlapping_transcripts[0], self.call_error + self.read_length - 1
        )
        w2 = TranscriptomeEvidence._expand_breakpoint_interval_by_exonic(
            self.break2, self.overlapping_transcripts[1], self.call_error + self.read_length - 1
        )

        self.inner_windows = (w1, w2)

    def compute_fragment_size(self, read, mate):
        all_fragments = []
        if read.reference_start > mate.reference_start:
            read, mate = mate, read
        for t in itertools.chain.from_iterable([ust.transcripts for ust in self.transcripts]):
            try:
                cs = t.convert_genomic_to_nearest_cdna(read.reference_start + 1)
                ct = t.convert_genomic_to_nearest_cdna(mate.reference_end)
                cs = cs[0] - cs[1] if cs[1] < 0 else cs[0]
                ct = ct[0] + ct[1] if ct[1] > 0 else ct[0]
                fragment_size = abs(ct - cs) + 1
                all_fragments.append(fragment_size)
            except IndexError:
                pass
        if len(all_fragments) == 0:
            return GenomeEvidence.compute_fragment_size(self, read, mate)
        else:
            return Interval(min(all_fragments), max(all_fragments))

    @staticmethod
    def _generate_window(breakpoint, transcripts, read_length, call_error, max_expected_fragment_size):
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            annotations (dict of str and list of Gene): the set of reference annotations: genes, transcripts, etc
            read_length (int): the read length
            median_fragment_size (int): the median insert size
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
            stdev_fragment_size:
                the standard deviation away from the median for regular (non STV) read pairs
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        window = GenomeEvidence._generate_window(
            breakpoint,
            max_expected_fragment_size=max_expected_fragment_size,
            call_error=call_error,
            read_length=read_length
        )

        tgt_left = breakpoint.start - window.start + 1  # amount to expand to the left
        tgt_right = window.end - breakpoint.end + 1  # amount to expand to the right
        if len(transcripts) == 0:  # case 1. no overlapping transcripts
            return window
        return TranscriptomeEvidence._expand_breakpoint_interval_by_exonic(
            breakpoint, transcripts, tgt_left, tgt_right)
    
    @staticmethod
    def _expand_breakpoint_interval_by_exonic(breakpoint, transcripts, tgt_left, tgt_right):
        intervals = [(breakpoint.start - tgt_left + 1, breakpoint.end + tgt_right - 1)]

        for ust in transcripts:
            if breakpoint not in ust.position:
                continue
            mapping = {}

            cdna_length = sum([len(e) for e in ust.exons])
            s = 1
            for ex in ust.exons:
                mapping[Interval(ex.start, ex.end)] = Interval(s, s + len(ex) - 1)
                s += len(ex)
            reverse_mapping = {}
            for k, v in mapping.items():
                reverse_mapping[v] = k
            bs = breakpoint.start
            be = breakpoint.end
            # follow the window either direction to the nearest exon boundary
            for ex1, ex2 in zip(ust.exons, ust.exons[1::]):
                if bs > ex1.end and bs < ex2.start:
                    bs = ex1.end
                if be > ex1.end and be < ex2.start:
                    be = ex2.start
            assert(bs <= be)
            # now convert the exonic coordinates to determine the window boundaries
            cs = Interval.convert_pos(mapping, bs)
            gs = ust.start - (tgt_left - cs)
            if cs > tgt_left:
                gs = Interval.convert_pos(reverse_mapping, cs - tgt_left + 1)

            ce = Interval.convert_pos(mapping, be)
            ge = ust.end + (tgt_right - cdna_length + ce)
            if cdna_length > ce + tgt_right:
                ge = Interval.convert_pos(reverse_mapping, ce + tgt_right - 1)

            intervals.append(Interval(gs, ge))
        return Interval.union(*intervals)
