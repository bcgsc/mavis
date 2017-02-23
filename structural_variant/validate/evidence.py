from .base import Evidence
from ..interval import Interval
from ..constants import ORIENT, PROTOCOL, STRAND
from ..annotate.variant import overlapping_transcripts
import itertools


class GenomeEvidence(Evidence):

    def __init__(self, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)
        self.protocol = PROTOCOL.GENOME

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
        print('GenomeEvidence outer {0}:{2}-{3} {1}:{4}-{5} inner {0}:{6}-{7}  {1}:{8}-{9}'.format(
            self.break1.chr, self.break2.chr, self.outer_window1[0], self.outer_window1[1],
            self.outer_window2[0], self.outer_window2[1], self.inner_window1[0], self.inner_window1[1],
            self.inner_window2[0], self.inner_window2[1]))

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
        return Interval(max([1, start]), max([end, 1]))

    def compute_fragment_size(self, read, mate=None):
        return Interval(abs(read.template_length))


class TranscriptomeEvidence(Evidence):
    def __init__(self, ANNOTATIONS, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)
        self.protocol = PROTOCOL.TRANS
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
        tgt = self.call_error + self.read_length - 1
        temp = TranscriptomeEvidence.traverse_exonic_distance(
            self.break1.start, tgt, ORIENT.LEFT, self.overlapping_transcripts[0])
        w1 = TranscriptomeEvidence.traverse_exonic_distance(
            self.break1.end, tgt, ORIENT.RIGHT, self.overlapping_transcripts[0])
        w1 = w1 | temp

        temp = TranscriptomeEvidence.traverse_exonic_distance(
            self.break2.start, tgt, ORIENT.LEFT, self.overlapping_transcripts[1])
        w2 = TranscriptomeEvidence.traverse_exonic_distance(
            self.break2.end, tgt, ORIENT.RIGHT, self.overlapping_transcripts[1])
        w2 = w2 | temp

        self.inner_windows = (w1, w2)
    
    def traverse_exonic_distance(start, distance, direction, transcripts):
        is_left = True if direction == ORIENT.LEFT else False
        input_distance = distance
        positions = [start - distance + 1 if is_left else start + distance - 1]

        for ust in transcripts:
            distance = input_distance  # reset the distance between transcripts
            pos = start
            if is_left:
                if start > ust.end:
                    available = start - ust.end
                    if available <= distance:
                        distance -= available
                        pos = ust.end + 1
                    else:
                        pos = start - distance + 1
                        distance = 0
                elif start < ust.start:
                    pos = start - distance + 1
                    distance = 0

                for i, ex in enumerate(ust.exons[::-1]):
                    if distance == 0:
                        break
                    if start >= ex.start and start <= ex.end:  # within this exon
                        available = start - ex.start + 1
                        if available <= distance:
                            pos = ex.start
                            distance -= available
                        else:  # more distance available than we want to travel
                            pos = start - distance + 1
                            distance = 0
                    elif start < ex.start:  # before this exon
                        continue
                    else:  # after this exon
                        available = len(ex)
                        if available <= distance:
                            pos = ex.start
                            distance -= len(ex)
                        else:  # more distance available than we want to travel
                            pos = ex.end - distance + 1
                            distance = 0
                if distance > 0:  # didn't consume all the distance, keep going
                    pos = ust.start - distance
            else:
                if start > ust.end:
                    pos = start + distance - 1
                    distance = 0
                elif start < ust.start:
                    available = ust.start - start
                    if available <= distance:
                        distance -= available
                        pos = ust.start - 1
                    else:
                        pos = start + distance - 1
                        distance = 0
                for i, ex in enumerate(ust.exons):
                    if distance == 0:
                        break
                    if start >= ex.start and start <= ex.end:  # within this exon
                        available = ex.end - start + 1
                        if available <= distance:
                            pos = ex.end
                            distance -= available
                        else:  # more distance available than we want to travel
                            pos = start + distance - 1
                            distance = 0
                    elif start < ex.start:  # before this exon
                        available = len(ex)
                        if available <= distance:
                            pos = ex.end
                            distance -= len(ex)
                        else:  # more distance available than we want to travel
                            pos = ex.start + distance - 1
                            distance = 0
                    else:  # after this exon
                        continue
                if distance > 0:  # didn't consume all the distance, keep going right
                    pos = ust.end + distance
            positions.append(pos)
        return Interval(min(positions), max(positions))

    def compute_fragment_size(self, read, mate):
        all_fragments = []
        if read.reference_start > mate.reference_start:
            read, mate = mate, read
        transcripts = self.overlapping_transcripts[0] & self.overlapping_transcripts[1]

        for t in itertools.chain.from_iterable([ust.transcripts for ust in transcripts]):
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

    def compute_distance(self, start, end):
        """
        give the current list of transcripts, computes the putative cdna distance
        given two genomic positions

        warning:
            only different from the genomic distance if both points land in the same transcript
        """
        all_fragments = []
        transcripts = self.overlapping_transcripts[0] | self.overlapping_transcripts[1]

        for t in itertools.chain.from_iterable([ust.transcripts for ust in transcripts]):
            try:
                cs = t.convert_genomic_to_nearest_cdna(start)
                ct = t.convert_genomic_to_nearest_cdna(end)
                cs = cs[0] - cs[1] if cs[1] < 0 else cs[0]
                ct = ct[0] + ct[1] if ct[1] > 0 else ct[0]
                fragment_size = abs(ct - cs) + 1
                all_fragments.append(fragment_size)
            except IndexError:
                pass
        if len(all_fragments) == 0:
            return Interval(start, end)
        else:
            return Interval(min(all_fragments), max(all_fragments))

    def get_genomic_offset(self, genomic_pos, cdna_distance):
        """
        given some genomic position and the amount of cdna we wish to travel computes the
        new genomic position based on the transcripts associated with this evidence 
        object
        """
        all_pos = []
        transcripts = self.overlapping_transcripts[0] | self.overlapping_transcripts[1]

        for t in itertools.chain.from_iterable([ust.transcripts for ust in transcripts]):
            try:
                pos, off = t.convert_genomic_to_nearest_cdna(genomic_pos)
                gpos = genomic_pos + cdna_distance
                if t.get_strand() == STRAND.NEG:
                    if pos - cdna_distance + 1 >= 1:
                        gpos = t.convert_cdna_to_genomic(pos - cdna_distance + 1)
                    else:
                        gpos = t.convert_cdna_to_genomic(1)
                        gpos += cdna_distance - pos
                else:
                    if pos + cdna_distance - 1 > t.end:
                        gpos = t.convert_cdna_to_genomic(t.end)
                        gpos += cdna_distance - (t.end - pos)
                    else:
                        gpos = t.convert_cdna_to_genomic(pos + cdna_distance - 1)
                all_pos.append(gpos)
            except IndexError:
                pass
        if len(all_pos) == 0:
            return Interval(genomic_pos + cdna_distance - 1)
        else:
            return Interval(min(all_pos), max(all_pos))

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

        w1 = TranscriptomeEvidence.traverse_exonic_distance(breakpoint.start, tgt_left, ORIENT.LEFT, transcripts)
        w2 = TranscriptomeEvidence.traverse_exonic_distance(breakpoint.end, tgt_right, ORIENT.RIGHT, transcripts)
        return w1 | w2
