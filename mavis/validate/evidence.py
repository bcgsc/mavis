import itertools

from .base import Evidence
from ..annotate.variant import overlapping_transcripts
from ..breakpoint import Breakpoint
from ..constants import ORIENT, PROTOCOL, STRAND, SVTYPE
from ..interval import Interval


class GenomeEvidence(Evidence):

    def __init__(self, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)
        self.protocol = PROTOCOL.GENOME

        self.outer_window1 = self.generate_window(self.break1)
        self.outer_window2 = self.generate_window(self.break2)
        self.inner_window1 = Interval(
            max([self.break1.start - self.call_error - self.read_length + 1, 1]),
            self.break1.end + self.call_error + self.read_length - 1)
        self.inner_window2 = Interval(
            max([self.break2.start - self.call_error - self.read_length + 1, 1]),
            self.break2.end + self.call_error + self.read_length - 1)

        if SVTYPE.INS in self.putative_event_types():
            comb = len(self.break1 | self.break2)
            if comb > len(self.break1) and comb > len(self.break2):
                compt_break1 = Breakpoint(
                    self.break1.chr, self.break1.start, self.break1.end, orient=ORIENT.RIGHT, strand=self.break1.strand)
                compt_break2 = Breakpoint(
                    self.break2.chr, self.break2.start, self.break2.end, orient=ORIENT.LEFT, strand=self.break2.strand)

                self.compatible_window1 = self.generate_window(compt_break1)
                self.compatible_window2 = self.generate_window(compt_break2)
        elif SVTYPE.DUP in self.putative_event_types():
            compt_break1 = Breakpoint(
                self.break1.chr, self.break1.start, self.break1.end, orient=ORIENT.LEFT, strand=self.break1.strand)
            compt_break2 = Breakpoint(
                self.break2.chr, self.break2.start, self.break2.end, orient=ORIENT.RIGHT, strand=self.break2.strand)

            self.compatible_window1 = self.generate_window(compt_break1)
            self.compatible_window2 = self.generate_window(compt_break2)

    def generate_window(self, breakpoint):
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
        start = breakpoint.start - self.max_expected_fragment_size - self.call_error + 1
        end = breakpoint.end + self.max_expected_fragment_size + self.call_error - 1

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + self.call_error + self.read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - self.call_error - self.read_length + 1
        return Interval(max([1, start]), max([end, 1]))

    def compute_fragment_size(self, read, mate=None):
        return Interval(abs(read.template_length))


class TranscriptomeEvidence(Evidence):

    def __init__(self, annotations, *pos, **kwargs):
        Evidence.__init__(self, *pos, **kwargs)
        self.protocol = PROTOCOL.TRANS
        # get the list of overlapping transcripts
        self.overlapping_transcripts = overlapping_transcripts(annotations, self.break1) | overlapping_transcripts(annotations, self.break2)

        self.outer_window1 = self.generate_window(self.break1)
        self.outer_window2 = self.generate_window(self.break2)
        tgt = self.call_error + self.read_length - 1

        self.inner_window1 = self.traverse(self.break1.end, tgt, ORIENT.RIGHT) | self.traverse(self.break1.start, tgt, ORIENT.LEFT)
        self.inner_window2 = self.traverse(self.break2.start, tgt, ORIENT.LEFT) | self.traverse(self.break2.end, tgt, ORIENT.RIGHT)

        if SVTYPE.INS in self.putative_event_types():
            comb = len(self.break1 | self.break2)
            if comb > len(self.break1) and comb > len(self.break2):
                compt_break1 = Breakpoint(
                    self.break1.chr, self.break1.start, self.break1.end, orient=ORIENT.RIGHT, strand=self.break1.strand)
                compt_break2 = Breakpoint(
                    self.break2.chr, self.break2.start, self.break2.end, orient=ORIENT.LEFT, strand=self.break2.strand)
                self.compatible_window1 = self.generate_window(compt_break1)
                self.compatible_window2 = self.generate_window(compt_break2)

        elif SVTYPE.DUP in self.putative_event_types():
            compt_break1 = Breakpoint(
                self.break1.chr, self.break1.start, self.break1.end, orient=ORIENT.LEFT, strand=self.break1.strand)
            compt_break2 = Breakpoint(
                self.break2.chr, self.break2.start, self.break2.end, orient=ORIENT.RIGHT, strand=self.break2.strand)

            self.compatible_window1 = self.generate_window(compt_break1)
            self.compatible_window2 = self.generate_window(compt_break2)

    def traverse(self, start, distance, direction, strand=STRAND.NS, chrom=None):
        """
        given some genomic position and a distance. Uses the input transcripts to
        compute all possible genomic end positions at that distance if intronic
        positions are ignored

        Args:
            start (int): the genomic start position
            distance (int): the amount of exonic/intergenic units to traverse
            direction (ORIENT): the direction wrt to the positive/forward reference strand to traverse
            transcripts (:class:`list` of :class:`UsTranscript`): list of transcripts to use
        """
        transcripts = self._select_transcripts(chrom, strand)
        is_left = True if direction == ORIENT.LEFT else False
        genomic_end_positions = set()
        normal_end = GenomeEvidence.traverse(start, distance, direction).start

        for transcript in itertools.chain.from_iterable([ust.transcripts for ust in transcripts]):
            # convert the start to cdna coordinates
            if any([
                start < transcript.reference_object.start and is_left,
                start > transcript.reference_object.end and not is_left
            ]):
                continue
            cdna_start, start_shift = transcript.convert_genomic_to_nearest_cdna(
                start, stick_direction=ORIENT.LEFT if is_left else ORIENT.RIGHT, allow_outside=True)
            if abs(start_shift) > distance:  # entirely within an intron
                continue
            if transcript.is_reverse:
                if is_left:
                    cdna_end = cdna_start + (distance - start_shift)
                else:
                    cdna_end = cdna_start - (distance + start_shift)
            else:
                if is_left:
                    cdna_end = cdna_start - (distance - start_shift)
                else:
                    cdna_end = cdna_start + (distance + start_shift)
            if cdna_end <= 0:
                cdna_end -= 1
            # convert the cdna end back to genomic coordinates
            genomic_end = transcript.convert_cdna_to_genomic(cdna_end)
            if genomic_end == normal_end:
                continue
            genomic_end_positions.add(genomic_end)

        if not genomic_end_positions:
            genomic_end_positions.add(normal_end)
        return Interval.from_iterable(genomic_end_positions)

    def compute_fragment_size(self, read, mate):
        if read.reference_start > mate.reference_start:
            read, mate = mate, read
        if read.reference_name == mate.reference_name:
            start, end = self.distance(read.reference_start + 1, mate.reference_end, chrom=read.reference_name)
            return Interval(start + 1, end + 1)
        return Interval(0)

    def _select_transcripts(self, chrom=None, strand=STRAND.NS):
        result = []
        for transcript in self.overlapping_transcripts:
            if (chrom is None or transcript.get_chr() == chrom) and STRAND.compare(transcript.strand, strand):
                result.append(transcript)
        return result

    def distance(self, start, end, strand=STRAND.NS, chrom=None):
        """
        give the current list of transcripts, computes the putative exonic/intergenic distance
        given two genomic positions. Intronic positions are ignored

        Intergenic calculations are only done if exonic only fails
        """
        exonic = []
        mixed = []
        inter = []
        transcripts = self._select_transcripts(chrom, strand)
        # try to calculate assuming the positions are exonic
        for transcript in itertools.chain.from_iterable([t.transcripts for t in transcripts]):
            if not transcript.reference_object.position & Interval(start, end):
                continue
            cdna_start, start_shift = transcript.convert_genomic_to_nearest_cdna(start, stick_direction=ORIENT.RIGHT)
            cdna_end, end_shift = transcript.convert_genomic_to_nearest_cdna(end, stick_direction=ORIENT.LEFT)
            dist = abs(cdna_end - cdna_start) + abs(start_shift) + abs(end_shift)
            if start_shift != 0 and end_shift != 0:
                inter.append(dist)
            elif start_shift != 0 or end_shift != 0:
                mixed.append(dist)
            else:
                exonic.append(dist)
        if exonic:
            return Interval.from_iterable(exonic)
        elif mixed:
            return Interval.from_iterable(mixed)
        elif inter:
            return Interval.from_iterable(inter)
        return Evidence.distance(start, end)

    def generate_window(self, breakpoint):
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
        window = GenomeEvidence.generate_window(self, breakpoint)
        tgt_left = Evidence.distance(window.start, breakpoint.start)  # amount to expand to the left
        tgt_right = Evidence.distance(breakpoint.end, window.end)  # amount to expand to the right
        window1 = self.traverse(breakpoint.start, tgt_left.end, ORIENT.LEFT, strand=breakpoint.strand, chrom=breakpoint.chr)
        window2 = self.traverse(breakpoint.end, tgt_right.end, ORIENT.RIGHT, strand=breakpoint.strand, chrom=breakpoint.chr)
        return window1 | window2
