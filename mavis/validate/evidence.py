import itertools

from .base import Evidence
from ..align import SplitAlignment, call_read_events
from ..bam import cigar as _cigar
from ..annotate.variant import overlapping_transcripts
from ..breakpoint import Breakpoint
from ..constants import ORIENT, PROTOCOL, STRAND, SVTYPE, CIGAR
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

        # set the transcriptome specific overrides
        if self.trans_min_mapping_quality is not None:
            self.min_mapping_quality = self.trans_min_mapping_quality
        if self.trans_fetch_reads_limit is not None:
            self.fetch_reads_limit = self.trans_fetch_reads_limit

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
            transcripts (:class:`list` of :class:`PreTranscript`): list of transcripts to use
        """
        transcripts = self._select_transcripts(chrom, strand)
        is_left = True if direction == ORIENT.LEFT else False
        genomic_end_positions = set()
        normal_end = GenomeEvidence.traverse(start, distance, direction).start

        for transcript in itertools.chain.from_iterable([pre_transcript.transcripts for pre_transcript in transcripts]):
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
        genomic_distance = Evidence.distance(start, end).end
        # try to calculate assuming the positions are exonic
        for transcript in itertools.chain.from_iterable([t.transcripts for t in transcripts]):
            if not transcript.reference_object.position & Interval(start, end):
                continue
            cdna_start, start_shift = transcript.convert_genomic_to_nearest_cdna(start)
            cdna_end, end_shift = transcript.convert_genomic_to_nearest_cdna(end)
            dist = abs(cdna_end - cdna_start) + abs(start_shift) + abs(end_shift)
            if cdna_start == cdna_end:
                dist = abs(start_shift - end_shift)
            if start_shift and end_shift:
                inter.append(dist)
            elif start_shift or end_shift:
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

    def min_cds_shift(self, pos, strand=STRAND.NS, chrom=None):
        exon_boundaries = set()
        for transcript in self._select_transcripts(chrom, strand):
            for exon in transcript.exons:
                exon_boundaries.update({exon.start, exon.end})
        return min(exon_boundaries, key=lambda x: abs(x - pos))

    def exon_boundary_shift_cigar(self, read):
        """
        given an input read, converts deletions to N when the deletion matches the exon boundaries. Also shifts alignments
        to correspond to the exon boundaries where possible
        """
        reference_pos = read.reference_start
        query_pos = 0
        new_cigar = []

        # collapsed transcript model
        exon_ends = set()
        exon_starts = set()
        for transcript in self._select_transcripts(read.reference_name):
            for exon in transcript.exons:
                exon_starts.add(exon.start)
                exon_ends.add(exon.end)
        refseq = self.reference_genome[read.reference_name].seq
        for i, (state, freq) in enumerate(read.cigar):
            # shift to coincide with exon boundaries if possible
            if new_cigar and i < len(read.cigar) - 1 and exon_ends and exon_starts:
                next_state, next_freq = read.cigar[i + 1]
                prev_state, prev_freq = new_cigar[-1]
                # compare deletions surrounded by exact alignments. Indels at exon boundaries will
                # be aligned the same as genome indels
                if state in {CIGAR.D, CIGAR.N} and {next_state} & {prev_state} & {CIGAR.EQ}:
                    nearest_end_boundary = min(exon_ends, key=lambda x: abs(x - reference_pos - 1))
                    prev_alignment_seq = refseq[reference_pos - prev_freq:reference_pos]
                    next_reference_pos = reference_pos + freq
                    next_alignment_seq = refseq[max(reference_pos, next_reference_pos - prev_freq):next_reference_pos]
                    shift = 0
                    for prev_base, next_base in zip(prev_alignment_seq[::-1], next_alignment_seq[::-1]):
                        if prev_base == next_base:
                            shift += 1
                        else:
                            break
                    shift = min(shift, reference_pos - nearest_end_boundary)
                    if shift > 0:
                        if shift == prev_freq:
                            del new_cigar[-1]
                        else:
                            new_cigar[-1] = (prev_state, prev_freq - shift)
                        new_cigar.extend([(state, freq), (CIGAR.EQ, shift)])
                        reference_pos += freq
                        query_pos += shift
                        continue
            if state in _cigar.REFERENCE_ALIGNED_STATES:
                reference_pos += freq
            if state in _cigar.QUERY_ALIGNED_STATES:
                query_pos += freq
            new_cigar.append((state, freq))
        # mark intron deletions as N instead of D
        reference_pos = read.reference_start
        for i, (state, freq) in enumerate(new_cigar):
            if state == CIGAR.D:
                dist = self.distance(reference_pos, reference_pos + freq + 1)
                if dist.start == 1:
                    state = CIGAR.N
            if state in _cigar.REFERENCE_ALIGNED_STATES:
                reference_pos += freq
            new_cigar[i] = (state, freq)
        return _cigar.join(new_cigar)

    def standardize_read(self, read):
        read.cigar = self.exon_boundary_shift_cigar(read)  # do this twice to avoid merging events accidentally
        read = Evidence.standardize_read(self, read)
        read.cigar = self.exon_boundary_shift_cigar(read)
        return read
