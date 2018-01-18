"""
Should take in a sam file from a aligner like bwa aln or bwa mem and convert it into a
"""
from copy import copy
import itertools
import os
import subprocess
import warnings

import pysam

from .bam import cigar as _cigar
from .bam import read as _read
from .breakpoint import BreakpointPair, Breakpoint
from .constants import CIGAR, COLUMNS, MavisNamespace, ORIENT, reverse_complement, STRAND, SVTYPE, NA_MAPPING_QUALITY
from .error import InvalidRearrangement
from .interval import Interval
from .util import devnull


SUPPORTED_ALIGNER = MavisNamespace(BWA_MEM='bwa mem', BLAT='blat', __name__='~mavis.align.SUPPORTED_ALIGNER')
""":class:`~mavis.constants.MavisNamespace`: supported aligners

- :term:`blat`
- :term:`bwa mem<BWA>`
"""


class SplitAlignment(BreakpointPair):

    def __init__(self, *pos, **kwargs):
        self.read1 = kwargs.pop('read1')
        self.read2 = kwargs.pop('read2', None)
        self.query_sequence = self.read1.query_sequence
        self.query_name = self.read1.query_name
        BreakpointPair.__init__(self, *pos, **kwargs)

    def query_coverage_read1(self):
        return query_coverage_interval(self.read1)

    def query_coverage_read2(self):
        seqlen = len(self.read1.query_sequence)
        qc1 = self.query_coverage_read1()
        qc2 = qc1
        if self.read2 is not None:
            qc2 = query_coverage_interval(self.read2)
            if self.read2.is_reverse != self.read1.is_reverse:
                qc2 = Interval(seqlen - qc2.end, seqlen - qc2.start)
        return qc2

    def query_coverage(self):
        """
        interval representing the total region of the input sequence that is covered by the combination of alignments
        """
        if self.read2 is None:
            return self.query_coverage_read1()
        return self.query_coverage_read1() | self.query_coverage_read2()

    def query_consumption(self):
        """
        fraction of the query sequence which is aligned (everything not soft-clipped) in either alignment
        """
        if self.read2 is None or Interval.overlaps(self.query_coverage_read1(), self.query_coverage_read2()):
            return len(self.query_coverage()) / self.read1.query_length
        return (len(self.query_coverage_read1()) + len(self.query_coverage_read2())) / self.read1.query_length

    def query_overlap_extension(self):
        if self.read2 is not None:
            max_init_overlap = max(len(self.query_coverage_read1()), len(self.query_coverage_read2()))
            total_overlap = len(self.query_coverage()) - max_init_overlap
            return total_overlap
        return 0

    def score(self, consec_bonus=10):
        """
        scores events between 0 and 1 penalizing events interrupting the alignment. Counts a split
        alignment as a single event
        """
        def score_matches(cigar):
            return sum([v + (v - 1) * consec_bonus for s, v in cigar if s == CIGAR.EQ])
        score = score_matches(_cigar.join([(s, f) for s, f in self.read1.cigar if s != CIGAR.N]))

        max_score = sum([v for s, v in self.read1.cigar if s in _cigar.ALIGNED_STATES])
        if self.read2:
            max_score += sum([v for s, v in self.read2.cigar if s in _cigar.ALIGNED_STATES])
        max_score = score_matches([(CIGAR.EQ, max_score)])
        if self.read2:
            score += score_matches([(s, f) for s, f in self.read2.cigar if s != CIGAR.N])
        return score / max_score

    def mapping_quality(self):
        if not self.read2:
            return Interval(self.read1.mapping_quality)
        return Interval(self.read1.mapping_quality) | Interval(self.read2.mapping_quality)

    @staticmethod
    def breakpoint_contig_remapped_depth(breakpoint, contig, read):
        if breakpoint.chr != read.reference_name:
            raise AssertionError('breakpoint chromosome does not match read reference', breakpoint, read.reference_name)
        if len(breakpoint) > 1:
            raise NotImplementedError('only applies to exact breakpoint calls')
        # get the reference positions for each breakpoint interval from the breakpointpair
        # convert this to the query intervals using the alignment
        # for each query interval calculate the read coverage as a pileup over the distance
        st = read.reference_start + 1
        end = read.reference_end
        if breakpoint.orient == ORIENT.LEFT:
            if breakpoint.start < st:
                return 0
            end = min(breakpoint.start, end)
        elif breakpoint.orient == ORIENT.RIGHT:
            if breakpoint.start > end:
                return 0
            st = max(st, breakpoint.start)
        qrange = _read.map_ref_range_to_query_range(read, Interval(st, end))
        return contig.remap_depth(qrange)


def query_coverage_interval(read):
    """
    Returns:
        :class:`~mavis.interval.Interval`: The portion of the original query sequence that is aligned by this read
    """
    seq = read.query_sequence
    st = 0
    end = len(seq) - 1
    if read.cigar[0][0] == CIGAR.S:
        st += read.cigar[0][1]
    if read.cigar[-1][0] == CIGAR.S:
        end -= read.cigar[-1][1]
    return Interval(st, end)


def convert_to_duplication(alignment, reference_genome):
    """
    Given a breakpoint call, tests if the untemplated sequences matches the preceding
    reference sequence. If it does this is annotated as a duplication and the new
    breakpoint pair is returned. If not, then the original breakpoint pair is returned
    """
    # assumes that events with deletions cannot be duplications
    if alignment.untemplated_seq and not alignment.interchromosomal and alignment.break1.end == alignment.break2.start - 1:
        # must be more than half the length or better to call it an insertion
        for dup_len in reversed(range(len(alignment.untemplated_seq) // 2 + 1, len(alignment.untemplated_seq) + 1)):
            refseq = reference_genome[alignment.break1.chr].seq[alignment.break1.start - dup_len:alignment.break1.start]
            refseq = str(refseq).upper()
            if refseq != alignment.untemplated_seq[:dup_len]:
                continue

            result = SplitAlignment(
                Breakpoint(alignment.break2.chr, alignment.break2.start - dup_len, orient=ORIENT.RIGHT, strand=alignment.break2.strand),
                Breakpoint(alignment.break1.chr, alignment.break1.start, orient=ORIENT.LEFT, strand=alignment.break1.strand),
                untemplated_seq=alignment.untemplated_seq[dup_len:],
                opposing_strands=alignment.opposing_strands,
                data=alignment.data,
                read1=alignment.read1,
                read2=alignment.read2
            )
            return result
    return alignment


def call_read_events(read, secondary_read=None):
    """
    Given a read, return breakpoint pairs representing all putative events
    """
    events = []
    reference_pos = read.reference_start
    query_pos = 0
    curr_event = None
    for state, freq in read.cigar:
        if state in _cigar.QUERY_ALIGNED_STATES & _cigar.REFERENCE_ALIGNED_STATES:
            query_pos += freq
            reference_pos += freq
            if curr_event:
                events.append(curr_event)
                curr_event = None
        elif state == CIGAR.S:
            query_pos += freq
            if curr_event:
                events.append(curr_event)
                curr_event = None
        elif state in _cigar.REFERENCE_ALIGNED_STATES:  # del
            if curr_event:
                ref_start, delsize, insseq = curr_event
                curr_event = (ref_start, delsize + freq, insseq)
            else:
                curr_event = (reference_pos, freq, '')
            reference_pos += freq
        elif state in _cigar.QUERY_ALIGNED_STATES:  # ins
            if curr_event:
                ref_start, delsize, insseq = curr_event
                curr_event = (ref_start, delsize, insseq + read.query_sequence[query_pos: query_pos + freq])
            else:
                curr_event = (reference_pos, 0, read.query_sequence[query_pos: query_pos + freq])
            query_pos += freq
        elif state != CIGAR.H:
            raise NotImplementedError('Should never happen. Invalid cigar state is not reference or query aligned', state)
    if curr_event:
        events.append(curr_event)
    result = []
    strand = STRAND.NEG if read.is_reverse else STRAND.POS
    for ref_start, delsize, insseq in events:
        bpp = SplitAlignment(
            Breakpoint(read.reference_name, ref_start, orient=ORIENT.LEFT, strand=strand),
            Breakpoint(read.reference_name, ref_start + delsize + 1, orient=ORIENT.RIGHT, strand=strand),
            untemplated_seq=insseq, read1=read, read2=secondary_read
        )
        result.append(bpp)
    return result


def read_breakpoint(read):
    """
    convert a given read to a single breakpoint
    """
    start_softclipping = read.cigar[0][1] if read.cigar[0][0] in _cigar.CLIPPING_STATE else 0
    end_softclipping = read.cigar[-1][1] if read.cigar[-1][0] in _cigar.CLIPPING_STATE else 0

    if start_softclipping == end_softclipping:
        raise AssertionError('softclipping is equal and therefore orientation cannot be determined')
    elif start_softclipping + end_softclipping == 0:
        raise AssertionError('no softclipping, therefore no orientation can be called')
    elif start_softclipping > end_softclipping:
        return Breakpoint(
            read.reference_name, read.reference_start + 1,
            orient=ORIENT.RIGHT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS,
            seq=read.query_alignment_sequence
        )
    else:
        return Breakpoint(
            read.reference_name, read.reference_end,
            orient=ORIENT.LEFT,
            strand=STRAND.NEG if read.is_reverse else STRAND.POS,
            seq=read.query_alignment_sequence
        )


def call_paired_read_event(read1, read2):
    """
    For a given pair of reads call all applicable events. Assume there is a major
    event from both reads and then call indels from the individual reads
    """
    # sort the reads so that we are calling consistently
    break1 = read_breakpoint(read1)
    break2 = read_breakpoint(read2)
    if break2.key < break1.key:
        break1, break2 = break2, break1
        read1, read2 = read2, read1
    r1_query_cover = query_coverage_interval(read1)
    r2_query_cover = query_coverage_interval(read2)

    # if the reads are on the same strand, then the query sequence should match
    if read1.is_reverse == read2.is_reverse:
        assert read1.query_sequence == read2.query_sequence
    else:
        assert read1.query_sequence == reverse_complement(read2.query_sequence)
        length = len(read1.query_sequence) - 1
        r2_query_cover = Interval(length - r2_query_cover.end, length - r2_query_cover.start)

    overlap = 0
    if Interval.overlaps(r1_query_cover, r2_query_cover):
        # adjust the second read to remove the overlapping query region
        overlap = len(r1_query_cover & r2_query_cover)

    if break2.orient == ORIENT.RIGHT:
        break2.start += overlap
        break2.end += overlap
        break2.seq = read2.query_alignment_sequence[overlap:]
        assert break2.start <= break2.end
    else:
        seq = read2.query_alignment_sequence
        if overlap > 0:
            seq = read2.query_alignment_sequence[:-1 * overlap]
        break2.start -= overlap
        break2.end -= overlap
        break2.seq = seq
    assert break2.start <= break2.end

    # now check for untemplated sequence
    untemplated_seq = ''
    dist = Interval.dist(r1_query_cover, r2_query_cover)
    if dist > 0:  # query coverage for read1 is after query coverage for read2
        untemplated_seq = read1.query_sequence[r2_query_cover[1] + 1:r1_query_cover[0]]
    elif dist < 0:  # query coverage for read2 is after query coverage for read1
        untemplated_seq = read1.query_sequence[r1_query_cover[1] + 1:r2_query_cover[0]]
    else:  # query coverage overlaps
        pass
    return SplitAlignment(break1, break2, untemplated_seq=untemplated_seq, read1=read1, read2=read2)


def align_sequences(
    sequences,
    input_bam_cache,
    reference_genome,
    aligner,
    aligner_reference,
    aligner_output_file='aligner_out.temp',
    aligner_fa_input_file='aligner_in.fa',
    aligner_output_log='aligner_out.log',
    blat_limit_top_aln=25,
    blat_min_identity=0.7,
    clean_files=True,
    log=devnull,
    **kwargs
):
    """
    calls the alignment tool and parses the return output for a set of sequences

    Args:
        sequences (dict of str to str): dictionary of sequences by name
        input_bam_cache (BamCache): bam cache to be used as a template for reading the alignments
        reference_genome: the reference genome
        aligner (SUPPORTED_ALIGNER): the name of the aligner to be used
        aligner_reference (str): path to the aligner reference file
    """
    try:
        # write the input sequences to a fasta file
        count = 1
        with open(aligner_fa_input_file, 'w') as fh:
            for name, seq in sorted(sequences.items()):
                fh.write('>' + name + '\n' + seq + '\n')
                count += 1
        if not sequences:
            return []

        log('will use', aligner, 'to align', len(sequences), 'unique sequences', time_stamp=False)

        # call the aligner using subprocess
        if aligner == SUPPORTED_ALIGNER.BLAT:
            from .blat import process_blat_output
            # call the aligner using subprocess
            blat_min_identity *= 100
            blat_options = kwargs.pop(
                'align_options', '-stepSize=5 -repMatch=2253 -minScore=0 -minIdentity={0}'.format(blat_min_identity))
            # call the blat subprocess
            # will raise subprocess.CalledProcessError if non-zero exit status
            # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
            command = ' '.join([
                SUPPORTED_ALIGNER.BLAT, aligner_reference, aligner_fa_input_file, aligner_output_file,
                '-out=pslx', '-noHead', blat_options])
            log('writing aligner logging to:', aligner_output_log, time_stamp=False)
            with open(aligner_output_log, 'w') as log_fh:
                log_fh.write('>>> {}\n'.format(command))
                subprocess.check_call(command, shell=True, stdout=log_fh, stderr=log_fh)
            return process_blat_output(
                input_bam_cache=input_bam_cache,
                query_id_mapping=sequences,
                reference_genome=reference_genome,
                aligner_output_file=aligner_output_file,
                blat_limit_top_aln=blat_limit_top_aln
            )

        elif aligner == SUPPORTED_ALIGNER.BWA_MEM:
            align_options = kwargs.get('align_options', '')
            command = '{} {} {} -Y {}'.format(aligner, aligner_reference, aligner_fa_input_file, align_options)
            log('writing aligner logging to:', aligner_output_log, time_stamp=False)
            with open(aligner_output_log, 'w') as log_fh, open(aligner_output_file, 'w') as aligner_output_fh:
                log_fh.write('>>> {}\n'.format(command))
                subprocess.check_call(command, stdout=aligner_output_fh, shell=True, stderr=log_fh)

            with pysam.AlignmentFile(aligner_output_file, 'r', check_sq=bool(len(sequences))) as samfile:
                reads_by_query = {}
                for read in samfile.fetch():
                    if read.is_unmapped:
                        continue
                    read = _read.SamRead.copy(read)
                    try:
                        read.reference_id = input_bam_cache.reference_id(read.reference_name)
                    except KeyError:
                        log('dropping alignment (unknown reference)', read.reference_name, time_stamp=False)
                    else:
                        if read.is_paired:
                            read.next_reference_id = input_bam_cache.reference_id(read.next_reference_name)
                        read.cigar = _cigar.recompute_cigar_mismatch(read, reference_genome[read.reference_name])
                        query_seq = sequences[read.query_name]
                        reads_by_query.setdefault(query_seq, []).append(read)
            return reads_by_query
        else:
            raise NotImplementedError('unsupported aligner', aligner)
    finally:
        # clean up
        if clean_files:
            for outputfile in [aligner_output_file, aligner_fa_input_file, aligner_output_log]:
                if os.path.exists(outputfile):
                    try:
                        os.remove(outputfile)
                    except OSError as err:
                        warnings.warn(repr(err))


def select_contig_alignments(evidence, reads_by_query):
    """
    standardize/simplify reads and filter bad/irrelevant alignments
    adds the contig alignments to the contigs
    """
    for contig in evidence.contigs:
        std_reads = set()
        alignments = []
        for raw_read in reads_by_query.get(contig.seq, []):
            if raw_read.reference_name not in {evidence.break1.chr, evidence.break2.chr}:
                continue
            read = evidence.standardize_read(raw_read)
            read.cigar = _cigar.merge_internal_events(
                read.cigar,
                inner_anchor=evidence.contig_aln_merge_inner_anchor,
                outer_anchor=evidence.contig_aln_merge_outer_anchor
            )
            read = evidence.standardize_read(read)  # genome needs to merge first, trans needs to standard first

            for single_alignment in call_read_events(read):
                single_alignment = convert_to_duplication(single_alignment, evidence.reference_genome)
                if single_alignment.break1.chr in {evidence.break1.chr, evidence.break2.chr}:
                    if any([
                        (single_alignment.break1 | single_alignment.break2) & evidence.outer_window1,
                        (single_alignment.break1 | single_alignment.break2) & evidence.outer_window2
                    ]):
                        alignments.append(single_alignment)

            std_reads.add(_read.convert_events_to_softclipping(
                read, evidence.break1.orient,
                max_event_size=evidence.contig_aln_max_event_size,
                min_anchor_size=evidence.contig_aln_min_anchor_size
            ))
            if evidence.break1.orient == evidence.break2.orient:
                continue
            std_reads.add(_read.convert_events_to_softclipping(
                read, evidence.break2.orient,
                max_event_size=evidence.contig_aln_max_event_size,
                min_anchor_size=evidence.contig_aln_min_anchor_size
            ))

        for read1, read2 in itertools.combinations(std_reads, 2):
            try:
                paired_event = call_paired_read_event(read1, read2)
            except AssertionError:
                continue
            if all([
                BreakpointPair.classify(paired_event) & BreakpointPair.classify(evidence),
                paired_event.break1.chr == evidence.break1.chr,
                paired_event.break2.chr == evidence.break2.chr,
                paired_event.break1 & evidence.outer_window1,
                paired_event.break2 & evidence.outer_window2
            ]):
                alignments.append(paired_event)
                for event in call_read_events(read1, read2) + call_read_events(read2, read1):
                    event = convert_to_duplication(event, evidence.reference_genome)
                    alignments.append(event)
        filtered_alignments = set()
        for alignment in sorted(alignments, key=lambda x: (x.read2 is None, x.score()), reverse=True):
            if any([
                alignment.query_consumption() < evidence.contig_aln_min_query_consumption,
                alignment.read2 is not None and alignment.query_overlap_extension() < evidence.contig_aln_min_extend_overlap,
                alignment in filtered_alignments,
                alignment.score() < evidence.contig_aln_min_score,
                alignment.mapping_quality() == Interval(0)
            ]):
                continue
            filtered_alignments.add(alignment)
        contig.alignments.update(filtered_alignments)
