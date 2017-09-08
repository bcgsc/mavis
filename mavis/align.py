"""
Should take in a sam file from a aligner like bwa aln or bwa mem and convert it into a

"""
import itertools
import pysam
import subprocess
import warnings
import os
from copy import copy
from .constants import COLUMNS, SVTYPE, CIGAR, reverse_complement, ORIENT
from .bam import cigar as cigar_tools
from .bam import read as read_tools
from .interval import Interval
from .util import devnull
from .breakpoint import BreakpointPair
from .error import InvalidRearrangement
from vocab import Vocab


SUPPORTED_ALIGNER = Vocab(BWA_MEM='bwa mem', BLAT='blat')


class SplitAlignment:
    def __init__(self, read1, read2=None):
        if read2 is not None and any([
            read1.reference_name > read2.reference_name,
            read1.reference_name == read2.reference_name and read1.reference_start > read2.reference_start
        ]):
            read1, read2 = read2, read1

        self.read1 = read1
        self.read2 = read2
        if self.read2 is not None:
            read2_seq = read2.query_sequence if not self.opposing_strands else reverse_complement(read2.query_sequence)
            if read1.query_sequence != read2_seq:
                raise ValueError('valid split alignments must share the same (or reverse complement) sequence')
        self.query_sequence = self.read1.query_sequence

        # TODO: if the reads 'share' query sequence overlap, remove this from the second alignment

    def __getitem__(self, index):
        if isinstance(index, int):
            if index == 0:
                return self.read1
            elif index == 1:
                return self.read2
            else:
                raise IndexError('index out of bounds', index)
        else:
            raise ValueError('index must be an integer')

    @property
    def opposing_strands(self):
        if self.read2 is None:
            return False
        else:
            return self.read1.is_reverse != self.read2.is_reverse

    def query_coverage_read1(self):
        qc1 = query_coverage_interval(self.read1)
        return qc1

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
        else:
            return self.query_coverage_read1() | self.query_coverage_read2()

    def query_consumption(self):
        """
        fraction of the query sequence which is aligned (everything not soft-clipped) in either alignment
        """
        if self.read2 is None or Interval.overlaps(self.query_coverage_read1(), self.query_coverage_read2()):
            return len(self.query_coverage()) / len(self.query_sequence)
        else:
            return (len(self.query_coverage_read1()) + len(self.query_coverage_read2())) / len(self.query_sequence)

    def query_overlap_extension(self):
        if self.read2 is not None:
            max_init_overlap = max(len(self.query_coverage_read1()), len(self.query_coverage_read2()))
            total_overlap = len(self.query_coverage()) - max_init_overlap
            return total_overlap
        else:
            return 0

    def score(self, consec_bonus=10):
        def score_matches(cigar):
            return sum([v + (v - 1) * consec_bonus for s, v in cigar if s == CIGAR.EQ])
        score = score_matches(self.read1.cigar)
        qlen = len(self.query_coverage())
        if self.read2 is not None:
            score += score_matches(self.read2.cigar)
            if Interval.overlaps(self.query_coverage_read1(), self.query_coverage_read2()):
                qlen += len(self.query_coverage_read1() & self.query_coverage_read2())
        return score / (qlen + (qlen - 1) * consec_bonus)

    def select_supporting_alignments(
        bpp, alignments,
        min_query_consumption,
        min_extend_overlap,
        max_event_size,
        min_anchor_size,
        merge_inner_anchor,
        merge_outer_anchor
    ):
        """
        give a breakpoint pair and a set of alignments for contigs associated with the given pair,
        alignments are paired (some events cannot be represented with a single bamfile alignment)
        and the most appropriate alignments supporting the breakpoint pair are selected and returned
        """
        # now for each bpp assign an alignment to each contig
        putative_alignments = []
        putative_event_types = set(bpp.putative_event_types())
        if {SVTYPE.INS, SVTYPE.DUP} & putative_event_types:
            putative_event_types = putative_event_types | {SVTYPE.INS, SVTYPE.DUP}

        # for events on the same template and strand we expect to find a single contig alignment
        if not bpp.interchromosomal and not bpp.opposing_strands:
            for read in alignments:
                try:
                    aln = SplitAlignment(read)
                except ValueError:
                    continue
                # if it covers both breakpoints add to putative alignments
                ref_cover = Interval(read.reference_start, read.reference_end - 1)
                if all([
                    aln.read1.reference_name == bpp.break1.chr,
                    Interval.overlaps(bpp.outer_window1, ref_cover),
                    Interval.overlaps(bpp.outer_window2, ref_cover)
                ]):
                    if aln.query_consumption() < min_query_consumption:
                        continue
                    # split the continuous alignment, assume ins/dup or indel
                    ins = sum([v for c, v in read.cigar if c == CIGAR.I] + [0])
                    dln = sum([v for c, v in read.cigar if c in [CIGAR.D, CIGAR.N]] + [0])

                    aln.read1 = copy(read)
                    aln.read1.cigar = cigar_tools.merge_internal_events(aln.read1.cigar, merge_inner_anchor, merge_outer_anchor)

                    for event_type in putative_event_types:
                        if event_type in {SVTYPE.INS, SVTYPE.DUP} and ins > 0 and ins > dln:
                            putative_alignments.append(aln)
                        elif event_type in {SVTYPE.DEL, SVTYPE.INV} and dln > 0 and dln > ins:
                            putative_alignments.append(aln)

        # don't use reads in combined alignments if they have already been assigned in a single alignment
        combo_prohibited = [x for x, y in putative_alignments]
        for read1, read2 in itertools.combinations([x for x in alignments if x not in combo_prohibited], 2):
            # do they overlap both breakpoints
            try:
                aln = SplitAlignment(read1, read2)
            except ValueError:
                continue
            if aln.read1.reference_name != bpp.break1.chr or aln.read2.reference_name != bpp.break2.chr:
                continue
            aln.read1 = read_tools.convert_events_to_softclipping(
                aln.read1, bpp.break1.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
            aln.read2 = read_tools.convert_events_to_softclipping(
                aln.read2, bpp.break2.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
            # check that the coverage intervals overlap the event windows
            if any([
                not Interval.overlaps((aln.read1.reference_start + 1, aln.read1.reference_end), bpp.outer_window1),
                not Interval.overlaps((aln.read2.reference_start + 1, aln.read2.reference_end), bpp.outer_window2)
            ]):
                continue
            # reads should have unique reference overlap
            if not bpp.interchromosomal and aln.read1.reference_end > aln.read2.reference_end:
                continue
            # check that the combination extends the amount of the initial query sequence we consume
            qc = len(aln.query_coverage())
            if any([
                len(aln.query_coverage_read1()) >= qc or len(aln.query_coverage_read2()) >= qc,
                aln.query_consumption() < min_query_consumption,
                aln.read2 is not None and aln.query_overlap_extension() < min_extend_overlap
            ]):
                continue

            try:
                call = BreakpointPair.call_breakpoint_pair(read1, read2)
                if not set(BreakpointPair.classify(call)) & putative_event_types:
                    continue
            except (InvalidRearrangement, AssertionError):
                continue
            putative_alignments.append(aln)
        return putative_alignments

    @staticmethod
    def breakpoint_contig_remapped_depth(breakpoint, contig, read):
        if breakpoint.chr != read.reference_name:
            raise AssertionError('breakpoint chromosome does not match read reference', breakpoint, read.reference_name)
        if len(breakpoint) > 1:
            raise NotImplementedError('only applies to exact breakpoint calls')
        # get the reference positions for each breakpoint interval from the breakpointpair
        # convert this to the query intervals using the alignment
        # for each query interval calculate the read coverage as a pileup over the distance
        s = read.reference_start + 1
        t = read.reference_end
        if breakpoint.orient == ORIENT.LEFT:
            if breakpoint.start < s:
                return 0
            t = min(breakpoint.start, t)
        elif breakpoint.orient == ORIENT.RIGHT:
            if breakpoint.start > t:
                return 0
            s = max(s, breakpoint.start)
        qrange = read_tools.map_ref_range_to_query_range(read, Interval(s, t))
        return contig.remap_depth(qrange)


def query_coverage_interval(read):
    """
    Returns:
        :class:`~mavis.interval.Interval`: The portion of the original query sequence that is aligned by this read
    """
    seq = read.query_sequence
    s = 0
    t = len(seq) - 1
    if read.cigar[0][0] == CIGAR.S:
        s += read.cigar[0][1]
    if read.cigar[-1][0] == CIGAR.S:
        t -= read.cigar[-1][1]
    return Interval(s, t)


def align_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome,
        aligner,
        aligner_reference,
        aligner_output_file='aligner_out.temp',
        aligner_fa_input_file='aligner_in.fa',
        blat_min_identity=0.7,
        contig_aln_min_query_consumption=0.5,
        contig_aln_max_event_size=50,
        contig_aln_min_anchor_size=50,
        contig_aln_merge_inner_anchor=20,
        contig_aln_merge_outer_anchor=20,
        blat_limit_top_aln=25,
        is_protein=False,
        min_extend_overlap=10,
        clean_files=True,
        log=devnull,
        **kwargs):
    """
    given a set of contigs, call the aligner from the command line and adds the results to the contigs
    associated with each Evidence object

    """
    if is_protein:
        raise NotImplementedError('currently does not support aligning protein sequences')

    try:
        # write the input sequences to a fasta file
        query_id_mapping = {}
        sequences = set()
        count = 1
        ev_by_seq = {}
        for e in evidence:
            for c in e.contigs:
                sequences.add(c.seq)
                ev_by_seq.setdefault(c.seq, []).append(e.data.get(COLUMNS.cluster_id, None))

        with open(aligner_fa_input_file, 'w') as fh:
            for seq in sequences:
                n = 'seq{}'.format(count)
                log(n, [x for x in ev_by_seq[seq] if x is not None])
                query_id_mapping[n] = seq
                fh.write('>' + n + '\n' + seq + '\n')
                count += 1
        if len(sequences) == 0:
            return

        log('will use', aligner, 'to align', len(sequences), 'unique sequences', time_stamp=False)

        # call the aligner using subprocess
        if aligner == SUPPORTED_ALIGNER.BLAT:
            from .blat import process_blat_output
            # call the aligner using subprocess
            blat_min_identity *= 100
            blat_options = kwargs.pop(
                'blat_options', ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(blat_min_identity)])
            # call the blat subprocess
            # will raise subprocess.CalledProcessError if non-zero exit status
            # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
            log([SUPPORTED_ALIGNER.BLAT, aligner_reference,
                 aligner_fa_input_file, aligner_output_file, '-out=pslx', '-noHead'] + blat_options)
            subprocess.check_output([
                SUPPORTED_ALIGNER.BLAT, aligner_reference,
                aligner_fa_input_file, aligner_output_file, '-out=pslx', '-noHead'] + blat_options)
            reads_by_query = process_blat_output(
                INPUT_BAM_CACHE=INPUT_BAM_CACHE,
                query_id_mapping=query_id_mapping,
                reference_genome=reference_genome,
                aligner_output_file=aligner_output_file,
                blat_limit_top_aln=blat_limit_top_aln,
                is_protein=is_protein
            )

        elif aligner == SUPPORTED_ALIGNER.BWA_MEM:
            command = '{} {} {} -Y'.format(aligner, aligner_reference, aligner_fa_input_file)
            log(command)  # for bwa
            with open(aligner_output_file, 'w') as f:
                subprocess.call(command, stdout=f, shell=True)

            with pysam.AlignmentFile(aligner_output_file, 'r') as samfile:
                reads_by_query = {}
                for read in samfile.fetch():
                    read = read_tools.SamRead.copy(read)
                    read.reference_id = INPUT_BAM_CACHE.reference_id(read.reference_name)
                    if read.is_paired:
                        read.next_reference_id = INPUT_BAM_CACHE.reference_id(read.next_reference_name)
                    read.cigar = cigar_tools.recompute_cigar_mismatch(read, reference_genome[read.reference_name])
                    query_seq = query_id_mapping[read.query_name]
                    reads_by_query.setdefault(query_seq, []).append(read)
        else:
            raise NotImplementedError('unsupported aligner', aligner)

        for e in evidence:
            for contig in e.contigs:
                aln = reads_by_query.get(contig.seq, [])
                putative_alignments = SplitAlignment.select_supporting_alignments(
                    e, aln,
                    min_extend_overlap=min_extend_overlap,
                    min_query_consumption=contig_aln_min_query_consumption,
                    min_anchor_size=contig_aln_min_anchor_size,
                    max_event_size=contig_aln_max_event_size,
                    merge_inner_anchor=contig_aln_merge_inner_anchor,
                    merge_outer_anchor=contig_aln_merge_outer_anchor
                )
                contig.alignments.extend(putative_alignments)
    finally:
        # clean up
        if clean_files:
            for f in [aligner_output_file, aligner_fa_input_file]:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                    except OSError as e:
                        warnings.warn(repr(e))
