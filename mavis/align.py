"""
Should take in a sam file from a aligner like bwa aln or bwa mem and convert it into a

"""
import itertools
import pysam
import subprocess
import warnings
import os
from copy import copy
from .constants import COLUMNS, SVTYPE, CIGAR, PYSAM_READ_FLAGS, reverse_complement
from .bam import cigar as cigar_tools
from .bam import read as read_tools
from .interval import Interval
from .util import devnull
from .breakpoint import BreakpointPair
from .error import InvalidRearrangement


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


def paired_alignment_score(read1, read2=None):
    score = read_tools.calculate_alignment_score(read1) * query_coverage_interval(read1).length()
    if read2 is not None:
        qci1 = query_coverage_interval(read1)
        qci2 = query_coverage_interval(read2)
        total_len = qci1.length() + qci2.length()
        s1 = read_tools.calculate_alignment_score(read1)
        s2 = read_tools.calculate_alignment_score(read2)
        avg = s1 * qci1.length() / total_len + s2 * qci2.length() / total_len
        score = avg * len(qci1 | qci2)
    return score


def select_paired_alignments(
    bpp, aligned_contigs,
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
        for read in aligned_contigs:
            # if it covers both breakpoints add to putative alignments
            read_cover = Interval(read.reference_start, read.reference_end - 1)
            if all([
                read.reference_name == bpp.break1.chr,
                Interval.overlaps(bpp.outer_window1, read_cover),
                Interval.overlaps(bpp.outer_window2, read_cover)
            ]):
                # split the continuous alignment, assume ins/dup or indel
                ins = sum([v for c, v in read.cigar if c == CIGAR.I] + [0])
                dln = sum([v for c, v in read.cigar if c in [CIGAR.D, CIGAR.N]] + [0])
                consume = len(query_coverage_interval(read)) / len(read.query_sequence)
                read = copy(read)
                read.cigar = cigar_tools.merge_internal_events(read.cigar, merge_inner_anchor, merge_outer_anchor)
                if consume < min_query_consumption:
                    continue

                for event_type in putative_event_types:
                    if event_type in {SVTYPE.INS, SVTYPE.DUP} and ins > 0 and ins > dln:
                        putative_alignments.append((read, None))
                    elif event_type in {SVTYPE.DEL, SVTYPE.INV} and dln > 0 and dln > ins:
                        putative_alignments.append((read, None))

    # don't use reads in combined alignments if they have already been assigned in a single alignment
    combo_prohibited = [x for x, y in putative_alignments]

    for read1, read2 in itertools.combinations([x for x in aligned_contigs if x not in combo_prohibited], 2):
        # do they overlap both breakpoints
        if any([
            read1.reference_name > read2.reference_name,
            read1.reference_name == read2.reference_name and read1.reference_start > read2.reference_start
        ]):
            read1, read2 = read2, read1

        if read1.reference_name != bpp.break1.chr or read2.reference_name != bpp.break2.chr:
            continue
        read1 = read_tools.convert_events_to_softclipping(
            read1, bpp.break1.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
        read2 = read_tools.convert_events_to_softclipping(
            read2, bpp.break2.orient, max_event_size=max_event_size, min_anchor_size=min_anchor_size)
        # check that the coverage intervals overlap
        if any([
            not Interval.overlaps((read1.reference_start + 1, read1.reference_end), bpp.outer_window1),
            not Interval.overlaps((read2.reference_start + 1, read2.reference_end), bpp.outer_window2)
        ]):
            continue
        # reads should have unique reference overlap
        if not bpp.interchromosomal and read1.reference_end > read2.reference_end:
            continue
        # check that the combination extends the amount of the initial query sequence we consume
        query_cover1 = query_coverage_interval(read1)
        query_cover2 = query_coverage_interval(read2)

        if read2.query_sequence != read1.query_sequence:  # alignments aligned to opposite strands
            if not bpp.opposing_strands:
                continue
            if read2.query_sequence != reverse_complement(read1.query_sequence):
                continue
            l = len(read1.query_sequence) - 1
            query_cover2 = Interval(l - query_cover2.end, l - query_cover2.start)
        elif bpp.opposing_strands:
            continue

        consume = len(query_cover1 | query_cover2)
        if not Interval.overlaps(query_cover1, query_cover2):
            consume = len(query_cover1) + len(query_cover2)

        if consume / len(read1.query_sequence) < min_query_consumption:
            continue
        if consume - len(query_cover1) < min_extend_overlap or consume - len(query_cover2) < min_extend_overlap:
            continue

        try:
            call = BreakpointPair.call_breakpoint_pair(read1, read2)
            if not set(BreakpointPair.classify(call)) & putative_event_types:
                continue
        except (InvalidRearrangement, AssertionError) as err:
            continue
        putative_alignments.append((read1, read2))
    return putative_alignments


def align_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome,
        blat_2bit_reference='/home/pubseq/genomes/Homo_sapiens/GRCh37/blat/hg19.2bit',
        blat_min_percent_of_max_score=0.8,
        blat_min_identity=0.7,
        aligner='blat',
        aligner_options='mem',
        aligner_output_file='blat_out.pslx',
        aligner_fa_input_file='blat_in.fa',
        aligner_reference='/projects/trans_scratch/references/genomes/transabyss/bwamem-0.7.10/hg19a.fa',
        contig_aln_min_query_consumption=0.5,
        contig_aln_max_event_size=50,
        contig_aln_min_anchor_size=50,
        contig_aln_merge_inner_anchor=20,
        contig_aln_merge_outer_anchor=20,
        is_protein=False,
        min_extend_overlap=10,
        pair_scoring_function=paired_alignment_score,
        clean_files=True,
        log=devnull,
        **kwargs):
    """
    given a set of contigs, call the aligner from the command line and adds the results to the contigs
    asscoated with each Evidence object

    """
    if is_protein:
        raise NotImplementedError('currently does not support aligning protein sequences')
#    aligner_options = kwargs.pop(        'align_options')

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
        if aligner == "blat":
            from .blat import Blat
            blat_min_identity *= 100
            blat_options = kwargs.pop(
                 'blat_options', ["-stepSize=5", "-repMatch=2253", "-minScore=0", "-minIdentity={0}".format(blat_min_identity)])
            # call the blat subprocess
            # will raise subprocess.CalledProcessError if non-zero exit status
            # parameters from https://genome.ucsc.edu/FAQ/FAQblat.html#blat4
            log(['blat', blat_2bit_reference,
                 aligner_fa_input_file, aligner_output_file, '-out=pslx', '-noHead'] + blat_options)
            subprocess.check_output([
                    'blat', blat_2bit_reference,
                    aligner_fa_input_file, aligner_output_file, '-out=pslx', '-noHead'] + blat_options)

            header, rows = Blat.read_pslx(aligner_output_file, query_id_mapping, is_protein=is_protein)

            # split the rows by query id
            rows_by_query = {}
            for row in rows:
                if row['qname'] not in rows_by_query:
                    rows_by_query[row['qname']] = []
                rows_by_query[row['qname']].append(row)

            reads_by_query = {}
            for s in sequences:
                reads_by_query[s] = []
            for query_id, rows in rows_by_query.items():
                query_seq = query_id_mapping[query_id]
                # filter on percent id
                filtered_rows = [row for row in rows if round(row['percent_ident'], 0) >= blat_min_identity]

                # filter on score
                scores = sorted([r['score'] for r in rows])
                max_score = scores[-1]
                filtered_rows.sort(key=lambda x: x['score'], reverse=True)
                reads = []
                for rank, row in enumerate(filtered_rows):
                    try:
                        read = Blat.pslx_row_to_pysam(row, INPUT_BAM_CACHE, reference_genome)
                        read.set_tag(PYSAM_READ_FLAGS.BLAT_SCORE, row['score'], value_type='i')
                        read.set_tag(PYSAM_READ_FLAGS.BLAT_ALIGNMENTS, len(filtered_rows), value_type='i')
                        read.set_tag(PYSAM_READ_FLAGS.BLAT_PMS, blat_min_percent_of_max_score, value_type='f')
                        read.set_tag(PYSAM_READ_FLAGS.BLAT_RANK, rank, value_type='i')
                        read.set_tag(PYSAM_READ_FLAGS.BLAT_PERCENT_IDENTITY, row['percent_ident'], value_type='f')
                        reads.append(read)
                    except KeyError as e:
                        warnings.warn(
                            'warning: reference template name not recognized {0}'.format(e))
                    except AssertionError as e:
                        warnings.warn('warning: invalid blat alignment: {}'.format(e))
                reads_by_query[query_seq] = reads

        else:
            log([aligner, aligner_options, aligner_reference, aligner_fa_input_file]) #for bwa
            with open(aligner_output_file, 'w') as f:
                subprocess.call([
                        aligner, aligner_options, aligner_reference, aligner_fa_input_file],
                                stdout=f)

            samfile = pysam.AlignmentFile(aligner_output_file, 'r')

            reads_by_query = {}
            for read in samfile.fetch():
                reads_by_query.setdefault(read.query_name, []).append(read)

        for e in evidence:
            for contig in e.contigs:
                aln = reads_by_query.get(contig.seq, [])
                putative_alignments = select_paired_alignments(
                    e, aln,
                    min_extend_overlap=min_extend_overlap,
                    min_query_consumption=contig_aln_min_query_consumption,
                    min_anchor_size=contig_aln_min_anchor_size,
                    max_event_size=contig_aln_max_event_size,
                    merge_inner_anchor=contig_aln_merge_inner_anchor,
                    merge_outer_anchor=contig_aln_merge_outer_anchor
                )
                if len(putative_alignments) == 0:
                    continue
                score_by_alignments = {}
                for read1, read2 in sorted(putative_alignments, key=lambda x: pair_scoring_function(x[0], x[1])):
                    score = pair_scoring_function(read1, read2)
                    score_by_alignments[(read1, read2)] = score
                max_score = max(list(score_by_alignments.values()))
                for pair in putative_alignments:
                    score = score_by_alignments[pair]
                    if score / max_score >= blat_min_percent_of_max_score:
                        contig.alignments.append(pair)
    finally:
        # clean up
        if clean_files:
            for f in [aligner_output_file, aligner_fa_input_file]:
                if os.path.exists(f):
                    try:
                        os.remove(f)
                    except OSError as e:
                        warnings.warn(repr(e))
