"""
Module which contains functions related to picking reads from a BAM file which support a putative event
"""
import pysam

from ..bam import cigar as _cigar
from ..bam import read as _read
from ..constants import (
    CIGAR,
    NA_MAPPING_QUALITY,
    ORIENT,
    PYSAM_READ_FLAGS,
    SVTYPE,
    reverse_complement,
)
from ..error import NotSpecifiedError
from ..interval import Interval
from ..util import logger
from .base import Evidence


def collect_split_read(evidence_bpp: Evidence, read: pysam.AlignedSegment, first_breakpoint: bool):
    """
    adds a split read if it passes the criteria filters and raises a warning if it does not

    Args:
        read: the read to add
        first_breakpoint: add to the first breakpoint (or second if false)
    Returns:
        bool:
            - True: the read was collected and stored in the current evidence object
            - False: the read was not collected
    Raises:
        NotSpecifiedError: if the breakpoint orientation is not specified
    """
    breakpoint = evidence_bpp.break1 if first_breakpoint else evidence_bpp.break2
    window = evidence_bpp.inner_window1 if first_breakpoint else evidence_bpp.inner_window2
    opposite_breakpoint = evidence_bpp.break2 if first_breakpoint else evidence_bpp.break1
    opposite_window = evidence_bpp.inner_window2 if first_breakpoint else evidence_bpp.inner_window1

    if read.cigar[0][0] != CIGAR.S and read.cigar[-1][0] != CIGAR.S:
        return False  # split read is not softclipped
    elif breakpoint.orient == ORIENT.LEFT and read.cigar[-1][0] != CIGAR.S:
        return False  # split read is not softclipped
    elif breakpoint.orient == ORIENT.RIGHT and read.cigar[0][0] != CIGAR.S:
        return False  # split read is not softclipped

    # the first breakpoint of a BreakpointPair is always the lower breakpoint
    # if this is being added to the second breakpoint then we'll need to check if the
    # read soft-clipping needs to be adjusted

    # check if the read falls within the evidence collection window
    # recall that pysam is 0-indexed and the window is 1-indexed
    if (
        read.reference_start > window.end - 1
        or read.reference_end < window.start
        or evidence_bpp.bam_cache.get_read_reference_name(read) != breakpoint.chr
    ):
        return False  # read not in breakpoint evidence window
    # can only enforce strand if both the breakpoint and the bam are stranded
    if evidence_bpp.stranded and evidence_bpp.bam_cache.stranded:
        strand = _read.sequenced_strand(
            read, strand_determining_read=evidence_bpp.strand_determining_read
        )
        if strand != breakpoint.strand:
            return False  # split read not on the appropriate strand
    unused = ''
    primary = ''
    clipped = ''
    if breakpoint.orient == ORIENT.LEFT:
        unused = read.query_sequence[: read.query_alignment_start]
        primary = read.query_sequence[read.query_alignment_start : read.query_alignment_end]
        # end is exclusive in pysam
        clipped = read.query_sequence[read.query_alignment_end :]
    elif breakpoint.orient == ORIENT.RIGHT:
        clipped = read.query_sequence[: read.query_alignment_start]
        primary = read.query_sequence[read.query_alignment_start : read.query_alignment_end]
        unused = read.query_sequence[read.query_alignment_end :]
    else:
        raise NotSpecifiedError(
            'cannot assign split reads to a breakpoint where the orientation has not been specified'
        )
    if len(primary) + len(clipped) + len(unused) != len(read.query_sequence):
        raise AssertionError(
            'unused, primary, and clipped sequences should make up the original sequence',
            unused,
            primary,
            clipped,
            read.query_sequence,
            len(read.query_sequence),
        )

    if (
        len(primary) < evidence_bpp.config['validate.min_anchor_exact']
        or len(clipped) < evidence_bpp.config['validate.min_softclipping']
    ):
        # split read does not meet the minimum anchor criteria
        return False
    if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(
        PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
    ):
        read = evidence_bpp.standardize_read(read)
    # data quality filters
    if (
        _cigar.alignment_matches(read.cigar)
        >= evidence_bpp.config['validate.min_sample_size_to_apply_percentage']
        and _cigar.match_percent(read.cigar) < evidence_bpp.config['validate.min_anchor_match']
    ):
        return False  # too poor quality of an alignment
    if (
        _cigar.longest_exact_match(read.cigar) < evidence_bpp.config['validate.min_anchor_exact']
        and _cigar.longest_fuzzy_match(
            read.cigar, evidence_bpp.config['validate.fuzzy_mismatch_number']
        )
        < evidence_bpp.config['validate.min_anchor_fuzzy']
    ):
        return False  # too poor quality of an alignment
    else:
        evidence_bpp.split_reads[0 if first_breakpoint else 1].add(read)

    # try mapping the soft-clipped portion to the other breakpoint
    w = (opposite_window[0], opposite_window[1])
    opposite_breakpoint_ref = evidence_bpp.reference_genome[opposite_breakpoint.chr].seq[
        w[0] - 1 : w[1]
    ]

    # figure out how much of the read must match when remapped
    min_match_tgt = read.cigar[-1][1] if breakpoint.orient == ORIENT.LEFT else read.cigar[0][1]
    min_match_tgt = min(
        min_match_tgt * evidence_bpp.config['validate.min_anchor_match'], min_match_tgt - 1
    ) / len(read.query_sequence)
    if not evidence_bpp.opposing_strands:  # same strand
        sc_align = _read.nsb_align(
            opposite_breakpoint_ref,
            read.query_sequence,
            min_consecutive_match=evidence_bpp.config['validate.min_anchor_exact'],
            min_match=min_match_tgt,
            min_overlap_percent=min_match_tgt,
        )  # split half to this side

        for alignment in sc_align:
            alignment.flag = read.flag
        putative_alignments = sc_align
    else:
        # should align opposite the current read
        revcomp_sc_align = reverse_complement(read.query_sequence)
        revcomp_sc_align = _read.nsb_align(
            opposite_breakpoint_ref,
            revcomp_sc_align,
            min_consecutive_match=evidence_bpp.config['validate.min_anchor_exact'],
            min_match=min_match_tgt,
            min_overlap_percent=min_match_tgt,
        )

        for alignment in revcomp_sc_align:
            alignment.flag = read.flag ^ PYSAM_READ_FLAGS.REVERSE  # EXOR
        putative_alignments = revcomp_sc_align

    scores = []
    for alignment in putative_alignments:  # loop over the alignments
        alignment.flag = alignment.flag | PYSAM_READ_FLAGS.SUPPLEMENTARY
        # set this flag so we don't recompute the cigar multiple
        alignment.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
        # add information from the original read
        alignment.reference_start = w[0] - 1 + alignment.reference_start
        alignment._reference_name = (
            opposite_breakpoint.chr
        )  # must be set since not associated with an alignment file
        alignment.reference_id = evidence_bpp.bam_cache.reference_id(opposite_breakpoint.chr)
        alignment.query_name = read.query_name
        alignment.set_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT, 1, value_type='i')
        alignment.next_reference_start = read.next_reference_start
        alignment.next_reference_id = read.next_reference_id
        alignment.mapping_quality = NA_MAPPING_QUALITY
        try:
            cigar, offset = _cigar.extend_softclipping(
                alignment.cigar, evidence_bpp.config['validate.min_anchor_exact']
            )
            alignment.cigar = cigar
            alignment.reference_start = alignment.reference_start + offset
        except AttributeError:
            # if the matches section is too small you can't extend the
            # softclipping
            pass
        s = _cigar.score(alignment.cigar)
        alignment.cigar = _cigar.join(alignment.cigar)
        if alignment.reference_id == alignment.next_reference_id:
            # https://samtools.github.io/hts-specs/SAMv1.pdf
            # unsigned observed template length equals the number of bases from the leftmost
            # mapped base to the rightmost mapped base
            tlen = abs(alignment.reference_start - alignment.next_reference_start) + 1
            if alignment.reference_start < alignment.next_reference_start:
                alignment.template_length = tlen
            else:
                alignment.template_length = -1 * tlen
        else:
            alignment.template_length = 0
        if (
            _cigar.alignment_matches(alignment.cigar)
            >= evidence_bpp.config['validate.min_sample_size_to_apply_percentage']
            and _cigar.match_percent(alignment.cigar)
            < evidence_bpp.config['validate.min_anchor_match']
        ):
            continue
        if (
            _cigar.longest_exact_match(alignment.cigar)
            < evidence_bpp.config['validate.min_anchor_exact']
            and _cigar.longest_fuzzy_match(
                alignment.cigar, evidence_bpp.config['validate.fuzzy_mismatch_number']
            )
            < evidence_bpp.config['validate.min_anchor_fuzzy']
        ):
            continue
        if evidence_bpp.config['validate.max_sc_preceeding_anchor'] is not None:
            if opposite_breakpoint.orient == ORIENT.LEFT:
                if (
                    alignment.cigar[0][0] == CIGAR.S
                    and alignment.cigar[0][1]
                    > evidence_bpp.config['validate.max_sc_preceeding_anchor']
                ):
                    continue
            elif opposite_breakpoint.orient == ORIENT.RIGHT:
                if (
                    alignment.cigar[-1][0] == CIGAR.S
                    and alignment.cigar[-1][1]
                    > evidence_bpp.config['validate.max_sc_preceeding_anchor']
                ):
                    continue
        alignment.set_key()  # set the hash key before we add the read as evidence
        scores.append((s, _cigar.match_percent(alignment.cigar), alignment))

    scores = sorted(scores, key=lambda x: (x[0], x[1]), reverse=True) if scores else []

    if len(scores) > 1:
        if scores[0][0] != scores[1][0] and scores[0][1] != scores[1][1]:
            # not multimap, pick highest scoring alignment
            clipped = scores[0][2]
            evidence_bpp.split_reads[1 if first_breakpoint else 0].add(
                clipped
            )  # add to the opposite breakpoint
    elif len(scores) == 1:
        clipped = scores[0][2]
        evidence_bpp.split_reads[1 if first_breakpoint else 0].add(
            clipped
        )  # add to the opposite breakpoint
    return True


def collect_half_mapped(
    evidence_bpp: Evidence, read: pysam.AlignedSegment, mate: pysam.AlignedSegment
):
    """
    Args:
        read: the read to add
        mate: the unmapped mate

    Returns:
        bool:
            - True: the read was collected and stored in the current evidence object
            - False: the read was not collected
    Raises:
        AssertionError: if the mate is not unmapped
    """
    if not mate.is_unmapped:
        raise AssertionError('expected the mate to be unmapped')
    read_itvl = Interval(read.reference_start + 1, read.reference_end)
    added = False
    if evidence_bpp.break1.chr == read.reference_name and Interval.overlaps(
        evidence_bpp.outer_window1, read_itvl
    ):
        evidence_bpp.half_mapped[0].add(mate)
        added = True
    if evidence_bpp.break2.chr == read.reference_name and Interval.overlaps(
        evidence_bpp.outer_window2, read_itvl
    ):
        evidence_bpp.half_mapped[1].add(mate)
        added = True
    return added


def collect_flanking_pair(
    evidence_bpp: Evidence, read: pysam.AlignedSegment, mate: pysam.AlignedSegment
):
    """
    checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
    this event

    Args:
        read: the read to add
        mate: the mate

    Returns:
        bool:
            - True: the pair was collected and stored in the current evidence object
            - False: the pair was not collected
    Raises:
        ValueError: if the input reads are not a valid pair

    see [theory - types of flanking evidence](/background/theory/#types-of-flanking-evidence)
    """
    if read.is_unmapped or mate.is_unmapped:
        raise ValueError(
            'input reads must be a mapped and mated pair. One or both of the reads is unmapped'
        )
    elif read.query_name != mate.query_name:
        raise ValueError(
            'input reads must be a mapped and mated pair. The query names do not match'
        )
    elif abs(read.template_length) != abs(mate.template_length):
        raise ValueError(
            'input reads must be a mapped and mated pair. The template lengths (abs value) do not match',
            abs(read.template_length),
            abs(mate.template_length),
        )
    elif (
        read.mapping_quality < evidence_bpp.min_mapping_quality
        or mate.mapping_quality < evidence_bpp.min_mapping_quality
    ):
        return False  # do not meet the minimum mapping quality requirements
    # check that the references are right for the pair
    if evidence_bpp.interchromosomal:
        if read.reference_id == read.next_reference_id:
            return False
    elif read.reference_id != read.next_reference_id:
        return False

    if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(
        PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
    ):
        read = evidence_bpp.standardize_read(read)
    if not mate.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not mate.get_tag(
        PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
    ):
        mate = evidence_bpp.standardize_read(mate)
    # order the read pairs so that they are in the same order that we expect for the breakpoints
    if read.reference_id != mate.reference_id:
        if evidence_bpp.bam_cache.get_read_reference_name(
            read
        ) > evidence_bpp.bam_cache.get_read_reference_name(mate):
            read, mate = mate, read
    elif read.reference_start > mate.reference_start:
        read, mate = mate, read

    if (
        evidence_bpp.bam_cache.get_read_reference_name(read) != evidence_bpp.break1.chr
        or evidence_bpp.bam_cache.get_read_reference_name(mate) != evidence_bpp.break2.chr
    ):
        return False
    if any(
        [
            evidence_bpp.break1.orient == ORIENT.LEFT and read.is_reverse,
            evidence_bpp.break1.orient == ORIENT.RIGHT and not read.is_reverse,
            evidence_bpp.break2.orient == ORIENT.LEFT and mate.is_reverse,
            evidence_bpp.break2.orient == ORIENT.RIGHT and not mate.is_reverse,
        ]
    ):
        return False
    # check if this read falls in the first breakpoint window
    iread = Interval(read.reference_start + 1, read.reference_end)
    imate = Interval(mate.reference_start + 1, mate.reference_end)

    if evidence_bpp.stranded:
        strand1 = _read.sequenced_strand(read, evidence_bpp.strand_determining_read)
        strand2 = _read.sequenced_strand(mate, evidence_bpp.strand_determining_read)

        if strand1 != evidence_bpp.break1.strand or strand2 != evidence_bpp.break2.strand:
            return False

    for event_type in evidence_bpp.putative_event_types():

        # check that the pair orientation is correct
        if not _read.orientation_supports_type(read, event_type):
            continue

        # check that the fragment size is reasonable
        fragment_size = evidence_bpp.compute_fragment_size(read, mate)

        if event_type == SVTYPE.DEL:
            if fragment_size.end <= evidence_bpp.max_expected_fragment_size:
                continue
        elif event_type == SVTYPE.INS:
            if fragment_size.start >= evidence_bpp.min_expected_fragment_size:
                continue

        # check that the positions of the reads and the strands make sense
        if Interval.overlaps(iread, evidence_bpp.outer_window1) and Interval.overlaps(
            imate, evidence_bpp.outer_window2
        ):
            evidence_bpp.flanking_pairs.add((read, mate))
            return True

    return False


def collect_compatible_flanking_pair(
    evidence_bpp: Evidence,
    read: pysam.AlignedSegment,
    mate: pysam.AlignedSegment,
    compatible_type: str,
) -> bool:
    """
    checks if a given read meets the minimum quality criteria to be counted as evidence as stored as support for
    this event

    Args:
        read: the read to add
        mate: the mate
        compatible_type (SVTYPE): the type we are collecting for

    Returns:
        bool:
            - True: the pair was collected and stored in the current evidence object
            - False: the pair was not collected
    Raises:
        ValueError: if the input reads are not a valid pair

    Note:
        see [theory - types of flanking evidence](/background/theory/#compatible-flanking-pairs)
    """
    if (
        read.is_unmapped
        or mate.is_unmapped
        or read.query_name != mate.query_name
        or read.is_read1 == mate.is_read1
    ):
        raise ValueError('input reads must be a mapped and mated pair')
    if not evidence_bpp.compatible_window1:
        raise ValueError('compatible windows were not given')
    if evidence_bpp.interchromosomal:
        raise NotImplementedError('interchromosomal events do not have compatible flanking pairs')
    # check that the references are right for the pair
    if read.reference_id != read.next_reference_id:
        return False
    elif (
        read.mapping_quality < evidence_bpp.min_mapping_quality
        or mate.mapping_quality < evidence_bpp.min_mapping_quality
    ):
        return False
    if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(
        PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
    ):
        read = evidence_bpp.standardize_read(read)
    if not mate.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not mate.get_tag(
        PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
    ):
        mate = evidence_bpp.standardize_read(mate)
    # order the read pairs so that they are in the same order that we expect for the breakpoints
    if read.reference_start > mate.reference_start:
        read, mate = mate, read

    if (
        evidence_bpp.bam_cache.get_read_reference_name(read) != evidence_bpp.break1.chr
        or evidence_bpp.bam_cache.get_read_reference_name(mate) != evidence_bpp.break2.chr
    ):
        return False

    # check if this read falls in the first breakpoint window
    iread = Interval(read.reference_start + 1, read.reference_end)
    imate = Interval(mate.reference_start + 1, mate.reference_end)

    if evidence_bpp.stranded:
        strand1 = _read.sequenced_strand(read, evidence_bpp.strand_determining_read)
        strand2 = _read.sequenced_strand(mate, evidence_bpp.strand_determining_read)
        if strand1 != evidence_bpp.break1.strand or strand2 != evidence_bpp.break2.strand:
            return False

    # check that the pair orientation is correct
    if not _read.orientation_supports_type(read, compatible_type):
        return False

    # check that the fragment size is reasonable
    fragment_size = evidence_bpp.compute_fragment_size(read, mate)

    if compatible_type == SVTYPE.DEL:
        if fragment_size.end <= evidence_bpp.max_expected_fragment_size:
            return False
    elif compatible_type == SVTYPE.INS:
        if fragment_size.start >= evidence_bpp.min_expected_fragment_size:
            return False

    # check that the positions of the reads and the strands make sense
    if Interval.overlaps(iread, evidence_bpp.compatible_window1) and Interval.overlaps(
        imate, evidence_bpp.compatible_window2
    ):
        evidence_bpp.compatible_flanking_pairs.add((read, mate))
        return True

    return False


def collect_spanning_read(evidence_bpp: Evidence, read: pysam.AlignedSegment):
    """
    spanning read: a read covering BOTH breakpoints

    This is only applicable to small events. Do not need to look for soft clipped reads
    here since they will be collected already

    Args:
        read: the putative spanning read

    Returns:
        bool:
            - True: the read was collected and stored in the current evidence object
            - False: the read was not collected
    """
    if evidence_bpp.interchromosomal:
        return False
    elif not Interval.overlaps(evidence_bpp.inner_window1, evidence_bpp.inner_window2):
        # too far apart to call spanning reads
        return False

    if evidence_bpp.stranded:
        strand = _read.sequenced_strand(read, evidence_bpp.strand_determining_read)
        if strand != evidence_bpp.break1.strand and strand != evidence_bpp.break2.strand:
            return False

    combined = evidence_bpp.inner_window1 & evidence_bpp.inner_window2
    read_interval = Interval(read.reference_start + 1, read.reference_end)

    if Interval.overlaps(combined, read_interval):

        if not read.has_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR) or not read.get_tag(
            PYSAM_READ_FLAGS.RECOMPUTED_CIGAR
        ):

            read = evidence_bpp.standardize_read(read)
        # in the correct position, now determine if it can support the event types
        for event_type in evidence_bpp.putative_event_types():
            if event_type in [SVTYPE.DUP, SVTYPE.INS]:
                if CIGAR.I in [c[0] for c in read.cigar]:
                    evidence_bpp.spanning_reads.add(read)
                    return True
            elif event_type == SVTYPE.DEL:
                if CIGAR.D in [c[0] for c in read.cigar]:
                    evidence_bpp.spanning_reads.add(read)
                    return True
            elif event_type == SVTYPE.INV:
                if CIGAR.X in [c[0] for c in read.cigar]:
                    return True
    return False


def load_evidence(evidence_bpp: Evidence):
    """
    open the associated bam file and read and store the evidence
    does some preliminary read-quality filtering
    """

    def cache_if_true(read):
        if read.is_unmapped or read.mate_is_unmapped:
            return True
        elif any(
            [
                evidence_bpp.config['validate.filter_secondary_alignments'] and read.is_secondary,
                read.mapping_quality < evidence_bpp.min_mapping_quality,
            ]
        ):
            return False
        elif set([x[0] for x in read.cigar]) & {CIGAR.S, CIGAR.H}:
            return True
        elif not read.is_proper_pair:
            if any(
                [
                    _read.orientation_supports_type(read, e)
                    for e in evidence_bpp.putative_event_types()
                ]
            ):
                return True
            elif evidence_bpp.compatible_type and _read.orientation_supports_type(
                read, evidence_bpp.compatible_type
            ):
                return True
        elif not evidence_bpp.interchromosomal and not evidence_bpp.opposing_strands:
            min_frag_est = (
                abs(read.reference_start - read.next_reference_start) - evidence_bpp.read_length
            )
            max_frag_est = min_frag_est + 3 * evidence_bpp.read_length
            if (
                min_frag_est < evidence_bpp.min_expected_fragment_size
                or max_frag_est > evidence_bpp.max_expected_fragment_size
            ):
                return True

        return False

    def filter_if_true(read):
        if not cache_if_true(read):
            if any(
                [
                    evidence_bpp.config['validate.filter_secondary_alignments']
                    and read.is_secondary,
                    read.mapping_quality < evidence_bpp.min_mapping_quality,
                ]
            ):
                return True
            elif not evidence_bpp.interchromosomal and set([x[0] for x in read.cigar]) & {
                CIGAR.I,
                CIGAR.D,
            }:
                return False
            return True
        return False

    flanking_pairs = set()  # collect putative pairs
    half_mapped_partners1 = set()
    half_mapped_partners2 = set()

    for read in evidence_bpp.bam_cache.fetch_from_bins(
        '{0}'.format(evidence_bpp.break1.chr),
        evidence_bpp.outer_window1[0],
        evidence_bpp.outer_window1[1],
        read_limit=evidence_bpp.config['validate.fetch_reads_limit'],
        sample_bins=evidence_bpp.config['validate.fetch_reads_bins'],
        min_bin_size=evidence_bpp.config['validate.fetch_min_bin_size'],
        cache=True,
        cache_if=cache_if_true,
        filter_if=filter_if_true,
    ):
        if read.mapping_quality < evidence_bpp.min_mapping_quality:
            continue
        evidence_bpp.counts[0] += 1
        if read.is_unmapped:
            continue
        if not collect_split_read(evidence_bpp, read, True):
            collect_spanning_read(evidence_bpp, read)
        if read.mate_is_unmapped:
            half_mapped_partners1.add(read)
        elif (
            any(
                [
                    _read.orientation_supports_type(read, et)
                    for et in evidence_bpp.putative_event_types()
                ]
            )
            and (read.reference_id != read.next_reference_id) == evidence_bpp.interchromosomal
        ):
            flanking_pairs.add(read)

    for read in evidence_bpp.bam_cache.fetch_from_bins(
        '{0}'.format(evidence_bpp.break2.chr),
        evidence_bpp.outer_window2[0],
        evidence_bpp.outer_window2[1],
        read_limit=evidence_bpp.config['validate.fetch_reads_limit'],
        sample_bins=evidence_bpp.config['validate.fetch_reads_bins'],
        min_bin_size=evidence_bpp.config['validate.fetch_min_bin_size'],
        cache=True,
        cache_if=cache_if_true,
        filter_if=filter_if_true,
    ):
        if read.mapping_quality < evidence_bpp.min_mapping_quality:
            continue

        evidence_bpp.counts[1] += 1

        if read.is_unmapped:
            continue
        if not collect_split_read(evidence_bpp, read, False):
            collect_spanning_read(evidence_bpp, read)
        if read.mate_is_unmapped:
            half_mapped_partners2.add(read)
        elif (
            any(
                [
                    _read.orientation_supports_type(read, et)
                    for et in evidence_bpp.putative_event_types()
                ]
            )
            and (read.reference_id != read.next_reference_id) == evidence_bpp.interchromosomal
        ):
            flanking_pairs.add(read)
    for flanking_read in sorted(flanking_pairs, key=lambda x: (x.query_name, x.reference_start)):
        # try and get the mate from the cache
        try:
            mates = evidence_bpp.bam_cache.get_mate(flanking_read, allow_file_access=False)
            for mate in mates:
                if mate.is_unmapped:
                    logger.debug(f'ignoring unmapped mate {mate.query_name}')
                    continue
                collect_flanking_pair(evidence_bpp, flanking_read, mate)
        except KeyError:
            pass

    if evidence_bpp.compatible_window1:
        compatible_type = SVTYPE.DUP
        if SVTYPE.DUP in evidence_bpp.putative_event_types():
            compatible_type = SVTYPE.INS

        compt_flanking = set()
        for read in evidence_bpp.bam_cache.fetch_from_bins(
            '{0}'.format(evidence_bpp.break1.chr),
            evidence_bpp.compatible_window1[0],
            evidence_bpp.compatible_window1[1],
            read_limit=evidence_bpp.config['validate.fetch_reads_limit'],
            sample_bins=evidence_bpp.config['validate.fetch_reads_bins'],
            min_bin_size=evidence_bpp.config['validate.fetch_min_bin_size'],
            cache=True,
            cache_if=cache_if_true,
            filter_if=filter_if_true,
        ):
            if _read.orientation_supports_type(read, compatible_type):
                compt_flanking.add(read)

        for read in evidence_bpp.bam_cache.fetch_from_bins(
            '{0}'.format(evidence_bpp.break2.chr),
            evidence_bpp.compatible_window2[0],
            evidence_bpp.compatible_window2[1],
            read_limit=evidence_bpp.config['validate.fetch_reads_limit'],
            sample_bins=evidence_bpp.config['validate.fetch_reads_bins'],
            min_bin_size=evidence_bpp.config['validate.fetch_min_bin_size'],
            cache=True,
            cache_if=cache_if_true,
            filter_if=filter_if_true,
        ):
            if _read.orientation_supports_type(read, compatible_type):
                compt_flanking.add(read)

        for flanking_read in compt_flanking:
            # try and get the mate from the cache
            try:
                mates = evidence_bpp.bam_cache.get_mate(flanking_read, allow_file_access=False)
                for mate in mates:
                    if mate.is_unmapped:
                        logger.debug(f'ignoring unmapped mate {mate.query_name}')
                        continue
                    try:
                        collect_compatible_flanking_pair(
                            evidence_bpp, flanking_read, mate, compatible_type
                        )
                    except ValueError:
                        pass
            except KeyError:
                pass

    # now collect the half mapped reads
    logger.info(
        f'collected {len(half_mapped_partners1 | half_mapped_partners2)} putative half mapped reads',
    )
    mates_found = 0
    for read in half_mapped_partners1 | half_mapped_partners2:
        # try and get the mate from the cache
        try:
            mates = evidence_bpp.bam_cache.get_mate(read, allow_file_access=False)
            mates_found += 1
            for mate in mates:
                collect_half_mapped(evidence_bpp, read, mate)
        except KeyError:
            pass
    logger.info(f'{mates_found} half-mapped mates found')
