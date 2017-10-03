from .constants import DEFAULTS, PAIRING_DISTANCES

from ..annotate.variant import determine_prime
from ..breakpoint import Breakpoint
from ..constants import CALL_METHOD, COLUMNS, ORIENT, PRIME, PROTOCOL, STRAND
from ..error import NotSpecifiedError
from ..interval import Interval


def predict_transcriptome_breakpoint(breakpoint, transcript):
    """
    for a given genomic breakpoint and the target transcript. Predicts the possible transcriptomic
    breakpoints that would be expected based on the splicing model for abrogated splice sites

    Args:
        breakpoint (Breakpoint): the genomic breakpoint
        transcript (UsTranscript): the transcript

    see :ref:`theory - pairing similar events <theory-pairing-similar-events>`
    """
    prime = determine_prime(transcript, breakpoint)
    exons = transcript.exons[:]
    if not Interval.overlaps(breakpoint, transcript):
        raise AssertionError('breakpoint does not overlap the transcript', breakpoint, transcript)
    if transcript.get_strand() == STRAND.NEG:
        exons.reverse()

    tbreaks = []

    for i, curr in enumerate(exons):
        temp = curr.acceptor_splice_site | curr.donor_splice_site

        if Interval.overlaps(breakpoint, temp):  # overlaps a splice site or exon
            if len(breakpoint) > 1:
                raise NotSpecifiedError(
                    'breakpoint overlaps an exon or splice site and is not specific (i.e. has '
                    'an interval greater than 1)')
            elif prime == PRIME.FIVE:
                if i > 0:
                    prev = exons[i - 1]
                    tbreaks.append(
                        Breakpoint(
                            breakpoint.chr,
                            prev.donor,
                            strand=breakpoint.strand,
                            orient=breakpoint.orient
                        ))
            else:  # prime == PRIME.THREE
                if i + 1 < len(exons):
                    nexxt = exons[i + 1]
                    tbreaks.append(
                        Breakpoint(
                            breakpoint.chr,
                            nexxt.acceptor,
                            strand=breakpoint.strand,
                            orient=breakpoint.orient
                        ))
            tbreaks.append(
                Breakpoint(
                    breakpoint.chr, breakpoint.start, breakpoint.end,
                    strand=breakpoint.strand,
                    orient=breakpoint.orient
                ))
        elif i > 0:  # look at the previous intron
            prev = exons[i - 1]
            try:
                intron_start, intron_end = sorted([prev.donor_splice_site.end, curr.acceptor_splice_site.start])
                intron = Interval(intron_start + 1, intron_end - 1)
                if Interval.overlaps(breakpoint, intron):
                    if prime == PRIME.FIVE:
                        tbreaks.append(
                            Breakpoint(
                                breakpoint.chr,
                                prev.donor,
                                strand=breakpoint.strand,
                                orient=breakpoint.orient
                            ))
                    else:
                        tbreaks.append(
                            Breakpoint(
                                breakpoint.chr,
                                curr.acceptor,
                                strand=breakpoint.strand,
                                orient=breakpoint.orient
                            ))
            except AttributeError:  # for introns that are smaller than this ignore (covered by exon check)
                pass

    if not tbreaks:
        raise AssertionError('could not predict breakpoint')
    return sorted(tbreaks)


def _equivalent_events(event1, event2):
    # basic checks
    if any([
        event1.break1.chr != event2.break1.chr,
        event1.break2.chr != event2.break2.chr,
        len(set([STRAND.NS, event1.break1.strand, event2.break1.strand])) > 2,
        len(set([STRAND.NS, event1.break2.strand, event2.break2.strand])) > 2,
        len(set([ORIENT.NS, event1.break1.orient, event2.break1.orient])) > 2,
        len(set([ORIENT.NS, event1.break2.orient, event2.break2.orient])) > 2,
        event1.opposing_strands != event2.opposing_strands
    ]):
        return False
    return True


def compairson_distance(event1, event2, input_distances=None):
    distances = {}
    distances.update(PAIRING_DISTANCES.items())
    if input_distances is not None:
        distances.update(input_distances)
    max_distance = max(
        distances[event1.data.get(COLUMNS.call_method, CALL_METHOD.CONTIG)],
        distances[event2.data.get(COLUMNS.call_method, CALL_METHOD.CONTIG)])
    return max_distance


def equivalent(event1, event2, distances=None):
    """
    compares two events by breakpoint position to see if they are equivalent
    """

    max_distance = compairson_distance(event1, event2, distances)

    if not _equivalent_events(event1, event2):
        return False
    # location comparison
    if any([
        not _equivalent_events(event1, event2),
        abs(Interval.dist(event1.break1, event2.break1)) > max_distance,
        abs(Interval.dist(event1.break2, event2.break2)) > max_distance,
        event1.data[COLUMNS.event_type] != event2.data[COLUMNS.event_type]
    ]):
        return False
    return True


def inferred_equivalent(event1, event2, reference_transcripts, distances=None, product_sequences=None):
    """
    comparison of events using product prediction and breakpoint prediction
    """
    product_sequences = dict() if product_sequences is None else product_sequences
    # basic checks
    if not _equivalent_events(event1, event2):
        return False

    if event1.data[COLUMNS.protocol] != PROTOCOL.GENOME:
        event1, event2 = event2, event1

    max_distance = compairson_distance(event1, event2, distances)

    if event1.data[COLUMNS.fusion_sequence_fasta_id] and event2.data[COLUMNS.fusion_sequence_fasta_id]:
        fusion1 = product_sequences[event1.data[COLUMNS.fusion_sequence_fasta_id]]
        fusion2 = product_sequences[event2.data[COLUMNS.fusion_sequence_fasta_id]]

        if fusion1 != fusion2:
            return False
        for col in [COLUMNS.fusion_cdna_coding_start, COLUMNS.fusion_cdna_coding_end]:
            if event1.data[col] != event2.data[col]:
                return False
        return True

    break1_match = False
    break2_match = False

    if event1.data[COLUMNS.protocol] != event2.data[COLUMNS.protocol]:  # mixed
        # predict genome breakpoints to compare by location
        transcript1 = reference_transcripts.get(event1.data[COLUMNS.transcript1], None)
        if transcript1:
            try:
                pbreaks = predict_transcriptome_breakpoint(event1.break1, transcript1)
                for breakpoint in pbreaks:
                    if abs(Interval.dist(breakpoint, event2.break1)) <= max_distance:
                        break1_match = True
                        break
            except NotSpecifiedError:
                pass
        transcript2 = reference_transcripts.get(event1.data[COLUMNS.transcript2], None)
        if transcript2:
            try:
                pbreaks = predict_transcriptome_breakpoint(event1.break2, transcript2)
                for breakpoint in pbreaks:
                    if abs(Interval.dist(breakpoint, event2.break2)) <= max_distance:
                        break2_match = True
                        break
            except NotSpecifiedError:
                pass
    elif event1.data[COLUMNS.event_type] != event2.data[COLUMNS.event_type]:
        return False

    if abs(Interval.dist(event1.break1, event2.break1)) <= max_distance:
        break1_match = True
    if abs(Interval.dist(event1.break2, event2.break2)) <= max_distance:
        break2_match = True

    return break1_match and break2_match
