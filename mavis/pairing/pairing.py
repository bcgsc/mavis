from .constants import DEFAULTS, PAIRING_DISTANCES

from ..annotate.variant import determine_prime
from ..breakpoint import Breakpoint
from ..constants import CALL_METHOD, COLUMNS, ORIENT, PRIME, PROTOCOL, STRAND
from ..error import NotSpecifiedError
from ..interval import Interval
from ..util import DEVNULL


def product_key(bpp):
    """
    unique id for the product row
    """
    return '_'.join([str(v) for v in [
        bpp.library,
        bpp.protocol,
        bpp.annotation_id,
        bpp.fusion_splicing_pattern,
        bpp.fusion_cdna_coding_start,
        bpp.fusion_cdna_coding_end
    ]])


def predict_transcriptome_breakpoint(breakpoint, transcript):
    """
    for a given genomic breakpoint and the target transcript. Predicts the possible transcriptomic
    breakpoints that would be expected based on the splicing model for abrogated splice sites

    Args:
        breakpoint (Breakpoint): the genomic breakpoint
        transcript (PreTranscript): the transcript

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


def comparison_distance(event1, event2, input_distances=None):
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

    max_distance = comparison_distance(event1, event2, distances)

    if not _equivalent_events(event1, event2):
        return False
    seqlen = sum([
        len(event1.untemplated_seq) if event1.untemplated_seq else 0,
        len(event2.untemplated_seq) if event2.untemplated_seq else 0
    ])
    max_distance += seqlen
    # location comparison
    if any([
        abs(Interval.dist(event1.break1, event2.break1)) > max_distance,
        abs(Interval.dist(event1.break2, event2.break2)) > max_distance,
        event1.data[COLUMNS.event_type] != event2.data[COLUMNS.event_type]
    ]):
        return False
    return True


def pair_by_distance(calls, distances, log=DEVNULL, against_self=False):
    """
    for a set of input calls, pair by distance
    """
    distance_pairings = {}
    break1_sorted = sorted(calls, key=lambda b: b.break1.start)
    break2_sorted = sorted(calls, key=lambda b: b.break2.start)
    lowest_resolution = max([len(b.break1) for b in calls] + [len(b.break2) for b in calls] + [1])
    max_distance = max(distances.values())
    max_useq = max([len(c.untemplated_seq) if c.untemplated_seq else 0 for c in calls] + [0])
    max_distance += max_useq * 2
    log('lowest_resolution', lowest_resolution, 'max_distance', max_distance, 'possible comparisons', len(break1_sorted) * len(break1_sorted), time_stamp=False)

    comparisons = 0
    for i in range(0, len(break1_sorted)):
        current = break1_sorted[i]
        distance_pairings.setdefault(product_key(current), set())
        for j in range(i + 1, len(break1_sorted)):
            other = break1_sorted[j]

            dist = abs(Interval.dist(current.break1, other.break1))
            if dist > max_distance + lowest_resolution:
                break
            comparisons += 1
            if not against_self and current.library == other.library and current.protocol == other.protocol:
                continue  # do not pair within a single library
            if equivalent(current, other, distances=distances):
                distance_pairings[product_key(current)].add(product_key(other))
                distance_pairings.setdefault(product_key(other), set()).add(product_key(current))
        current = break2_sorted[i]
        for j in range(i + 1, len(break2_sorted)):
            other = break2_sorted[j]
            dist = abs(Interval.dist(current.break2, other.break2))
            if dist > max_distance + lowest_resolution:
                break
            comparisons += 1
            if not against_self and current.library == other.library and current.protocol == other.protocol:
                continue  # do not pair within a single library
            if equivalent(current, other, distances=distances):
                distance_pairings.setdefault(product_key(current), set()).add(product_key(other))
                distance_pairings.setdefault(product_key(other), set()).add(product_key(current))
    log('computed {} comparisons'.format(comparisons), time_stamp=False)
    return distance_pairings


def inferred_equivalent(event1, event2, reference_transcripts, distances=None):
    """
    comparison of events using product prediction and breakpoint prediction
    """
    # basic checks
    if not _equivalent_events(event1, event2):
        return False

    if event1.data[COLUMNS.protocol] != PROTOCOL.GENOME:
        event1, event2 = event2, event1

    max_distance = comparison_distance(event1, event2, distances)

    if event1.data[COLUMNS.fusion_sequence_fasta_id] and event2.data[COLUMNS.fusion_sequence_fasta_id]:
        if event1.fusion_sequence_fasta_id != event2.fusion_sequence_fasta_id:
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
