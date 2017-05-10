from ..annotate.variant import determine_prime
from ..interval import Interval
from ..constants import STRAND, PRIME, CALL_METHOD, COLUMNS, ORIENT, PROTOCOL
from ..error import NotSpecifiedError
from ..breakpoint import Breakpoint


def predict_transcriptome_breakpoint(breakpoint, transcript):
    """
    for a given genomic breakpoint and the target transcript. Predicts the possible transcriptomic
    breakpoints that would be expected based on the splicing model for abrogated splice sites

    Args:
        breakpoint (Breakpoint): the genomic breakpoint
        transcript (usTranscript): the transcript

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
                s, t = sorted([prev.donor_splice_site.end, curr.acceptor_splice_site.start])
                intron = Interval(s + 1, t - 1)
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

    if len(tbreaks) == 0:
        raise AssertionError('could not predict breakpoint')
    return sorted(tbreaks)


def equivalent_events(ev1, ev2, reference_transcripts, DISTANCES=None, product_sequences=None):
    temp = {CALL_METHOD.CONTIG: 0, CALL_METHOD.SPLIT: 10, CALL_METHOD.FLANK: 0}
    temp.update(DISTANCES if DISTANCES else {})
    DISTANCES = temp

    product_sequences = dict() if product_sequences is None else product_sequences
    # basic checks
    if any([
        ev1.break1.chr != ev2.break1.chr,
        ev1.break2.chr != ev2.break2.chr,
        len(set([STRAND.NS, ev1.break1.strand, ev2.break1.strand])) > 2,
        len(set([STRAND.NS, ev1.break2.strand, ev2.break2.strand])) > 2,
        len(set([ORIENT.NS, ev1.break1.orient, ev2.break1.orient])) > 2,
        len(set([ORIENT.NS, ev1.break2.orient, ev2.break2.orient])) > 2,
        ev1.opposing_strands != ev2.opposing_strands
    ]):
        return False

    if ev1.data[COLUMNS.protocol] != PROTOCOL.GENOME:
        ev1, ev2 = ev2, ev1

    methods = set([
        ev1.data[COLUMNS.break1_call_method],
        ev1.data[COLUMNS.break2_call_method],
        ev2.data[COLUMNS.break1_call_method],
        ev2.data[COLUMNS.break2_call_method]
    ])
    max_distance = max([DISTANCES[m] for m in methods])

    if ev1.data[COLUMNS.fusion_sequence_fasta_id] and ev2.data[COLUMNS.fusion_sequence_fasta_id]:
        fusion1 = product_sequences[ev1.data[COLUMNS.fusion_sequence_fasta_id]]
        fusion2 = product_sequences[ev2.data[COLUMNS.fusion_sequence_fasta_id]]

        if fusion1 != fusion2:
            return False
        for col in [COLUMNS.fusion_cdna_coding_start, COLUMNS.fusion_cdna_coding_end]:
            if ev1.data[col] != ev2.data[col]:
                return False
        return True

    break1_match = False
    break2_match = False

    if ev1.data[COLUMNS.protocol] != ev2.data[COLUMNS.protocol]:  # mixed
        # predict genome breakpoints to compare by location
        t1 = reference_transcripts.get(ev1.data[COLUMNS.transcript1], None)
        if t1:
            try:
                pbreaks = predict_transcriptome_breakpoint(ev1.break1, t1)
                for b in pbreaks:
                    if abs(Interval.dist(b, ev2.break1)) <= max_distance:
                        break1_match = True
                        break
            except NotSpecifiedError:
                pass
        t2 = reference_transcripts.get(ev1.data[COLUMNS.transcript2], None)
        if t2:
            try:
                pbreaks = predict_transcriptome_breakpoint(ev1.break2, t2)
                for b in pbreaks:
                    if abs(Interval.dist(b, ev2.break2)) <= max_distance:
                        break2_match = True
                        break
            except NotSpecifiedError:
                pass
    elif ev1.data[COLUMNS.event_type] != ev2.data[COLUMNS.event_type]:
        return False
    # location comparison
    if abs(Interval.dist(ev1.break1, ev2.break1)) <= max_distance:
        break1_match = True
    if abs(Interval.dist(ev1.break2, ev2.break2)) <= max_distance:
        break2_match = True
    return break1_match and break2_match
