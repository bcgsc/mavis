from .annotate.variant import determine_prime
from .interval import Interval
from .constants import STRAND, PRIME, CALL_METHOD, COLUMNS, ORIENT, PROTOCOL
from .error import NotSpecifiedError
from .breakpoint import Breakpoint


def predict_transcriptome_breakpoint(breakpoint, transcript):
    """
    for a given genomic breakpoint and the target transcript. Predicts the possible transcriptomic
    breakpoints that would be expected based on the splicing model for abrogated splice sites

    Args:
        breakpoint (Breakpoint): the genomic breakpoint
        transcript (usTranscript): the transcript
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
                raise NotSpecifiedError('breakpoint overlaps an exon or splice site and is not specific')
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
                intron = Interval(prev.donor_splice_site.end + 1, curr.acceptor_splice_site.start - 1)
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


def equivalent_events(ev1, ev2, TRANSCRIPTS, DISTANCES=None, SEQUENCES=None):
    if DISTANCES is None:
        DISTANCES = {CALL_METHOD.CONTIG: 0, CALL_METHOD.SPLIT: 10, CALL_METHOD.FLANK: 0}
    SEQUENCES = dict() if SEQUENCES is None else SEQUENCES
    # basic checks
    if ev1.break1.chr != ev2.break1.chr or ev1.break2.chr != ev2.break2.chr or \
            len(set([STRAND.NS, ev1.break1.strand, ev2.break1.strand])) > 2 or \
            len(set([STRAND.NS, ev1.break2.strand, ev2.break2.strand])) > 2 or \
            len(set([ORIENT.NS, ev1.break1.orient, ev2.break1.orient])) > 2 or \
            len(set([ORIENT.NS, ev1.break2.orient, ev2.break2.orient])) > 2 or \
            ev1.opposing_strands != ev2.opposing_strands:
        return False

    if ev1.data[COLUMNS.protocol] != PROTOCOL.GENOME:
        ev1, ev2 = ev2, ev1

    methods = set([
        ev1.data[COLUMNS.break1_call_method],
        ev1.data[COLUMNS.break2_call_method],
        ev2.data[COLUMNS.break1_call_method],
        ev2.data[COLUMNS.break2_call_method]
    ])
    max_distance = max([DISTANCES.get(m, 0) for m in methods])

    fusion1 = SEQUENCES.get(ev1.data[COLUMNS.fusion_sequence_fasta_id], None)
    fusion2 = SEQUENCES.get(ev2.data[COLUMNS.fusion_sequence_fasta_id], None)

    break1_match = False
    break2_match = False

    if ev1.data[COLUMNS.protocol] != ev2.data[COLUMNS.protocol]:  # mixed
        if fusion1 and fusion2:
            # compare product
            if fusion1 != fusion2:
                return False
            for col in [COLUMNS.fusion_cdna_coding_start, COLUMNS.fusion_cdna_coding_end]:
                if ev1.data[col] != ev2.data[col]:
                    return False
            return True

        # predict genome breakpoints to compare by location
        t1 = TRANSCRIPTS.get(ev1.data[COLUMNS.transcript1], None)
        if t1:
            pbreaks = predict_transcriptome_breakpoint(ev1.break1, t1)
            for b in pbreaks:
                if abs(Interval.dist(b, ev2.break1)) <= max_distance:
                    break1_match = True
                    break
        t2 = TRANSCRIPTS.get(ev1.data[COLUMNS.transcript2], None)
        if t2:
            pbreaks = predict_transcriptome_breakpoint(ev1.break2, t2)
            for b in pbreaks:
                if abs(Interval.dist(b, ev2.break2)) <= max_distance:
                    break2_match = True
                    break

    elif ev1.data[COLUMNS.event_type] != ev2.data[COLUMNS.event_type]:
        return False

    # location comparison
    if abs(Interval.dist(ev1.break1, ev2.break1)) <= max_distance:
        break1_match = True
    if abs(Interval.dist(ev1.break2, ev2.break2)) <= max_distance:
        break2_match = True
    return break1_match and break2_match
