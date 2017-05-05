from ..constants import COLUMNS, STRAND

def alphanumeric_choice(bpp1, bpp2):
    """
    """
    chosen1 = sorted([bpp1.data[COLUMNS.transcript1], bpp2.data[COLUMNS.transcript1]])[0]
    chosen2 = sorted([bpp1.data[COLUMNS.transcript2], bpp2.data[COLUMNS.transcript2]])[0]
    if bpp1.data[COLUMNS.transcript1] == chosen1 and \
            bpp1.data[COLUMNS.transcript1] != bpp2.data[COLUMNS.transcript1]:
        return bpp1
    elif bpp1.data[COLUMNS.transcript2] == chosen2:
        return bpp1
    else:
        return bpp2


def compare_bpp_annotations(bpp1, bpp2, best_transcripts):
    """
    Args:
        bpp1 (BreakPointPair): 
        bpp2 (BreakpointPair):
        best_transcripts (dict):
    """
    # By priority
    # Case 1 an event has 2 genes and transcripts and a fusion cdna (orf)
    if bpp1.data[COLUMNS.fusion_cdna_coding_start] or bpp2.data[COLUMNS.fusion_cdna_coding_start]:
        # take the one with the longest cdna length
        if bpp1.data[COLUMNS.fusion_cdna_coding_start] is None:
            return bpp2
        elif bpp2.data[COLUMNS.fusion_cdna_coding_start] is None:
            return bpp1
        else:
            bpp1_cdna_len = int(bpp1.data[COLUMNS.fusion_cdna_coding_end]) - \
                int(bpp1.data[COLUMNS.fusion_cdna_coding_start])
            bpp2_cdna_len = int(bpp2.data[COLUMNS.fusion_cdna_coding_end]) - \
                int(bpp2.data[COLUMNS.fusion_cdna_coding_start])
            return bpp1 if bpp1_cdna_len >= bpp2_cdna_len else bpp2

    # Case 2 an event has 2 genes and transcripts
    elif bpp1.data[COLUMNS.gene1] and bpp1.data[COLUMNS.gene2] or bpp2.data[COLUMNS.gene1] and bpp2.data[COLUMNS.gene2]:
        # take the one with transcripts that are in best transcript, or the highest alphanumeric name
        if bpp1.data[COLUMNS.gene1] is None or bpp1.data[COLUMNS.gene2] is None:
            return bpp2
        elif bpp2.data[COLUMNS.gene1] is None or bpp2.data[COLUMNS.gene2] is None:
            return bpp1
        else:
            bpp1_t1, bpp1_t2 = (bpp1.data[COLUMNS.transcript1], bpp1.data[COLUMNS.transcript2])
            bpp2_t1, bpp2_t2 = (bpp2.data[COLUMNS.transcript1], bpp2.data[COLUMNS.transcript2])
            # both in best transcripts
            if bpp1_t1 in best_transcripts and bpp1_t2 in best_transcripts and bpp2_t1 in best_transcripts \
                    and bpp2_t2 in best_transcripts:
                return alphanumeric_choice(bpp1, bpp2)
            elif bpp1_t1 in best_transcripts and bpp1_t2 in best_transcripts:
                return bpp1
            elif bpp2_t1 in best_transcripts and bpp2_t2 in best_transcripts:
                return bpp2
            elif bpp1_t1 in best_transcripts or bpp1_t2 in best_transcripts:
                return bpp1
            elif bpp2_t1 in best_transcripts or bpp2_t2 in best_transcripts:
                return bpp2
            else:
                return alphanumeric_choice(bpp1, bpp2)

    # Case 3 an event has 1 gene and transcript
    elif bpp1.data[COLUMNS.gene1] or bpp1.data[COLUMNS.gene2] or bpp2.data[COLUMNS.gene1] or bpp2.data[COLUMNS.gene2]:
        # take the one with transcripts that are in best transcript, or the highest alphanumeric name
        if bpp1.data[COLUMNS.gene1] is None and bpp1.data[COLUMNS.gene2] is None:
            return bpp2
        elif bpp2.data[COLUMNS.gene1] is None and bpp2.data[COLUMNS.gene2] is None:
            return bpp1
        else:
            bpp1_t1, bpp1_t2 = (bpp1.data[COLUMNS.transcript1], bpp1.data[COLUMNS.transcript2])
            bpp2_t1, bpp2_t2 = (bpp2.data[COLUMNS.transcript1], bpp2.data[COLUMNS.transcript2])

            if bpp1_t1 in best_transcripts or bpp1_t2 in best_transcripts:
                return bpp1
            elif bpp2_t1 in best_transcripts or bpp2_t2 in best_transcripts:
                return bpp2
            else:
                return alphanumeric_choice(bpp1, bpp2)

    # Case 4 both have no genes present - will keep the positive strand event
    else:
        if bpp1.break1.strand == STRAND.POS:
            return bpp1
        else:
            return bpp2


def combine_evidence(bpp_to_keep, bpp_to_add):
    # combine the untemplated sequences
    if bpp_to_add.data[COLUMNS.untemplated_seq] is not None:
        bpp_to_keep.data[COLUMNS.untemplated_seq] = bpp_to_add.data[COLUMNS.untemplated_seq] if \
            bpp_to_keep.data[COLUMNS.untemplated_seq] is None else \
            ';'.join(sorted(list(set(bpp_to_add.data[COLUMNS.untemplated_seq],
                                     bpp_to_keep.data[COLUMNS.untemplated_seq]))))

    # combine the contig sequences
    if bpp_to_add.data[COLUMNS.contig_seq] is not None:
        bpp_to_keep.data[COLUMNS.contig_seq] = bpp_to_add.data[COLUMNS.contig_seq] if \
            bpp_to_keep.data[COLUMNS.contig_seq] is None else \
            ';'.join(sorted(list(set(bpp_to_add.data[COLUMNS.contig_seq],
                                     bpp_to_keep.data[COLUMNS.contig_seq]))))


    return bpp_to_keep

def group_events(bpp1, bpp2):
    pass

def filter_events(bpp1, bpp2, THRESHOLDS):
    pass


def annotate_aliases(bpp, reference_transcripts):
    # Should add the getting the alias to annotate instead of here?
    if bpp.data[COLUMNS.transcript1] in reference_transcripts:
        bpp.data['gene1_aliases'] = ",".join(reference_transcripts[bpp.data[COLUMNS.transcript1]].gene.aliases)
    if bpp.data[COLUMNS.transcript2] in reference_transcripts:
        bpp.data['gene2_aliases'] = ",".join(reference_transcripts[bpp.data[COLUMNS.transcript2]].gene.aliases)
    return(bpp)
