# /projects/trans_scratch/validations/workspace/cchoo/overlays/COLO-829/two_tools_to_validated/pairing/mavis_paired_A36971_A36973.tab

import TSV
from mavis.constants import COLUMNS, STRAND, sort_columns, CALL_METHOD
from mavis.annotate import load_reference_genes
from mavis.util import read_inputs, log, output_tabbed_file
import os
import argparse
from pprint import pprint

# set up a config file


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
    """
    annotation_info = [COLUMNS.gene1,
                       COLUMNS.gene2,
                       COLUMNS.transcript1,
                       COLUMNS.transcript2,
                       COLUMNS.fusion_cdna_coding_start,
                       COLUMNS.fusion_cdna_coding_end,
                       COLUMNS.break1_strand,
                       COLUMNS.break2_strand]

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
                pass
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
            pass
        pass

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


def main():
    # threshold parameters to be defined in config file
    min_contig_alignment_score = 5
    min_flanking_pairs = 4  # flanking evidence required to call an event by flanking
    min_linking_split_reads = 2
    min_spanning_reads = 3
    min_split_reads = 3  # includes forced split reads

    RANKING = {
        CALL_METHOD.FLANK: 1,
        CALL_METHOD.SPLIT: 3,
        CALL_METHOD.CONTIG: 4,
        CALL_METHOD.SPAN: 2
    }

    THRESHOLDS = {
        CALL_METHOD.FLANK: min_flanking_pairs,
        CALL_METHOD.SPLIT: min_split_reads,
        CALL_METHOD.CONTIG: min_contig_alignment_score,
        CALL_METHOD.SPAN: min_spanning_reads,
    }

    SUPPORT = {
        CALL_METHOD.FLANK: (COLUMNS.flanking_pairs, COLUMNS.flanking_pairs),
        CALL_METHOD.SPLIT: (COLUMNS.break1_split_reads, COLUMNS.break2_split_reads),
        CALL_METHOD.CONTIG: (COLUMNS.contig_alignment_score, COLUMNS.contig_alignment_score),
        CALL_METHOD.SPAN: (COLUMNS.spanning_reads, COLUMNS.spanning_reads)
    }

    parser = argparse.ArgumentParser(
        description='generate MAVIS summary file',
        add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-o', '--output', help='path to the output file', required=True)
    required.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('--annotations', help='path to the annotation file')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print('error: input file {0} does not exist'.format(args.input))
        sys.exit()
    print('reading:', args.input)
    inputs = [args.input]

    bpps = []
    bpps.extend(read_inputs(
        inputs,
        require=[COLUMNS.break1_chromosome,
                 COLUMNS.break1_homologous_seq,
                 COLUMNS.break1_orientation,
                 COLUMNS.break1_position_end,
                 COLUMNS.break1_position_start,
                 COLUMNS.break2_chromosome,
                 COLUMNS.break2_homologous_seq,
                 COLUMNS.break2_orientation,
                 COLUMNS.break2_position_end,
                 COLUMNS.break2_position_start,
                 COLUMNS.contig_seq,
                 COLUMNS.event_type,
                 COLUMNS.fusion_cdna_coding_end,
                 COLUMNS.fusion_cdna_coding_start,
                 COLUMNS.fusion_mapped_domains,
                 COLUMNS.gene1,
                 COLUMNS.gene1_direction,
                 COLUMNS.gene2,
                 COLUMNS.gene2_direction,
                 COLUMNS.gene_product_type,
                 COLUMNS.genes_encompassed,
                 COLUMNS.genes_overlapping_break1,
                 COLUMNS.genes_overlapping_break2,
                 COLUMNS.genes_proximal_to_break1,
                 COLUMNS.genes_proximal_to_break2,
                 COLUMNS.library,
                 COLUMNS.protocol,
                 COLUMNS.transcript1,
                 COLUMNS.transcript2,
                 COLUMNS.untemplated_seq,
                 COLUMNS.tools,
                 # evidence_columns
                 COLUMNS.break1_call_method,
                 COLUMNS.break1_split_reads,
                 COLUMNS.break1_split_reads_forced,
                 COLUMNS.break2_call_method,
                 COLUMNS.break2_split_reads,
                 COLUMNS.break2_split_reads_forced,
                 COLUMNS.linking_split_reads,
                 COLUMNS.flanking_pairs,
                 COLUMNS.contigs_aligned,
                 COLUMNS.contigs_assembled,
                 COLUMNS.contig_alignment_score,
                 COLUMNS.contig_remap_score,
                 COLUMNS.raw_break1_split_reads,
                 COLUMNS.raw_break2_split_reads,
                 COLUMNS.raw_break1_half_mapped_reads,
                 COLUMNS.raw_break2_half_mapped_reads,
                 # id
                 COLUMNS.product_id,
                 COLUMNS.annotation_id,
                 COLUMNS.validation_id
                 ],
        explicit_strand=True,
        expand_ns=False
    ))

    log('read {} breakpoint pairs'.format(len(bpps)))
    args.annotations = load_reference_genes(args.annotations, best_transcripts_only=True)
    # load all transcripts
    reference_transcripts = dict()
    for chr, genes in args.annotations.items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, t)
                reference_transcripts[t.name] = t

    # filter low evidence and give a evidence score?
    filtered_out_bpps = []
    bpp_by_validate_id = dict()
    bpp_by_pos = []
    bpps_to_keep = {}
    for bpp in bpps:
        event_id = '_'.join([bpp.data[COLUMNS.product_id], bpp.data[COLUMNS.annotation_id]])
        break1_call_method = bpp.data[COLUMNS.break1_call_method]
        break2_call_method = bpp.data[COLUMNS.break2_call_method]

        if int(bpp.data[SUPPORT[break1_call_method][0]]) < THRESHOLDS[break1_call_method] or \
                int(bpp.data[SUPPORT[break2_call_method][1]]) < THRESHOLDS[break2_call_method]:
            filtered_out_bpps.append(bpp)
        else:
            pos = (bpp.break1.chr,
                   bpp.break1.start,
                   bpp.break1.end,
                   bpp.break2.chr,
                   bpp.break2.start,
                   bpp.break2.end)

            # todo: add filtering here based on positions
            if pos in bpps_to_keep:
                bpp_call_score = RANKING[bpp.data[COLUMNS.break1_call_method]] + \
                    RANKING[bpp.data[COLUMNS.break2_call_method]]
                existing_call_score = RANKING[bpps_to_keep[pos].data[COLUMNS.break1_call_method]] + \
                    RANKING[bpps_to_keep[pos].data[COLUMNS.break2_call_method]]

                if bpp_call_score < existing_call_score:
                    filtered_out_bpps.extend(bpp)
                elif bpp_call_score > existing_call_score:
                    bpps_to_keep[pos] = bpp
                else:
                    # aggregate contigs and untemplated if present?
                    pass
                    bpps_to_keep[pos] = compare_bpp_annotations(bpp, bpps_to_keep[pos], args.annotations)
            else:
                bpps_to_keep[pos] = bpp
            # Also aggregate the contig and annotation information?
            # if bpp.data[COLUMNS.validation_id] not in bpp_by_validate_id and bpp in bpp_by_pos:
            #     pass
            # else:
            #     bpp_by_pos.append(bpp)
            #     bpp_by_validate_id[bpp.data[COLUMNS.validation_id]] = bpp
            #     bpp_by_product_key[event_id] = bpp

    # Filter based on the annotations, need to check
    fname = args.output
#    output_tabbed_file(bpps_to_keep.values(), fname)

    columns = [COLUMNS.break1_chromosome,
               COLUMNS.break1_homologous_seq,
               COLUMNS.break1_orientation,
               COLUMNS.break1_position_end,
               COLUMNS.break1_position_start,
               COLUMNS.break2_chromosome,
               COLUMNS.break2_homologous_seq,
               COLUMNS.break2_orientation,
               COLUMNS.break2_position_end,
               COLUMNS.break2_position_start,
               COLUMNS.break1_call_method,
               COLUMNS.break2_call_method,
               COLUMNS.contig_seq,
               COLUMNS.event_type,
               COLUMNS.fusion_cdna_coding_end,
               COLUMNS.fusion_cdna_coding_start,
               COLUMNS.fusion_mapped_domains,
               COLUMNS.gene1,
               COLUMNS.gene1_direction,
               COLUMNS.gene2,
               COLUMNS.gene2_direction,
               COLUMNS.gene_product_type,
               COLUMNS.genes_encompassed,
               COLUMNS.genes_overlapping_break1,
               COLUMNS.genes_overlapping_break2,
               COLUMNS.genes_proximal_to_break1,
               COLUMNS.genes_proximal_to_break2,
               COLUMNS.library,
               COLUMNS.protocol,
               COLUMNS.transcript1,
               COLUMNS.transcript2,
               COLUMNS.untemplated_seq,
               COLUMNS.tools]
    rows = []
    for row in bpps_to_keep.values():
        try:
            row = row.flatten()
        except AttributeError:
            pass
        rows.append(row)

    header = sort_columns(columns)

    with open(fname, 'w') as fh:
        log('writing', fname)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

    # elements = sort_columns(output[0].keys())
    # header = "\t".join(elements)
    # with open(args.output, 'w') as fh:
    #     fh.write(header + "\n")
    #     for event in output:
    #         line = []
    #         for element in elements:
    #             line.append(str(event[element]))
    #         fh.write("\t".join(line) + "\n")
    # print("Wrote {} gene fusion events to {}".format(len(output), args.output))


if __name__ == '__main__':
    main()
