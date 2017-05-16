from Bio import SeqIO
from ..constants import COLUMNS, STRAND, sort_columns, CALL_METHOD
from ..util import read_inputs, log
from ..pairing import equivalent_events
from ..pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .constants import DEFAULTS
import itertools

from .summary import group_events, compare_bpp_annotations, annotate_aliases


def main(
    inputs, output, annotations,
    product_sequence_files=None,
    filter_min_contig_alignment_score=DEFAULTS.filter_min_contig_alignment_score,
    filter_min_flanking_pairs=DEFAULTS.filter_min_flanking_pairs,
    filter_min_linking_split_reads=DEFAULTS.filter_min_linking_split_reads,
    filter_min_spanning_reads=DEFAULTS.filter_min_spanning_reads,
    filter_min_split_reads=DEFAULTS.filter_min_split_reads,
    flanking_call_distance=PAIRING_DEFAULTS.flanking_call_distance,
    split_call_distance=PAIRING_DEFAULTS.split_call_distance,
    contig_call_distance=PAIRING_DEFAULTS.contig_call_distance,
    spanning_call_distance=PAIRING_DEFAULTS.spanning_call_distance,
    **kwargs
):
    # pairing threshold parameters to be defined in config file
    DISTANCES = {
        CALL_METHOD.FLANK: flanking_call_distance,
        CALL_METHOD.SPLIT: split_call_distance,
        CALL_METHOD.CONTIG: contig_call_distance,
        CALL_METHOD.SPAN: spanning_call_distance
    }

    # ranking scores of the methods (more is better)
    RANKING = {
        CALL_METHOD.FLANK: 1,
        CALL_METHOD.SPLIT: 3,
        CALL_METHOD.CONTIG: 4,
        CALL_METHOD.SPAN: 2
    }

    THRESHOLDS = {
        CALL_METHOD.FLANK: filter_min_flanking_pairs,
        CALL_METHOD.SPLIT: filter_min_split_reads,
        CALL_METHOD.CONTIG: filter_min_contig_alignment_score,
        CALL_METHOD.SPAN: filter_min_spanning_reads,
    }

    SUPPORT = {
        CALL_METHOD.FLANK: (COLUMNS.flanking_pairs, COLUMNS.flanking_pairs),
        CALL_METHOD.SPLIT: (COLUMNS.break1_split_reads, COLUMNS.break2_split_reads),
        CALL_METHOD.CONTIG: (COLUMNS.contig_alignment_score, COLUMNS.contig_alignment_score),
        CALL_METHOD.SPAN: (COLUMNS.spanning_reads, COLUMNS.spanning_reads)

        # check on linking_split_reads, # contigs assembled
    }

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
        add={
            'gene1_aliases': None,
            'gene2_aliases': None,
            'summary_pairing': None},
        explicit_strand=True,
        expand_ns=False
    ))

    # load all transcripts
    reference_transcripts = dict()
    best_transcripts = dict()
    for chr, genes in annotations.items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, t)
                reference_transcripts[t.name] = t
                if t.is_best_transcript:
                    best_transcripts[t.name] = t

    # TODO: give an evidence score to the events based on call method and evidence levels
    filtered_out_bpps = []
    bpps_to_keep = dict()
    bpp_by_product_key = dict()
    product_sequences = dict()
    pairings = dict()
    product_sequence_files = set()

    for bpp in bpps:
        break1_call_method = bpp.data[COLUMNS.break1_call_method]
        break2_call_method = bpp.data[COLUMNS.break2_call_method]

        # filter low evidence
        if int(bpp.data[SUPPORT[break1_call_method][0]]) < THRESHOLDS[break1_call_method] or \
                int(bpp.data[SUPPORT[break2_call_method][1]]) < THRESHOLDS[break2_call_method]:
            filtered_out_bpps.append(bpp)
        else:
            lib = bpp.data[COLUMNS.library]
            # info needed for pairing
            if bpp.fusion_sequence_fasta_id:
                product_sequences[bpp.fusion_sequence_fasta_id] = None
            if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
                product_sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
            if lib not in pairings:
                pairings[lib] = dict()

            pos = (bpp.break1.chr,
                   bpp.break1.start,
                   bpp.break1.end,
                   bpp.break2.chr,
                   bpp.break2.start,
                   bpp.break2.end)

            bpp = annotate_aliases(bpp, reference_transcripts)

            if lib not in bpps_to_keep:
                bpps_to_keep[lib] = dict()

            if pos in bpps_to_keep[lib]:
                bpp_call_score = RANKING[bpp.data[COLUMNS.break1_call_method]] + \
                    RANKING[bpp.data[COLUMNS.break2_call_method]]
                existing_call_score = RANKING[bpps_to_keep[lib][pos].data[COLUMNS.break1_call_method]] + \
                    RANKING[bpps_to_keep[lib][pos].data[COLUMNS.break2_call_method]]

                if bpp_call_score < existing_call_score:
                    filtered_out_bpps.extend(bpp)
                elif bpp_call_score > existing_call_score:
                    filtered_out_bpps.extend(bpps_to_keep[lib][pos])
                    bpps_to_keep[lib][pos] = bpp
                else:
                    # Filter based on the annotations
                    better_bpp = compare_bpp_annotations(bpp, bpps_to_keep[lib][pos], best_transcripts)
                    bpps_to_keep[lib][pos] = better_bpp
            else:
                bpps_to_keep[lib][pos] = bpp

    # info needed for pairing
    for fname in sorted(list(product_sequence_files)):
        log('loading:', fname)
        try:
            with open(fname, 'rU') as fh:
                temp = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
                for fid, fseq in temp.items():
                    if fid in product_sequences and product_sequences[fid] is not None and \
                            product_sequences[fid] != fseq:
                        raise AssertionError('sequence identifiers are not unique', fid, fseq, product_sequences[fid])
                    product_sequences[fid] = fseq
        except IOError as err:
            log('failed for open input file', err)

    # ensure that all sequences have been found
    for seqid, seq in product_sequences.items():
        if seq is None:
            raise KeyError('failed to find sequence for the product', seqid, seq)
        product_sequences[seqid] = str(seq.seq)

    for lib in bpps_to_keep.keys():
        log("pairing", len(bpps_to_keep[lib].values())*len(bpps_to_keep[lib].values()), ' events for lib ', lib)
        for bpp in bpps_to_keep[lib].values():
            pairings[bpp.data[COLUMNS.library]][bpp.data[COLUMNS.product_id]] = set()
            bpp_by_product_key[bpp.data[COLUMNS.product_id]] = bpp

        for bpp1, bpp2 in itertools.product(bpps_to_keep[lib].values(), bpps_to_keep[lib].values()):
            if equivalent_events(
                bpp1,
                bpp2,
                DISTANCES=DISTANCES,
                reference_transcripts=reference_transcripts,
                product_sequences=product_sequences
            ):
                pairings[lib][bpp1.data[COLUMNS.product_id]].add(bpp2.data[COLUMNS.product_id])
                pairings[lib][bpp2.data[COLUMNS.product_id]].add(bpp1.data[COLUMNS.product_id])

    log('filtering based on transcript')
    # todo: add actual pairing information (i.e somatic, germline)
    bpp_to_keep = set()
    for lib in pairings:
        log(len(pairings[lib].items()), ' pairings found for lib ', lib)
        for product_key, paired_product_keys in pairings[lib].items():
            bpp = bpp_by_product_key[product_key]

        # filter any matches where genes match but transcripts do not
            filtered = []
            for paired_product_key in paired_product_keys:
                paired_bpp = bpp_by_product_key[paired_product_key]

                if bpp.data[COLUMNS.gene1] and bpp.data[COLUMNS.gene1] == paired_bpp.data[COLUMNS.gene1]:
                    if bpp.data[COLUMNS.transcript1] != paired_bpp.data[COLUMNS.transcript1]:
                        continue
                if bpp.data[COLUMNS.gene2] and bpp.data[COLUMNS.gene2] == paired_bpp.data[COLUMNS.gene2]:
                    if bpp.data[COLUMNS.transcript2] != paired_bpp.data[COLUMNS.transcript2]:
                        continue
                filtered.append(paired_product_key)
                bpp = group_events(bpp, paired_bpp)
            bpp.data['summary_pairing'] = ';'.join(sorted(filtered))
            bpp_by_product_key[product_key] = bpp
            bpp_to_keep.add(bpp)

    full_output_columns = [COLUMNS.break1_chromosome,
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
#                           COLUMNS.gene1_direction,
                           COLUMNS.gene2,
#                           COLUMNS.gene2_direction,
                           COLUMNS.gene_product_type,
#                           COLUMNS.genes_encompassed,
#                           COLUMNS.genes_overlapping_break1,
#                           COLUMNS.genes_overlapping_break2,
#                           COLUMNS.genes_proximal_to_break1,
#                           COLUMNS.genes_proximal_to_break2,
                           COLUMNS.library,
                           COLUMNS.protocol,
                           COLUMNS.transcript1,
                           COLUMNS.transcript2,
                           COLUMNS.untemplated_seq,
                           COLUMNS.tools,
#                           COLUMNS.gene1_aliases,
#                           COLUMNS.gene2_aliases,

                           # For debugging
                           COLUMNS.pairing,
                           COLUMNS.flanking_pairs,
                           COLUMNS.break1_split_reads,
                           COLUMNS.break2_split_reads,
                           COLUMNS.contig_alignment_score,
                           COLUMNS.spanning_reads,
                           'summary_pairing']

    rows = []
    for row in list(bpp_to_keep):
        try:
            row = row.flatten()
        except AttributeError:
            pass
        rows.append(row)

    header = sort_columns(full_output_columns)
    with open(output, 'w') as fh:
        log('writing', output)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')
    log("Wrote {} gene fusion events to {}".format(len(rows), output))


if __name__ == '__main__':
    main()
