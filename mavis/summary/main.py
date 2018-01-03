from functools import partial
import itertools
import os
import time

import tab

from .constants import DEFAULTS
from .summary import annotate_dgv, filter_by_annotations, filter_by_call_method, filter_by_evidence, get_pairing_state, group_by_distance
from ..constants import CALL_METHOD, COLUMNS
from ..pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from ..util import generate_complete_stamp, log, output_tabbed_file, read_inputs, soft_cast


def soft_cast_null(value):
    try:
        return tab.cast_null(value)
    except TypeError:
        return value


def main(
    inputs, output, annotations,
    dgv_annotation=None,
    filter_cdna_synon=DEFAULTS.filter_cdna_synon,
    filter_protein_synon=DEFAULTS.filter_protein_synon,
    filter_min_remapped_reads=DEFAULTS.filter_min_remapped_reads,
    filter_min_spanning_reads=DEFAULTS.filter_min_spanning_reads,
    filter_min_flanking_reads=DEFAULTS.filter_min_flanking_reads,
    filter_min_split_reads=DEFAULTS.filter_min_split_reads,
    filter_min_linking_split_reads=DEFAULTS.filter_min_linking_split_reads,
    flanking_call_distance=PAIRING_DEFAULTS.flanking_call_distance,
    split_call_distance=PAIRING_DEFAULTS.split_call_distance,
    contig_call_distance=PAIRING_DEFAULTS.contig_call_distance,
    spanning_call_distance=PAIRING_DEFAULTS.spanning_call_distance,
    start_time=int(time.time()),
    **kwargs
):
    # pairing threshold parameters to be defined in config file
    distances = {
        CALL_METHOD.FLANK: flanking_call_distance,
        CALL_METHOD.SPLIT: split_call_distance,
        CALL_METHOD.CONTIG: contig_call_distance,
        CALL_METHOD.SPAN: spanning_call_distance
    }

    bpps = []
    bpps.extend(read_inputs(
        inputs,
        require=[
            COLUMNS.event_type,
            COLUMNS.fusion_cdna_coding_end,
            COLUMNS.fusion_cdna_coding_start,
            COLUMNS.fusion_splicing_pattern,
            COLUMNS.fusion_mapped_domains,
            COLUMNS.gene1,
            COLUMNS.gene1_direction,
            COLUMNS.gene2,
            COLUMNS.gene2_direction,
            COLUMNS.gene_product_type,
            COLUMNS.genes_encompassed,
            COLUMNS.library,
            COLUMNS.protocol,
            COLUMNS.transcript1,
            COLUMNS.transcript2,
            COLUMNS.untemplated_seq,
            COLUMNS.tools,
            COLUMNS.exon_last_5prime,
            COLUMNS.exon_first_3prime,
            COLUMNS.disease_status
        ],
        add_default={**{k: None for k in [
            COLUMNS.contig_remapped_reads,
            COLUMNS.contig_seq,
            COLUMNS.break1_split_reads,
            COLUMNS.break1_split_reads_forced,
            COLUMNS.break2_split_reads,
            COLUMNS.break2_split_reads_forced,
            COLUMNS.linking_split_reads,
            COLUMNS.flanking_pairs,
            COLUMNS.contigs_aligned,
            COLUMNS.contigs_assembled,
            COLUMNS.contig_alignment_score,
            COLUMNS.contig_remap_score,
            COLUMNS.spanning_reads,
            COLUMNS.annotation_figure,
            COLUMNS.gene1_aliases,
            COLUMNS.gene2_aliases,
            COLUMNS.protein_synon,
            COLUMNS.cdna_synon,
            COLUMNS.net_size,
            COLUMNS.tracking_id,
            'dgv',
            'summary_pairing']
        }, COLUMNS.call_method: CALL_METHOD.INPUT},
        expand_strand=False, expand_orient=False, expand_svtype=False,
        cast={
            COLUMNS.break1_split_reads: partial(soft_cast, cast_type=int),
            COLUMNS.break2_split_reads: partial(soft_cast, cast_type=int),
            COLUMNS.contig_remapped_reads: partial(soft_cast, cast_type=int),
            COLUMNS.spanning_reads: partial(soft_cast, cast_type=int),
            COLUMNS.break1_split_reads_forced: partial(soft_cast, cast_type=int),
            COLUMNS.break2_split_reads_forced: partial(soft_cast, cast_type=int),
            COLUMNS.flanking_pairs: partial(soft_cast, cast_type=int),
            COLUMNS.linking_split_reads: partial(soft_cast, cast_type=int),
            COLUMNS.protein_synon: soft_cast_null,
            COLUMNS.cdna_synon: soft_cast_null
        }
    ))
    # load all transcripts
    reference_transcripts = dict()
    best_transcripts = dict()
    for chr, genes in annotations.items():
        for gene in genes:
            for t in gene.transcripts:
                reference_transcripts[t.name] = t
                if t.is_best_transcript:
                    best_transcripts[t.name] = t

    # filter by synonymous
    if filter_cdna_synon or filter_protein_synon:
        temp = []
        for bpp in bpps:
            if filter_protein_synon and bpp.protein_synon:
                continue
            elif filter_cdna_synon and bpp.cdna_synon:
                continue
            temp.append(bpp)
        bpps = temp

    # filter based on minimum evidence levels
    bpps, _ = filter_by_evidence(
        bpps, filter_min_remapped_reads=filter_min_remapped_reads,
        filter_min_spanning_reads=filter_min_spanning_reads,
        filter_min_flanking_reads=filter_min_flanking_reads,
        filter_min_split_reads=filter_min_split_reads,
        filter_min_linking_split_reads=filter_min_linking_split_reads
    )

    bpps_by_library = {}  # split the input pairs by library
    libraries = {}
    for bpp in bpps:
        bpps_by_library.setdefault(bpp.library, []).append(bpp)
        libraries[bpp.library] = (bpp.protocol, bpp.disease_status)

    # collapse identical calls with different call methods
    for library in bpps_by_library:
        uncollapsed = dict()
        for bpp in bpps_by_library[library]:
            group = (
                bpp,
                bpp.transcript1,
                bpp.transcript2,
                bpp.fusion_sequence_fasta_id,
                bpp.fusion_splicing_pattern,
                bpp.fusion_cdna_coding_start,
                bpp.fusion_cdna_coding_end
            )
            uncollapsed.setdefault(group, []).append(bpp)
        collapsed = []
        for bpp_set in uncollapsed.values():
            collapsed.extend(filter_by_call_method(bpp_set))
        bpps_by_library[library] = collapsed

    # collapse similar annotations for breakpoints with the same call position
    for library in bpps_by_library:
        uncollapsed = dict()
        for bpp in bpps_by_library[library]:
            uncollapsed.setdefault(bpp, []).append(bpp)

        collapsed = []
        for bpp_set in uncollapsed.values():
            collapsed.extend(filter_by_annotations(bpp_set, best_transcripts))
        bpps_by_library[library] = collapsed

    # group close split read calls with identical annotations
    for library in bpps_by_library:
        uncollapsed = dict()
        for bpp in bpps_by_library[library]:
            uncollapsed.setdefault((
                bpp.event_type,
                bpp.break1.chr, bpp.break2.chr,
                bpp.break1.orient, bpp.break2.orient,
                bpp.opposing_strands,
                bpp.break1.strand, bpp.break2.strand,
                bpp.transcript1 if bpp.gene1 else None,
                bpp.transcript2 if bpp.gene2 else None,
                bpp.fusion_sequence_fasta_id,  # id is a hash of the sequence
                bpp.fusion_cdna_coding_start,
                bpp.fusion_cdna_coding_end
            ), []).append(bpp)

        collapsed = []
        for bpp_set in uncollapsed.values():
            collapsed.extend([b for b in bpp_set if b.call_method != CALL_METHOD.SPLIT])
            collapsed.extend(group_by_distance([b for b in bpp_set if b.call_method == CALL_METHOD.SPLIT], distances))
        bpps_by_library[library] = collapsed

    # TODO: give an evidence score to the events based on call method and evidence levels
    # TODO: report the pairings so that germline and somatic etc can be determined properly
    output_columns = {
        COLUMNS.annotation_id,
        COLUMNS.break1_chromosome,
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
        COLUMNS.fusion_protein_hgvs,
        COLUMNS.fusion_mapped_domains,
        COLUMNS.gene1,
        COLUMNS.gene1_direction,
        COLUMNS.gene2,
        COLUMNS.gene2_direction,
        COLUMNS.gene_product_type,
        COLUMNS.genes_encompassed,
        COLUMNS.library,
        COLUMNS.protocol,
        COLUMNS.transcript1,
        COLUMNS.transcript2,
        COLUMNS.untemplated_seq,
        COLUMNS.tools,
        COLUMNS.break1_strand,
        COLUMNS.break2_strand,
        COLUMNS.gene1_aliases,
        COLUMNS.gene2_aliases,
        COLUMNS.annotation_figure,
        COLUMNS.exon_last_5prime,
        COLUMNS.exon_first_3prime,

        # For debugging
        COLUMNS.call_method,
        COLUMNS.flanking_pairs,
        COLUMNS.break1_split_reads,
        COLUMNS.break2_split_reads,
        COLUMNS.linking_split_reads,
        COLUMNS.contig_alignment_score,
        COLUMNS.spanning_reads,
        COLUMNS.contig_remapped_reads,
        COLUMNS.tracking_id,
        COLUMNS.net_size,
        'dgv'}

    rows = []
    for lib in bpps_by_library:
        log('annotating dgv for', lib)
        if dgv_annotation:
            annotate_dgv(bpps_by_library[lib], dgv_annotation, distance=10)  # TODO make distance a parameter
        log('adding pairing states for', lib)
        for row in bpps_by_library[lib]:
            # in case no pairing was done, add default (applicable to single library summaries)
            row.data.setdefault(COLUMNS.inferred_pairing, '')
            row.data.setdefault(COLUMNS.pairing, '')
            row.data.setdefault(COLUMNS.library, lib)
            # filter pairing ids based on what is still kept?
            paired_libraries = set([p.split('_')[0] for p in row.pairing.split(';')])
            inferred_paired_libraries = set([p.split('_')[0] for p in row.inferred_pairing.split(';')])
            for other_lib, (other_protocol, other_disease_state) in libraries.items():
                column_name = '{}_{}_{}'.format(other_lib, other_disease_state, other_protocol)
                if other_lib != row.library:
                    pairing_state = get_pairing_state(
                        *libraries[row.library],
                        other_protocol=other_protocol, other_disease_state=other_disease_state,
                        is_matched=other_lib in paired_libraries,
                        inferred_is_matched=other_lib in inferred_paired_libraries)
                else:
                    pairing_state = 'Not Applicable'
                row.data[column_name] = pairing_state
                output_columns.add(column_name)

            rows.append(row)
    fname = os.path.join(
        output,
        'mavis_summary_{}.tab'.format('_'.join(sorted(list(libraries.keys()))))
    )
    rows = sorted(rows, key=lambda bpp: (bpp.break1, bpp.break2))
    output_tabbed_file(rows, fname, header=output_columns)
    log('Wrote {} gene fusion events to {}'.format(len(rows), fname))
    generate_complete_stamp(output, log, start_time=start_time)
