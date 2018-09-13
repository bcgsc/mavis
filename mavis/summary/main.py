from functools import partial
import os
import re
import time

import tab

from .constants import DEFAULTS, HOMOPOLYMER_MIN_LENGTH
from .summary import annotate_dgv, filter_by_annotations, filter_by_call_method, filter_by_evidence, get_pairing_state, group_by_distance
from ..constants import CALL_METHOD, COLUMNS, PROTOCOL, SVTYPE
from ..pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from ..util import generate_complete_stamp, LOG, output_tabbed_file, read_inputs, soft_cast


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
    filter_trans_homopolymers=DEFAULTS.filter_trans_homopolymers,
    filter_min_linking_split_reads=DEFAULTS.filter_min_linking_split_reads,
    filter_min_complexity=DEFAULTS.filter_min_complexity,
    flanking_call_distance=PAIRING_DEFAULTS.flanking_call_distance,
    split_call_distance=PAIRING_DEFAULTS.split_call_distance,
    contig_call_distance=PAIRING_DEFAULTS.contig_call_distance,
    spanning_call_distance=PAIRING_DEFAULTS.spanning_call_distance,
    start_time=int(time.time()),
    **kwargs
):
    annotations.load()
    if dgv_annotation:
        dgv_annotation.load()
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
            COLUMNS.product_id,
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
            COLUMNS.assumed_untemplated,
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
    for chr, genes in annotations.content.items():
        for gene in genes:
            for t in gene.transcripts:
                reference_transcripts[t.name] = t
                if t.is_best_transcript:
                    best_transcripts[t.name] = t

    filtered_pairs = []
    temp = []  # store the bpps while we filter out

    for bpp in bpps:
        # filter by synonymous and RNA homopolymers
        if filter_protein_synon and bpp.protein_synon:
            bpp.data[COLUMNS.filter_comment] = 'synonymous protein'
            filtered_pairs.append(bpp)
            continue
        elif filter_cdna_synon and bpp.cdna_synon:
            bpp.data[COLUMNS.filter_comment] = 'synonymous cdna'
            filtered_pairs.append(bpp)
            continue
        elif all([
            filter_trans_homopolymers,
            bpp.protocol == PROTOCOL.TRANS,
            bpp.data.get(COLUMNS.repeat_count, None),
            bpp.event_type in [SVTYPE.DUP, SVTYPE.INS, SVTYPE.DEL]
        ]):
            # a transcriptome event in a repeat region
            match = re.match(r'^(-?\d+)-(-?\d+)$', str(bpp.data[COLUMNS.net_size]))
            if match:
                netsize_min = abs(int(match.group(1)))
                netsize_max = abs(int(match.group(2)))

                if all([
                    int(bpp.repeat_count) + 1 >= HOMOPOLYMER_MIN_LENGTH,  # repeat count is 1 less than the length of the repeat
                    netsize_min == netsize_max and netsize_min == 1,
                    PROTOCOL.GENOME not in bpp.data.get(COLUMNS.pairing, '')
                ]):
                    bpp.data[COLUMNS.filter_comment] = 'homopolymer filter'
                    filtered_pairs.append(bpp)
                    continue
        # filter based on the sequence call complexity
        sc = str(bpp.data.get(COLUMNS.call_sequence_complexity, 'none')).lower()
        if sc != 'none' and float(sc) < filter_min_complexity:
            bpp.data[COLUMNS.filter_comment] = 'low complexity'
            filtered_pairs.append(bpp)
            continue
        temp.append(bpp)
    bpps = temp  # reassign the filtered result

    # filter based on minimum evidence levels
    bpps, filtered = filter_by_evidence(
        bpps, filter_min_remapped_reads=filter_min_remapped_reads,
        filter_min_spanning_reads=filter_min_spanning_reads,
        filter_min_flanking_reads=filter_min_flanking_reads,
        filter_min_split_reads=filter_min_split_reads,
        filter_min_linking_split_reads=filter_min_linking_split_reads
    )
    for pair in filtered:
        pair.data[COLUMNS.filter_comment] = 'low evidence'
        filtered_pairs.append(pair)

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
            result, removed = filter_by_call_method(bpp_set)
            collapsed.extend(result)
            for bpp in removed:
                bpp.data[COLUMNS.filter_comment] = 'collapsed into another call'
                filtered_pairs.append(bpp)
        bpps_by_library[library] = collapsed

    # collapse similar annotations for breakpoints with the same call position
    for library in bpps_by_library:
        uncollapsed = dict()
        for bpp in bpps_by_library[library]:
            uncollapsed.setdefault(bpp, []).append(bpp)

        collapsed = []
        for bpp_set in uncollapsed.values():
            result, removed = filter_by_annotations(bpp_set, best_transcripts)
            collapsed.extend(result)
            for bpp in removed:
                bpp.data[COLUMNS.filter_comment] = 'collapsed into another call'
                filtered_pairs.append(bpp)
        bpps_by_library[library] = collapsed

    # group close calls with identical annotations
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
            grouped, removed = group_by_distance(bpp_set, distances)
            collapsed.extend(grouped)
            for bpp in removed:
                bpp.data[COLUMNS.filter_comment] = 'collapsed into another call'
                filtered_pairs.append(bpp)
        bpps_by_library[library] = collapsed

    # TODO: give an evidence score to the events based on call method and evidence levels
    # TODO: report the pairings so that germline and somatic etc can be determined properly
    output_columns = {
        COLUMNS.annotation_id,
        COLUMNS.product_id,
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
        COLUMNS.fusion_splicing_pattern,
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
        COLUMNS.supplementary_call,
        COLUMNS.protein_synon,
        COLUMNS.cdna_synon,
        COLUMNS.net_size,
        COLUMNS.assumed_untemplated,
        'dgv'}

    rows = []
    for lib in bpps_by_library:
        LOG('annotating dgv for', lib)
        if dgv_annotation:
            annotate_dgv(bpps_by_library[lib], dgv_annotation.content, distance=10)  # TODO make distance a parameter
        LOG('adding pairing states for', lib)
        for row in bpps_by_library[lib]:
            # in case no pairing was done, add default (applicable to single library summaries)
            row.data.setdefault(COLUMNS.inferred_pairing, '')
            row.data.setdefault(COLUMNS.pairing, '')
            row.data.setdefault(COLUMNS.library, lib)
            # filter pairing ids based on what is still kept?
            paired_libraries = set()
            for product_id in row.pairing.split(';'):
                for lib in bpps_by_library:
                    if product_id.startswith(lib):
                        paired_libraries.add(lib)
            inferred_paired_libraries = set()
            for product_id in row.inferred_pairing.split(';'):
                for lib in bpps_by_library:
                    if product_id.startswith(lib):
                        inferred_paired_libraries.add(lib)
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

            rows.append(row.flatten())
    fname = os.path.join(
        output,
        'mavis_summary_all_{}.tab'.format('_'.join(sorted(list(libraries.keys()))))
    )
    output_tabbed_file(rows, fname, header=output_columns)
    LOG('wrote {} structural variants to {}'.format(len(rows), fname))
    output_tabbed_file(filtered_pairs, os.path.join(output, 'filtered_pairs.tab'))
    # output by library non-synon protein-product
    for lib in bpps_by_library:
        filename = os.path.join(output, 'mavis_summary_{}_non-synonymous_coding_variants.tab'.format(lib))
        lib_rows = []
        for row in rows:
            if all([
                not row.get(COLUMNS.protein_synon, ''),
                not row.get(COLUMNS.cdna_synon, ''),
                str(row.get(COLUMNS.fusion_cdna_coding_start, None)) != 'None',
                row[COLUMNS.library] == lib,
                str(row.get(COLUMNS.supplementary_call, False)) != 'True'
            ]):
                lib_rows.append(row)
        output_tabbed_file(lib_rows, filename, header=output_columns)
