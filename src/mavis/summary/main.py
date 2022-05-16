import os
import re
import time
from typing import Dict, List, Tuple

import pandas as pd

from ..annotate.file_io import ReferenceFile
from ..breakpoint import BreakpointPair
from ..constants import CALL_METHOD, COLUMNS, PROTOCOL, SPLICE_TYPE, SVTYPE
from ..util import generate_complete_stamp, logger, output_tabbed_file, read_inputs
from .constants import HOMOPOLYMER_MIN_LENGTH
from .summary import (
    annotate_dgv,
    filter_by_annotations,
    filter_by_call_method,
    filter_by_evidence,
    get_pairing_state,
    group_by_distance,
)


def main(inputs: List[str], output: str, config: Dict, start_time=int(time.time())):
    annotations = ReferenceFile.load_from_config(config, 'annotations', eager_load=True)
    dgv_annotation = ReferenceFile.load_from_config(config, 'dgv_annotation')
    if not dgv_annotation.is_empty():
        dgv_annotation.load()
    # pairing threshold parameters to be defined in config file
    distances = {
        CALL_METHOD.FLANK: config['pairing.flanking_call_distance'],
        CALL_METHOD.SPLIT: config['pairing.split_call_distance'],
        CALL_METHOD.CONTIG: config['pairing.contig_call_distance'],
        CALL_METHOD.SPAN: config['pairing.spanning_call_distance'],
    }

    bpps = []
    bpps.extend(
        read_inputs(
            inputs,
            required_columns=[
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
                COLUMNS.disease_status,
            ],
            add_default={
                COLUMNS.call_method: CALL_METHOD.INPUT,
            },
            apply={
                COLUMNS.event_type: lambda x: SVTYPE.enforce(x),
                COLUMNS.fusion_splicing_pattern: lambda x: SPLICE_TYPE.enforce(x)
                if not pd.isnull(x)
                else x,
            },
            expand_strand=False,
            expand_orient=False,
            expand_svtype=False,
        )
    )
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
        if config['summary.filter_protein_synon'] and bpp.protein_synon:
            bpp.data[COLUMNS.filter_comment] = 'synonymous protein'
            filtered_pairs.append(bpp)
            continue
        elif config['summary.filter_cdna_synon'] and bpp.cdna_synon:
            bpp.data[COLUMNS.filter_comment] = 'synonymous cdna'
            filtered_pairs.append(bpp)
            continue
        elif all(
            [
                config['summary.filter_trans_homopolymers'],
                bpp.protocol == PROTOCOL.TRANS,
                bpp.data.get(COLUMNS.repeat_count, None),
                bpp.event_type in [SVTYPE.DUP, SVTYPE.INS, SVTYPE.DEL],
            ]
        ):
            # a transcriptome event in a repeat region
            match = re.match(r'^(-?\d+)-(-?\d+)$', str(bpp.data[COLUMNS.net_size]))
            if match:
                netsize_min = abs(int(match.group(1)))
                netsize_max = abs(int(match.group(2)))

                if all(
                    [
                        int(bpp.repeat_count) + 1
                        >= HOMOPOLYMER_MIN_LENGTH,  # repeat count is 1 less than the length of the repeat
                        netsize_min == netsize_max and netsize_min == 1,
                        PROTOCOL.GENOME not in bpp.data.get(COLUMNS.pairing, ''),
                    ]
                ):
                    bpp.data[COLUMNS.filter_comment] = 'homopolymer filter'
                    filtered_pairs.append(bpp)
                    continue
        # filter based on the sequence call complexity
        sc = str(bpp.data.get(COLUMNS.call_sequence_complexity, 'none')).lower()
        if sc != 'none' and float(sc) < config['summary.filter_min_complexity']:
            bpp.data[COLUMNS.filter_comment] = 'low complexity'
            filtered_pairs.append(bpp)
            continue
        temp.append(bpp)
    bpps = temp  # reassign the filtered result

    # filter based on minimum evidence levels
    bpps, filtered = filter_by_evidence(
        bpps,
        filter_min_remapped_reads=config['summary.filter_min_remapped_reads'],
        filter_min_spanning_reads=config['summary.filter_min_spanning_reads'],
        filter_min_flanking_reads=config['summary.filter_min_flanking_reads'],
        filter_min_split_reads=config['summary.filter_min_split_reads'],
        filter_min_linking_split_reads=config['summary.filter_min_linking_split_reads'],
    )
    for pair in filtered:
        pair.data[COLUMNS.filter_comment] = 'low evidence'
        filtered_pairs.append(pair)

    bpps_by_library: Dict[str, List[BreakpointPair]] = {}  # split the input pairs by library
    libraries = {}
    for bpp in bpps:
        bpps_by_library.setdefault(bpp.library, []).append(bpp)
        libraries[bpp.library] = (bpp.protocol, bpp.disease_status)

    # collapse identical calls with different call methods
    for library in bpps_by_library:
        uncollapsed: Dict[Tuple, List[BreakpointPair]] = dict()
        for bpp in bpps_by_library[library]:
            group: Tuple[BreakpointPair, str, str, str, str, int, int] = (
                bpp,
                bpp.data.get(COLUMNS.transcript1),
                bpp.data.get(COLUMNS.transcript2),
                bpp.fusion_sequence_fasta_id,
                bpp.fusion_splicing_pattern,
                bpp.fusion_cdna_coding_start,
                bpp.fusion_cdna_coding_end,
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
            uncollapsed.setdefault(
                (
                    bpp.event_type,
                    bpp.break1.chr,
                    bpp.break2.chr,
                    bpp.break1.orient,
                    bpp.break2.orient,
                    bpp.opposing_strands,
                    bpp.break1.strand,
                    bpp.break2.strand,
                    bpp.data.get(COLUMNS.transcript1) if bpp.data.get(COLUMNS.gene1) else None,
                    bpp.data.get(COLUMNS.transcript2) if bpp.data.get(COLUMNS.gene2) else None,
                    bpp.fusion_sequence_fasta_id,  # id is a hash of the sequence
                    bpp.fusion_cdna_coding_start,
                    bpp.fusion_cdna_coding_end,
                ),
                [],
            ).append(bpp)

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
        'dgv',
    }

    rows = []
    for lib in bpps_by_library:
        logger.info(f'annotating dgv for {lib}')
        if not dgv_annotation.is_empty():
            annotate_dgv(
                bpps_by_library[lib], dgv_annotation.content, distance=10
            )  # TODO make distance a parameter
        logger.info(f'adding pairing states for {lib}')
        for row in bpps_by_library[lib]:
            # in case no pairing was done, add default (applicable to single library summaries)
            row.data.setdefault(COLUMNS.inferred_pairing, '')
            row.data.setdefault(COLUMNS.pairing, '')
            row.data.setdefault(COLUMNS.library, lib)
            # filter pairing ids based on what is still kept?
            paired_libraries = set()
            for product_id in (row.pairing or '').split(';'):
                for lib in bpps_by_library:
                    if product_id.startswith(lib):
                        paired_libraries.add(lib)
            inferred_paired_libraries = set()
            for product_id in (row.inferred_pairing or '').split(';'):
                for lib in bpps_by_library:
                    if product_id.startswith(lib):
                        inferred_paired_libraries.add(lib)
            for other_lib, (other_protocol, other_disease_state) in libraries.items():
                column_name = '{}_{}_{}'.format(other_lib, other_disease_state, other_protocol)
                if other_lib != row.library:
                    pairing_state = get_pairing_state(
                        *libraries[row.library],
                        other_protocol=other_protocol,
                        other_disease_state=other_disease_state,
                        is_matched=other_lib in paired_libraries,
                        inferred_is_matched=other_lib in inferred_paired_libraries,
                    )
                else:
                    pairing_state = 'Not Applicable'
                row.data[column_name] = pairing_state
                output_columns.add(column_name)

            rows.append(row.flatten())
    fname = os.path.join(
        output, 'mavis_summary_all_{}.tab'.format('_'.join(sorted(list(libraries.keys()))))
    )
    output_tabbed_file(rows, fname, header=output_columns)
    logger.info(f'wrote {len(rows)} structural variants to {fname}')
    output_tabbed_file(filtered_pairs, os.path.join(output, 'filtered_pairs.tab'))
    # output by library non-synon protein-product
    for lib in bpps_by_library:
        filename = os.path.join(
            output, 'mavis_summary_{}_non-synonymous_coding_variants.tab'.format(lib)
        )
        lib_rows = []
        for row in rows:
            if all(
                [
                    not row.get(COLUMNS.protein_synon, ''),
                    not row.get(COLUMNS.cdna_synon, ''),
                    str(row.get(COLUMNS.fusion_cdna_coding_start, None)) != 'None',
                    row[COLUMNS.library] == lib,
                    str(row.get(COLUMNS.supplementary_call, False)) != 'True',
                ]
            ):
                lib_rows.append(row)
        output_tabbed_file(lib_rows, filename, header=output_columns)
    generate_complete_stamp(output, start_time=start_time)
