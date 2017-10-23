from functools import partial
import itertools
import os
import time

from Bio import SeqIO

from .constants import DEFAULTS
from .summary import annotate_dgv, filter_by_annotations, filter_by_call_method, filter_by_evidence, get_pairing_state, group_events
from ..constants import CALL_METHOD, COLUMNS
from ..pairing.pairing import equivalent
from ..pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from ..util import generate_complete_stamp, log, output_tabbed_file, read_inputs, soft_cast


def main(
    inputs, output, annotations,
    dgv_annotation=None,
    product_sequence_files=None,
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
            COLUMNS.protein_synon: partial(soft_cast, cast_type=bool)
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

    bpps, _ = filter_by_evidence(
        bpps, filter_min_remapped_reads=filter_min_remapped_reads,
        filter_min_spanning_reads=filter_min_spanning_reads,
        filter_min_flanking_reads=filter_min_flanking_reads,
        filter_min_split_reads=filter_min_split_reads,
        filter_min_linking_split_reads=filter_min_linking_split_reads)

    bpps_to_keep = dict()
    bpp_by_product_key = dict()
    product_sequences = dict()
    pairings = dict()
    product_sequence_files = set()
    libraries = dict()

    for bpp in bpps:
        lib = bpp.library
        if lib not in libraries:
            # get library info, (disease status, protocol)
            libraries[lib] = (bpp.protocol, bpp.disease_status)
        # info needed for pairing
        if bpp.fusion_sequence_fasta_id:
            product_sequences[bpp.fusion_sequence_fasta_id] = None
        if bpp.fusion_sequence_fasta_file:
            product_sequence_files.add(bpp.fusion_sequence_fasta_file)
        if lib not in pairings:
            pairings[lib] = dict()

        pos = (bpp.break1.chr,
               bpp.break1.start,
               bpp.break1.end,
               bpp.break2.chr,
               bpp.break2.start,
               bpp.break2.end,
               bpp.event_type)

        if lib not in bpps_to_keep:
            bpps_to_keep[lib] = dict()

        if pos in bpps_to_keep[lib]:
            try:
                # Filter based on call method
                bpps_to_keep[lib][pos] = filter_by_call_method(bpp, bpps_to_keep[lib][pos])
            except AssertionError:
                # Filter based on the annotations
                bpps_to_keep[lib][pos] = filter_by_annotations(bpp, bpps_to_keep[lib][pos], best_transcripts)
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

    total_comparisons = 0
    for lib in bpps_to_keep.keys():
        calls_by_lib = dict()
        categories = set()

        for bpp in bpps_to_keep[lib].values():
            category = (bpp.break1.chr, bpp.break2.chr, bpp.break1.strand, bpp.break2.strand)
            categories.add(category)
            pairings[bpp.library][bpp.product_id] = set()
            bpp_by_product_key[bpp.product_id] = bpp
            calls_by_lib.setdefault(lib, {})
            calls_by_lib[lib].setdefault(category, set())
            calls_by_lib[lib][category].add(bpp)
        for category in sorted(list(categories)):

            c = len(calls_by_lib[lib][category]) * len(calls_by_lib[lib][category])
            total_comparisons += c
            if c > 10000:
                log(c, 'comparisons for lib ', lib, 'for', category)

            for bpp1, bpp2 in itertools.product(calls_by_lib[lib][category], calls_by_lib[lib][category]):
                if bpp1 is not None and equivalent(
                    bpp1,
                    bpp2,
                    distances=distances
                ):
                    pairings[lib][bpp1.product_id].add(bpp2.product_id)
                    pairings[lib][bpp2.product_id].add(bpp1.product_id)
    log('checked', total_comparisons, 'total comparisons')

    log('filtering pairings based on transcript')
    bpp_to_keep = dict()
    annotation_ids_to_keep = []
    for lib in pairings:
        if lib not in bpp_to_keep:
            bpp_to_keep[lib] = set()
        log(len(pairings[lib].items()), 'pairings found for lib ', lib)
        for product_key, paired_product_keys in pairings[lib].items():
            bpp = bpp_by_product_key[product_key]

        # filter any matches where genes match but transcripts do not
            filtered = []
            for paired_product_key in paired_product_keys:
                paired_bpp = bpp_by_product_key[paired_product_key]

                if bpp.gene1 and bpp.gene1 == paired_bpp.gene1:
                    if bpp.transcript1 != paired_bpp.transcript1:
                        continue
                if bpp.gene2 and bpp.gene2 == paired_bpp.gene2:
                    if bpp.transcript2 != paired_bpp.transcript2:
                        continue
                filtered.append(paired_product_key)

                # take only the one with the higher ranking call method
                try:
                    bpp = filter_by_call_method(bpp, paired_bpp)
                except AssertionError:
                    # call method is the same, so just group them
                    bpp = group_events(bpp, paired_bpp)

            bpp.data['summary_pairing'] = ';'.join(sorted(filtered))
            bpp_by_product_key[product_key] = bpp
            bpp_to_keep[lib].add(bpp)
            annotation_ids_to_keep.extend(bpp.annotation_id.split(';'))

    # TODO: give an evidence score to the events based on call method and evidence levels
    # TODO: report the pairings so that germline and somatic etc can be determined properly
    output_columns = {
        COLUMNS.annotation_id,
        COLUMNS.pairing,
        COLUMNS.inferred_pairing,
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
        'summary_pairing',
        'dgv'}

    rows = []
    for lib in bpp_to_keep:
        log('annotating dgv for', lib)
        if dgv_annotation:
            annotated = annotate_dgv(list(bpp_to_keep[lib]), dgv_annotation, distance=10)  # TODO make distance a parameter
        log('adding pairing states for', lib)
        for row in annotated:
            # filter pairing ids based on what is still kept?
            paired_libraries = set([p.split('_')[0] for p in row.pairing.split(';')])
            inferred_paired_libraries = set([p.split('_')[0] for p in row.inferred_pairing.split(';')])
            for other_lib in libraries:
                other_protocol, other_disease_state = libraries[other_lib]
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
