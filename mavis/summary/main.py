from Bio import SeqIO
from ..constants import COLUMNS, sort_columns, CALL_METHOD
from ..util import read_inputs, log
from ..pairing import equivalent_events
from ..pairing.constants import DEFAULTS as PAIRING_DEFAULTS
from .constants import DEFAULTS
import itertools

from .summary import filter_by_evidence, group_events, filter_by_annotations, annotate_aliases,  filter_by_call_method

def main(
    inputs, output, annotations,
    product_sequence_files=None,
    filter_min_remapped_reads=DEFAULTS.filter_min_remapped_reads,
    filter_min_spanning_reads=DEFAULTS.filter_min_spanning_reads,
    filter_min_flanking_reads=DEFAULTS.filter_min_flanking_reads,
    filter_min_flanking_only_reads=DEFAULTS.filter_min_flanking_only_reads,
    filter_min_split_reads=DEFAULTS.filter_min_split_reads,
    filter_min_linking_split_reads=DEFAULTS.filter_min_linking_split_reads,
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

    bpps = []
    bpps.extend(read_inputs(
        inputs,
        require=[COLUMNS.break1_chromosome,
                 COLUMNS.break1_orientation,
                 COLUMNS.break1_position_end,
                 COLUMNS.break1_position_start,
                 COLUMNS.break2_chromosome,
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
                 ],
        add={
            COLUMNS.gene1_aliases: None,
            COLUMNS.gene2_aliases: None,
            'summary_pairing': None},
        explicit_strand=True,
        expand_ns=False,
        cast={COLUMNS.break1_split_reads: int,
              COLUMNS.break2_split_reads: int,
              COLUMNS.contig_remapped_reads: lambda x: 0 if x == 'None' or x is None else int(x),
              COLUMNS.spanning_reads: int,
              COLUMNS.break1_split_reads_forced: int,
              COLUMNS.break2_split_reads_forced: int,
              COLUMNS.flanking_pairs: int,
              COLUMNS.linking_split_reads: int}
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

    bpps = filter_by_evidence(bpps, filter_min_remapped_reads=filter_min_remapped_reads,
                              filter_min_spanning_reads=filter_min_spanning_reads,
                              filter_min_flanking_reads=filter_min_flanking_reads,
                              filter_min_flanking_only_reads=filter_min_flanking_only_reads,
                              filter_min_split_reads=filter_min_split_reads,
                              filter_min_linking_split_reads=filter_min_linking_split_reads)

    bpps_to_keep = dict()
    bpp_by_product_key = dict()
    product_sequences = dict()
    pairings = dict()
    product_sequence_files = set()

    for bpp in bpps:
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
               bpp.break2.end,
               bpp.event_type)
#               bpp.break1.strand,
#               bpp.break2.strand)

        bpp = annotate_aliases(bpp, reference_transcripts)

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
            pairings[bpp.data[COLUMNS.library]][bpp.data[COLUMNS.product_id]] = set()
            bpp_by_product_key[bpp.data[COLUMNS.product_id]] = bpp
            calls_by_lib.setdefault(lib, {})
            calls_by_lib[lib].setdefault(category, set())
            calls_by_lib[lib][category].add(bpp)
        for category in sorted(list(categories)):

            c = len(calls_by_lib[lib][category])*len(calls_by_lib[lib][category])
            total_comparisons += c
            if c > 10000:
                log(c,'comparisons for lib ', lib,'for', category)

            for bpp1, bpp2 in itertools.product(calls_by_lib[lib][category], calls_by_lib[lib][category]):
                if bpp1 is not None and equivalent_events(
                    bpp1,
                    bpp2,
                    DISTANCES=DISTANCES,
                    reference_transcripts=reference_transcripts,
                    product_sequences=product_sequences
                    ):
                    pairings[lib][bpp1.data[COLUMNS.product_id]].add(bpp2.data[COLUMNS.product_id])
                    pairings[lib][bpp2.data[COLUMNS.product_id]].add(bpp1.data[COLUMNS.product_id])
    log('checked', total_comparisons, 'total comparisons')

    log('filtering pairings based on transcript')
    bpp_to_keep = dict()
    for lib in pairings:
        if lib not in bpp_to_keep:
            bpp_to_keep[lib] = set()
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

                # take only the one with the higher ranking call method
                try:
                    bpp = filter_by_call_method(bpp, paired_bpp)
                except AssertionError:
                    # call method is the same, so just group them
                    bpp = group_events(bpp, paired_bpp)

            bpp.data['summary_pairing'] = ';'.join(sorted(filtered))
            bpp_by_product_key[product_key] = bpp
            bpp_to_keep[lib].add(bpp)

    # TODO: give an evidence score to the events based on call method and evidence levels
    # TODO: report the pairings so that germline and somatic etc can be determined properly
    output_columns = [
        COLUMNS.pairing,
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

        # For debugging
        COLUMNS.break1_call_method,
        COLUMNS.break2_call_method,
        COLUMNS.flanking_pairs,
        COLUMNS.break1_split_reads,
        COLUMNS.break2_split_reads,
        COLUMNS.contig_alignment_score,
        COLUMNS.spanning_reads,
        COLUMNS.flanking_pairs_read_names,
        COLUMNS.contig_remapped_reads,
        COLUMNS.contig_remapped_read_names,
        'summary_pairing']

    rows = []
    for lib in bpp_to_keep:
        for row in list(bpp_to_keep[lib]):
            try:
                row = row.flatten()
            except AttributeError:
                pass
            rows.append(row)

    header = sort_columns(output_columns)
    with open(output, 'w') as fh:
        log('writing', output)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in sorted(rows, key= lambda k: (k[COLUMNS.break1_chromosome],
                                                int(k[COLUMNS.break1_position_start]),
                                                k[COLUMNS.break2_chromosome],
                                                int(k[COLUMNS.break2_position_start]))):
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')
    log("Wrote {} gene fusion events to {}".format(len(rows), output))


if __name__ == '__main__':
    main()
