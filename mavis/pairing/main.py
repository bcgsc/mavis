import os
import itertools
from Bio import SeqIO
from ..constants import PROTOCOL, COLUMNS, CALL_METHOD, SVTYPE, STRAND
from ..annotate.constants import SPLICE_TYPE
from .constants import DEFAULTS
from ..util import read_inputs, output_tabbed_file, log, generate_complete_stamp
from . import equivalent_events


def main(
    inputs, output, annotations, product_sequence_files=None,
    flanking_call_distance=DEFAULTS.flanking_call_distance,
    split_call_distance=DEFAULTS.split_call_distance,
    contig_call_distance=DEFAULTS.contig_call_distance,
    spanning_call_distance=DEFAULTS.spanning_call_distance,
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        flanking_call_distance (int): pairing distance for pairing with an event called by :term:`flanking read pair`
        split_call_distance (int): pairing distance for pairing with an event called by :term:`split read`
        contig_call_distance (int): pairing distance for pairing with an event called by contig or :term:`spanning read`
    """
    product_sequence_files = set() if product_sequence_files is None else set(product_sequence_files)
    # load the file
    DISTANCES = {
        CALL_METHOD.FLANK: flanking_call_distance,
        CALL_METHOD.SPLIT: split_call_distance,
        CALL_METHOD.CONTIG: contig_call_distance,
        CALL_METHOD.SPAN: spanning_call_distance
    }

    bpps = []
    bpps.extend(read_inputs(
        inputs,
        require=[
            COLUMNS.annotation_id,
            COLUMNS.library,
            COLUMNS.fusion_cdna_coding_start,
            COLUMNS.fusion_cdna_coding_end,
            COLUMNS.fusion_sequence_fasta_id,
            COLUMNS.fusion_sequence_fasta_file
        ],
        in_={
            COLUMNS.protocol: PROTOCOL,
            COLUMNS.event_type: SVTYPE,
            COLUMNS.break1_call_method: CALL_METHOD,
            COLUMNS.break2_call_method: CALL_METHOD,
            COLUMNS.fusion_splicing_pattern: SPLICE_TYPE.values() + [None, 'None']
        },
        add={
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.fusion_sequence_fasta_file: None,
            COLUMNS.fusion_splicing_pattern: None
        },
        explicit_strand=False,
        expand_ns=False
    ))
    log('read {} breakpoint pairs'.format(len(bpps)))
    libraries = set()

    product_sequences = dict()
    # get all sequence ids and all files to read from
    for bpp in bpps:
        if bpp.fusion_sequence_fasta_file:
            product_sequence_files.add(bpp.fusion_sequence_fasta_file)
        if bpp.fusion_sequence_fasta_id:
            product_sequences[bpp.fusion_sequence_fasta_id] = None
        libraries.add(bpp.library)

    # load sequences from all files detected
    log('detected', len(libraries), 'libraries')
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

    # load all transcripts
    reference_transcripts = dict()
    for chr, genes in annotations.items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, t)
                reference_transcripts[t.name] = t

    # map the calls by library and ensure there are no name/key conflicts
    calls_by_lib = dict()
    bpp_by_product_key = dict()
    pairings = dict()
    categories = set()

    for bpp in bpps:
        product_key = (
            bpp.library,
            bpp.protocol,
            bpp.annotation_id,
            bpp.fusion_splicing_pattern,
            bpp.fusion_cdna_coding_start,
            bpp.fusion_cdna_coding_end
        )
        category = (bpp.break1.chr, bpp.break2.chr, bpp.opposing_strands)
        categories.add(category)
        bpp.data[COLUMNS.product_id] = product_key
        calls_by_lib.setdefault(bpp.library, {})
        calls_by_lib[bpp.library].setdefault(category, set())
        calls_by_lib[bpp.library][category].add(product_key)

        pairings[product_key] = set()
        if product_key in bpp_by_product_key:
            raise KeyError('duplicate bpp is not unique within lib', bpp.library, product_key, bpp, bpp.data)

        bpp_by_product_key[product_key] = bpp
    total_comparisons = 0
    # pairwise comparison of breakpoints between all libraries
    for category in sorted(list(categories)):
        for lib, other_lib in itertools.combinations(calls_by_lib.keys(), 2):
            # create combinations from other libraries in the same category
            pairs = calls_by_lib[lib].get(category, set())
            other_pairs = calls_by_lib[other_lib].get(category, set())
            c = len(pairs) * len(other_pairs)
            total_comparisons += c
            # for each two libraries pair all calls
            if c > 10000:
                log(c, 'comparison(s) between', lib, 'and', other_lib, 'for', category)

            for product_key1, product_key2 in itertools.product(pairs, other_pairs):
                if product_key1[:-3] == product_key2[:-3]:
                    continue
                if equivalent_events(
                    bpp_by_product_key[product_key1],
                    bpp_by_product_key[product_key2],
                    DISTANCES=DISTANCES,
                    reference_transcripts=reference_transcripts,
                    product_sequences=product_sequences
                ):
                    pairings[product_key1].add(product_key2)
                    pairings[product_key2].add(product_key1)
    log('checked', total_comparisons, 'total comparisons')
    for product_key, paired_product_keys in pairings.items():
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
        bpp.data[COLUMNS.pairing] = ';'.join(['_'.join([str(v) for v in key]) for key in sorted(filtered)])

    fname = os.path.join(
        output,
        'mavis_paired_{}.tab'.format('_'.join(sorted(list(libraries))))
    )
    output_tabbed_file(bpps, fname)
    generate_complete_stamp(output, log)
