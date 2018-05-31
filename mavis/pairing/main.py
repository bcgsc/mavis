import itertools
import os
import time

from .pairing import inferred_equivalent, product_key, pair_by_distance
from .constants import DEFAULTS
from ..annotate.constants import SPLICE_TYPE
from ..constants import CALL_METHOD, COLUMNS, PROTOCOL, SVTYPE
from ..util import generate_complete_stamp, LOG, output_tabbed_file, read_inputs


def main(
    inputs, output, annotations,
    flanking_call_distance=DEFAULTS.flanking_call_distance,
    split_call_distance=DEFAULTS.split_call_distance,
    contig_call_distance=DEFAULTS.contig_call_distance,
    spanning_call_distance=DEFAULTS.spanning_call_distance,
    start_time=int(time.time()),
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
    annotations.load()
    # load the file
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
            COLUMNS.annotation_id,
            COLUMNS.library,
            COLUMNS.fusion_cdna_coding_start,
            COLUMNS.fusion_cdna_coding_end,
            COLUMNS.fusion_sequence_fasta_id
        ],
        in_={
            COLUMNS.protocol: PROTOCOL.values(),
            COLUMNS.event_type: SVTYPE.values(),
            COLUMNS.fusion_splicing_pattern: SPLICE_TYPE.values() + [None, 'None']
        },
        add_default={
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.fusion_splicing_pattern: None
        },
        expand_strand=False, expand_orient=False, expand_svtype=False
    ))
    LOG('read {} breakpoint pairs'.format(len(bpps)))

    # load all transcripts
    reference_transcripts = dict()
    for genes in annotations.content.values():
        for gene in genes:
            for unspliced_t in gene.transcripts:
                if unspliced_t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, unspliced_t)
                reference_transcripts[unspliced_t.name] = unspliced_t

    # map the calls by library and ensure there are no name/key conflicts
    calls_by_cat = dict()
    calls_by_ann = dict()
    bpp_by_product_key = dict()
    libraries = set()

    # initialize the pairing mappings
    for bpp in bpps:
        libraries.add(bpp.library)
        category = (bpp.break1.chr, bpp.break2.chr, bpp.opposing_strands, bpp.event_type)
        bpp.data[COLUMNS.product_id] = product_key(bpp)
        calls_by_cat.setdefault(category, []).append(bpp)
        if bpp.gene1 or bpp.gene2:
            calls_by_ann.setdefault((bpp.transcript1, bpp.transcript2), []).append(bpp)
        bpp.data[COLUMNS.pairing] = ''
        bpp.data[COLUMNS.inferred_pairing] = ''

        if product_key(bpp) in bpp_by_product_key:
            raise KeyError('duplicate bpp is not unique within lib', product_key(bpp))
        bpp_by_product_key[product_key(bpp)] = bpp

    distance_pairings = {}
    product_pairings = {}
    LOG('computing distance based pairings')
    # pairwise comparison of breakpoints between all libraries
    for set_num, (category, calls) in enumerate(sorted(calls_by_cat.items(), key=lambda x: (len(x[1]), x[0]), reverse=True)):
        LOG('comparing set {} of {} with {} items'.format(set_num + 1, len(calls_by_cat), len(calls)))
        for node, adj_list in pair_by_distance(calls, distances, against_self=False).items():
            distance_pairings.setdefault(node, set()).update(adj_list)

    LOG('computing inferred (by product) pairings')
    for calls in calls_by_ann.values():
        calls_by_lib = {}
        for call in calls:
            calls_by_lib.setdefault(call.library, []).append(call)

        for lib, other_lib in itertools.combinations(calls_by_lib.keys(), 2):
            # create combinations from other libraries in the same category
            pairs = calls_by_lib[lib]
            other_pairs = calls_by_lib[other_lib]

            for current, other in itertools.product(pairs, other_pairs):
                if inferred_equivalent(
                    current,
                    other,
                    distances=distances,
                    reference_transcripts=reference_transcripts
                ):
                    product_pairings.setdefault(product_key(current), set()).add(product_key(other))
                    product_pairings.setdefault(product_key(other), set()).add(product_key(current))

    for pkey, pkeys in distance_pairings.items():
        bpp = bpp_by_product_key[pkey]
        bpp.data[COLUMNS.pairing] = ';'.join(sorted(pkeys))

    for pkey, pkeys in product_pairings.items():
        bpp = bpp_by_product_key[pkey]
        bpp.data[COLUMNS.inferred_pairing] = ';'.join(sorted(pkeys))

    fname = os.path.join(
        output,
        'mavis_paired_{}.tab'.format('_'.join(sorted(list(libraries))))
    )
    output_tabbed_file(bpps, fname)
