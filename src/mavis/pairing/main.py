import itertools
import os
import time
from typing import Dict, List, Set, Tuple

import pandas as pd

from ..annotate.file_io import ReferenceFile
from ..breakpoint import BreakpointPair
from ..constants import CALL_METHOD, COLUMNS, SPLICE_TYPE, SVTYPE
from ..util import generate_complete_stamp, logger, output_tabbed_file, read_inputs
from .pairing import inferred_equivalent, pair_by_distance, product_key


def main(
    inputs: List[str],
    output: str,
    config: Dict,
    start_time=int(time.time()),
):
    """
    Args:
        inputs (List[str]): list of input files to read
        output (str): path to the output directory
    """
    annotations = ReferenceFile.load_from_config(config, 'annotations', eager_load=True)

    # load the file
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
                COLUMNS.annotation_id,
                COLUMNS.library,
                COLUMNS.fusion_cdna_coding_start,
                COLUMNS.fusion_cdna_coding_end,
                COLUMNS.fusion_sequence_fasta_id,
            ],
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
    logger.info(f'read {len(bpps)} breakpoint pairs')

    # load all transcripts
    reference_transcripts = dict()
    for genes in annotations.content.values():
        for gene in genes:
            for unspliced_t in gene.transcripts:
                if unspliced_t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, unspliced_t)
                reference_transcripts[unspliced_t.name] = unspliced_t

    # map the calls by library and ensure there are no name/key conflicts
    calls_by_cat: Dict[Tuple[str, str, bool, str], List[BreakpointPair]] = dict()
    calls_by_ann: Dict[Tuple[str, str], List[BreakpointPair]] = dict()
    bpp_by_product_key: Dict[str, BreakpointPair] = dict()
    libraries = set()

    # initialize the pairing mappings
    for bpp in bpps:
        libraries.add(bpp.library)
        category = (bpp.break1.chr, bpp.break2.chr, bpp.opposing_strands, bpp.event_type)
        bpp.data[COLUMNS.product_id] = product_key(bpp)
        calls_by_cat.setdefault(category, []).append(bpp)
        if bpp.data.get(COLUMNS.gene1) or bpp.data.get(COLUMNS.gene2):
            calls_by_ann.setdefault(
                (bpp.data.get(COLUMNS.transcript1), bpp.data.get(COLUMNS.transcript2)), []
            ).append(bpp)
        bpp.data[COLUMNS.pairing] = ''
        bpp.data[COLUMNS.inferred_pairing] = ''

        if product_key(bpp) in bpp_by_product_key:
            diffs = {}
            other = bpp_by_product_key[product_key(bpp)]
            for key in (set(other.data.keys()) | set(bpp.data.keys())) - {'line_no'}:
                if bpp.data.get(key) != other.data.get(key):
                    diffs[key] = (bpp.data.get(key), other.data.get(key))
            if diffs:
                raise KeyError(
                    f'duplicate bpp ({product_key(bpp)}) is not unique within lib (diffs: {diffs})'
                )
        bpp_by_product_key[product_key(bpp)] = bpp

    distance_pairings: Dict[str, Set[str]] = {}
    product_pairings: Dict[str, Set[str]] = {}
    logger.info('computing distance based pairings')
    # pairwise comparison of breakpoints between all libraries
    for set_num, (category, calls) in enumerate(
        sorted(calls_by_cat.items(), key=lambda x: (len(x[1]), x[0]), reverse=True)
    ):
        logger.info(f'comparing set {set_num + 1} of {len(calls_by_cat)} with {len(calls)} items')
        for node, adj_list in pair_by_distance(calls, distances, against_self=False).items():
            distance_pairings.setdefault(node, set()).update(adj_list)

    logger.info('computing inferred (by product) pairings')
    for calls in calls_by_ann.values():
        calls_by_lib: Dict[str, List[BreakpointPair]] = {}
        for call in calls:
            calls_by_lib.setdefault(call.library, []).append(call)

        for lib, other_lib in itertools.combinations(calls_by_lib.keys(), 2):
            # create combinations from other libraries in the same category
            pairs = calls_by_lib[lib]
            other_pairs = calls_by_lib[other_lib]

            for current, other in itertools.product(pairs, other_pairs):
                if inferred_equivalent(
                    current, other, distances=distances, reference_transcripts=reference_transcripts
                ):
                    product_pairings.setdefault(product_key(current), set()).add(product_key(other))
                    product_pairings.setdefault(product_key(other), set()).add(product_key(current))

    for pkey, pkeys in distance_pairings.items():
        bpp = bpp_by_product_key[pkey]
        bpp.data[COLUMNS.pairing] = ';'.join(sorted(pkeys))

    for pkey, pkeys in product_pairings.items():
        bpp = bpp_by_product_key[pkey]
        bpp.data[COLUMNS.inferred_pairing] = ';'.join(sorted(pkeys))

    fname = os.path.join(output, 'mavis_paired.tab')
    output_tabbed_file(bpps, fname)
    generate_complete_stamp(output, start_time=start_time)
