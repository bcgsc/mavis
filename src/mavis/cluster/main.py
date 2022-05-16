import itertools
import os
import time
from typing import Dict, List

from mavis_config.constants import SUBCOMMAND
from shortuuid import uuid

from ..annotate.file_io import ReferenceFile
from ..breakpoint import BreakpointPair
from ..constants import COLUMNS
from ..util import (
    filter_on_overlap,
    filter_uninformative,
    generate_complete_stamp,
    logger,
    mkdirp,
    output_tabbed_file,
    read_inputs,
    write_bed_file,
)
from .cluster import merge_breakpoint_pairs

SECTION = SUBCOMMAND.CLUSTER


def split_clusters(
    clusters: List[BreakpointPair],
    outputdir: str,
    total_batches: int,
    write_bed_summary: bool = True,
) -> List[str]:
    """
    For a set of clusters creates a bed file representation of all clusters.
    Also splits the clusters evenly into multiple files based on the user parameters (max_files)

    Returns:
        list of output file names (not including the bed file)
    """
    if write_bed_summary:
        bedfile = os.path.join(outputdir, 'clusters.bed')
        write_bed_file(
            bedfile, itertools.chain.from_iterable([b.get_bed_repesentation() for b in clusters])
        )

    jobs: List[List[BreakpointPair]] = [[] for j in range(0, total_batches)]
    clusters = sorted(
        clusters, key=lambda x: (x.break1.chr, x.break1.start, x.break2.chr, x.break2.start)
    )

    # split up consecutive clusters
    for i, cluster in enumerate(clusters):
        jobs[i % len(jobs)].append(cluster)

    assert sum([len(j) for j in jobs]) == len(clusters)
    output_files = []
    for i, job in enumerate(jobs):
        # generate an output file
        filename = os.path.join(outputdir, 'batch-{}.tab'.format(i + 1))
        output_files.append(filename)
        output_tabbed_file(job, filename)
    return output_files


def main(
    inputs: List[str],
    output: str,
    library: str,
    config: Dict,
    start_time=int(time.time()),
    **kwargs,
):
    """
    Args:
        inputs: list of input files to read
        output: path to the output directory
        library: the library to look for in each of the input files
        masking (ReferenceFile): see :func:`mavis.annotate.file_io.load_masking_regions`
        annotations (ReferenceFile): see :func:`mavis.annotate.file_io.load_annotations`
    """
    masking = ReferenceFile.load_from_config(config, 'masking', eager_load=True)
    annotations = ReferenceFile.load_from_config(config, 'annotations')

    if config[f'{SECTION}.uninformative_filter'] and not annotations.is_empty():
        annotations.load()
    if not masking.is_empty():
        masking.load()

    lib_config = config['libraries'][library]

    # output files
    filtered_output = os.path.join(output, 'filtered_pairs.tab')
    cluster_assign_output = os.path.join(output, 'cluster_assignment.tab')

    # load the input files
    breakpoint_pairs = read_inputs(
        inputs,
        apply={
            COLUMNS.tools: lambda x: set(x.split(';'))
            if x
            else set()
            if not config[f'{SECTION}.split_only']
            else x
        },
        add_default={
            COLUMNS.library: library,
            COLUMNS.protocol: lib_config['protocol'],
            COLUMNS.tools: '',
            COLUMNS.disease_status: lib_config['disease_status'],
            COLUMNS.stranded: False,
            COLUMNS.tracking_id: '',
        },
        expand_strand=False,
        expand_orient=True,
        expand_svtype=True,
    )
    # filter any breakpoint pairs where the library and protocol don't match
    other_libs = set()
    other_chr = set()
    unfiltered_breakpoint_pairs = []
    filtered_pairs = []
    logger.info('filtering by library and chr name')
    for bpp in breakpoint_pairs:
        if bpp.library is None:
            bpp.library = library
        if bpp.library != library:
            other_libs.add(bpp.library)
            bpp.data[COLUMNS.filter_comment] = 'Not the target library name'
            filtered_pairs.append(bpp)
        elif not config[f'{SECTION}.limit_to_chr'] or (
            bpp.break1.chr in config[f'{SECTION}.limit_to_chr']
            and bpp.break2.chr in config[f'{SECTION}.limit_to_chr']
        ):
            unfiltered_breakpoint_pairs.append(bpp)
        else:
            other_chr.update({bpp.break1.chr, bpp.break2.chr})
            bpp.data[COLUMNS.filter_comment] = 'Non standard chromosome name'
            filtered_pairs.append(bpp)
    if config[f'{SECTION}.limit_to_chr']:
        other_chr -= set(config[f'{SECTION}.limit_to_chr'])
    breakpoint_pairs = unfiltered_breakpoint_pairs
    if other_libs:
        logger.info(
            f'warning: ignoring breakpoints found for other libraries: {sorted([lib for lib in other_libs])}',
        )
    if other_chr:
        logger.info(f'warning: filtered events on chromosomes {other_chr}')
    # filter by masking file
    breakpoint_pairs, masked_pairs = filter_on_overlap(breakpoint_pairs, masking.content)
    for bpp in masked_pairs:
        filtered_pairs.append(bpp)
    # filter by informative
    if config[f'{SECTION}.uninformative_filter']:
        logger.info(
            f'filtering from {len(breakpoint_pairs)} breakpoint pairs using informative filter'
        )
        pass_clusters, uninformative_clusters = filter_uninformative(
            annotations.content, breakpoint_pairs, max_proximity=config[f'{SECTION}.max_proximity']
        )
        logger.info(
            f'filtered from {len(breakpoint_pairs)} down to {len(pass_clusters)} (removed {len(uninformative_clusters)})'
        )
        breakpoint_pairs = pass_clusters
        for bpp in uninformative_clusters:
            bpp.data[COLUMNS.filter_comment] = 'Uninformative'
            filtered_pairs.append(bpp)
    else:
        logger.info('did not apply uninformative filter')

    mkdirp(output)
    output_tabbed_file(filtered_pairs, filtered_output)

    if not config[f'{SECTION}.split_only']:
        logger.info('computing clusters')
        clusters = merge_breakpoint_pairs(
            breakpoint_pairs,
            cluster_radius=config[f'{SECTION}.cluster_radius'],
            cluster_initial_size_limit=config[f'{SECTION}.cluster_initial_size_limit'],
        )

        hist: Dict[int, int] = {}
        length_hist: Dict[float, int] = {}

        for cluster in clusters:
            input_pairs = clusters[cluster]
            hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
            cluster1 = round(len(cluster[0]), -2)
            cluster2 = round(len(cluster[1]), -2)
            length_hist[cluster1] = length_hist.get(cluster1, 0) + 1
            length_hist[cluster2] = length_hist.get(cluster2, 0) + 1
            cluster.data[COLUMNS.cluster_id] = str(uuid())
            cluster.data[COLUMNS.cluster_size] = len(input_pairs)
            temp = set()
            data_items = set()
            combined_tracking_id = set()  # group the tracking ids
            for pair in input_pairs:
                temp.update(pair.data[COLUMNS.tools])
                data_items.update(pair.data.keys())
                if COLUMNS.tracking_id in pair.data and pair.tracking_id:
                    combined_tracking_id.update(pair.tracking_id.split(';'))
            cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
            cluster.data[COLUMNS.tracking_id] = ';'.join(sorted(list(combined_tracking_id)))

            data_items -= {COLUMNS.tools, COLUMNS.tracking_id}
            # retain all data where data is consistent between the input pairs
            for item in data_items:
                common_data = [p.data.get(item, None) for p in input_pairs]
                common_data = set(common_data)
                if len(common_data) == 1:
                    cluster.data[item] = list(common_data)[0]
        logger.info(f'computed {len(clusters)} clusters')
        logger.info(f'cluster input pairs distribution {sorted(hist.items())}')
        logger.info(f'cluster intervals lengths {sorted(length_hist.items())}')
        # map input pairs to cluster ids
        # now create the mapping from the original input files to the cluster(s)

        rows = {}
        for cluster, input_pairs in clusters.items():
            for pair in input_pairs:
                if pair not in rows:
                    rows[pair] = pair.flatten()
                rows[pair][COLUMNS.tools].update(pair.data[COLUMNS.tools])
                rows[pair].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id])
        for row in rows.values():
            row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
            row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools])))
        output_tabbed_file(rows.values(), cluster_assign_output)
        breakpoint_pairs = list(clusters.keys())

    output_files = split_clusters(
        breakpoint_pairs,
        output,
        total_batches=lib_config['total_batches'],
        write_bed_summary=True,
    )

    generate_complete_stamp(output, start_time=start_time)
    return output_files
