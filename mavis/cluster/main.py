import inspect
import itertools
import os
from shortuuid import uuid
import time

from .cluster import merge_breakpoint_pairs
from .constants import DEFAULTS
from ..constants import COLUMNS
from ..util import filter_on_overlap, filter_uninformative, generate_complete_stamp, LOG, log_arguments, mkdirp, output_tabbed_file, read_inputs, write_bed_file


def split_clusters(clusters, outputdir, batch_id, min_clusters_per_file=0, max_files=1, write_bed_summary=True):
    """
    For a set of clusters creates a bed file representation of all clusters.
    Also splits the clusters evenly into multiple files based on the user parameters (min_clusters_per_file, max_files)

    Returns:
        list: of output file names (not including the bed file)
    """
    if write_bed_summary:
        bedfile = os.path.join(outputdir, 'clusters.bed')
        write_bed_file(bedfile, itertools.chain.from_iterable([b.get_bed_repesentation() for b in clusters]))

    number_of_jobs = len(clusters) // min_clusters_per_file
    if number_of_jobs > max_files:
        number_of_jobs = max_files
    elif number_of_jobs == 0:
        number_of_jobs = 1

    jobs = [[] for j in range(0, number_of_jobs)]
    clusters = sorted(clusters, key=lambda x: (x.break1.chr, x.break1.start, x.break2.chr, x.break2.start))

    # split up consecutive clusters
    for i, cluster in enumerate(clusters):
        jobs[i % len(jobs)].append(cluster)

    assert sum([len(j) for j in jobs]) == len(clusters)
    output_files = []
    for i, job in enumerate(jobs):
        # generate an output file
        filename = os.path.join(outputdir, '{}-{}.tab'.format(batch_id, i + 1))
        output_files.append(filename)
        output_tabbed_file(job, filename)
    return output_files


def main(
    inputs, output, strand_specific, library, protocol, disease_status, masking, annotations,
    limit_to_chr=DEFAULTS.limit_to_chr,
    cluster_initial_size_limit=DEFAULTS.cluster_initial_size_limit,
    cluster_radius=DEFAULTS.cluster_radius,
    uninformative_filter=DEFAULTS.uninformative_filter,
    max_proximity=DEFAULTS.max_proximity,
    min_clusters_per_file=DEFAULTS.min_clusters_per_file,
    max_files=DEFAULTS.max_files,
    batch_id=None,
    split_only=False,
    start_time=int(time.time()),
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        strand_specific (bool): is the bam using a strand specific protocol
        library (str): the library to look for in each of the input files
        protocol (PROTOCOL): the sequence protocol (genome or transcriptome)
        masking (object): see :func:`~mavis.annotate.file_io.load_masking_regions`
        cluster_clique_size (int): the maximum size of cliques to search for using the exact algorithm
        cluster_radius (int): distance (in breakpoint pairs) used in deciding to join bpps in a cluster
        uninformative_filter (bool): if True then clusters should be filtered out if they are not
          within a specified (max_proximity) distance to any annotation
        max_proximity (int): the maximum distance away an annotation can be before the uninformative_filter
          is applied
        annotations (ReferenceFile): see :func:`~mavis.annotate.file_io.load_reference_genes`
        min_clusters_per_file (int): the minimum number of clusters to output to a file
        max_files (int): the maximum number of files to split clusters into
    """
    if uninformative_filter:
        annotations.load()
    if masking:
        masking.load()

    # output files
    batch_id = 'batch-' + str(uuid()) if batch_id is None else batch_id
    filtered_output = os.path.join(output, 'filtered_pairs.tab')
    cluster_assign_output = os.path.join(output, 'cluster_assignment.tab')

    # load the input files
    breakpoint_pairs = read_inputs(
        inputs,
        cast={COLUMNS.tools: lambda x: set(x.split(';')) if x else set() if not split_only else x},
        add_default={
            COLUMNS.library: library,
            COLUMNS.protocol: protocol,
            COLUMNS.tools: '',
            COLUMNS.disease_status: disease_status,
            COLUMNS.stranded: False,
            COLUMNS.tracking_id: ''
        },
        expand_strand=False, expand_orient=True, expand_svtype=True
    )
    # filter any breakpoint pairs where the library and protocol don't match
    other_libs = set()
    other_chr = set()
    unfiltered_breakpoint_pairs = []
    filtered_pairs = []
    LOG('filtering by library and chr name')
    for bpp in breakpoint_pairs:
        if bpp.library is None:
            bpp.library = library
        if bpp.library != library:
            other_libs.add(bpp.library)
            bpp.data[COLUMNS.filter_comment] = 'Not the target library name'
            filtered_pairs.append(bpp)
        elif None in limit_to_chr or (bpp.break1.chr in limit_to_chr and bpp.break2.chr in limit_to_chr):
            unfiltered_breakpoint_pairs.append(bpp)
        else:
            other_chr.update({bpp.break1.chr, bpp.break2.chr})
            bpp.data[COLUMNS.filter_comment] = 'Non standard chromosome name'
            filtered_pairs.append(bpp)
    other_chr -= set(limit_to_chr)
    breakpoint_pairs = unfiltered_breakpoint_pairs
    if other_libs:
        LOG('warning: ignoring breakpoints found for other libraries:', sorted([l for l in other_libs]))
    if other_chr:
        LOG('warning: filtered events on chromosomes', other_chr)
    # filter by masking file
    breakpoint_pairs, masked_pairs = filter_on_overlap(breakpoint_pairs, masking.content)
    for bpp in masked_pairs:
        filtered_pairs.append(bpp)
    # filter by informative
    if uninformative_filter:
        LOG('filtering from', len(breakpoint_pairs), 'breakpoint pairs using informative filter')
        pass_clusters, uninformative_clusters = filter_uninformative(annotations.content, breakpoint_pairs, max_proximity=max_proximity)
        LOG(
            'filtered from', len(breakpoint_pairs),
            'down to', len(pass_clusters),
            '(removed {})'.format(len(uninformative_clusters))
        )
        breakpoint_pairs = pass_clusters
        for bpp in uninformative_clusters:
            bpp.data[COLUMNS.filter_comment] = 'Uninformative'
            filtered_pairs.append(bpp)
    else:
        LOG('did not apply uninformative filter')

    output_tabbed_file(filtered_pairs, filtered_output)
    mkdirp(output)

    if not split_only:
        LOG('computing clusters')
        clusters = merge_breakpoint_pairs(
            breakpoint_pairs, cluster_radius=cluster_radius, cluster_initial_size_limit=cluster_initial_size_limit)

        hist = {}
        length_hist = {}
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
        LOG('computed', len(clusters), 'clusters', time_stamp=False)
        LOG('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
        LOG('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
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
        batch_id,
        min_clusters_per_file=min_clusters_per_file,
        max_files=max_files,
        write_bed_summary=True
    )

    generate_complete_stamp(output, LOG, start_time=start_time, prefix='MAVIS-{}.'.format(batch_id))
    return output_files
