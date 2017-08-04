import os
import itertools
from .cluster import merge_breakpoint_pairs
from ..constants import COLUMNS
from .constants import DEFAULTS
from ..util import read_inputs, output_tabbed_file, write_bed_file, generate_complete_stamp
from ..util import filter_on_overlap, log, mkdirp, filter_uninformative, log_arguments
import uuid
import inspect


def main(
    inputs, output, stranded_bam, library, protocol, disease_status, masking, annotations,
    limit_to_chr=DEFAULTS.limit_to_chr,
    cluster_initial_size_limit=DEFAULTS.cluster_initial_size_limit,
    cluster_radius=DEFAULTS.cluster_radius,
    uninformative_filter=DEFAULTS.uninformative_filter,
    max_proximity=DEFAULTS.max_proximity,
    min_clusters_per_file=DEFAULTS.min_clusters_per_file,
    max_files=DEFAULTS.max_files,
    fetch_method_individual=True,
    log_args=False,
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        stranded_bam (bool): is the bam using a strand specific protocol
        library (str): the library to look for in each of the input files
        protocol (PROTOCOL): the sequence protocol (genome or transcriptome)
        masking (object): see :func:`~mavis.annotate.file_io.load_masking_regions`
        cluster_clique_size (int): the maximum size of cliques to search for using the exact algorithm
        cluster_radius (int): distance (in breakpoint pairs) used in deciding to join bpps in a cluster
        uninformative_filter (bool): if True then clusters should be filtered out if they are not
          within a specified (max_proximity) distance to any annotation
        max_proximity (int): the maximum distance away an annotation can be before the uninformative_filter
          is applied
        annotations (object): see :func:`~mavis.annotate.file_io.load_reference_genes`
        min_clusters_per_file (int): the minimum number of clusters to output to a file
        max_files (int): the maximum number of files to split clusters into
    """
    if log_args:
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)
        args = {arg: values[arg] for arg in args if arg != 'log_args'}
        log_arguments(args)

    # output files
    cluster_batch_id = 'batch-' + str(uuid.uuid4())
    UNINFORM_OUTPUT = os.path.join(output, 'uninformative_clusters.txt')
    CLUSTER_ASSIGN_OUTPUT = os.path.join(output, 'cluster_assignment.tab')
    # TODO: CLUSTER_BED_OUTPUT = os.path.join(output, 'clusters.bed')

    def split_file_name_func(x):
        return os.path.join(output, '{}-{}.tab'.format(cluster_batch_id, x))
    # load the input files
    breakpoint_pairs = read_inputs(
        inputs,
        cast={COLUMNS.tools: lambda x: set(x.split(';')) if x else set()},
        add={
            COLUMNS.library: library,
            COLUMNS.protocol: protocol,
            COLUMNS.tools: '',
            COLUMNS.disease_status: disease_status
        },
        expand_ns=True, explicit_strand=False
    )
    # filter against chr and ignore other library inputs
    other_libs = set()
    other_chr = set()
    unfiltered_breakpoint_pairs = []
    log('filtering by library and chr name')
    for bpp in breakpoint_pairs:
        if bpp.library is None:
            bpp.library = library
        if bpp.library != library:
            other_libs.add(bpp.library)
        elif bpp.break1.chr in limit_to_chr and bpp.break2.chr in limit_to_chr:
            unfiltered_breakpoint_pairs.append(bpp)
        else:
            other_chr.update({bpp.break1.chr, bpp.break2.chr})
    other_chr -= set(limit_to_chr)
    breakpoint_pairs = unfiltered_breakpoint_pairs
    if len(other_libs) > 0:
        log('warning: ignoring breakpoints found for other libraries:', sorted([l for l in other_libs]))
    if len(other_chr) > 0:
        log('warning: filtered events on chromosomes not found in "limit_to_chr"', other_chr)
    # filter by masking file
    breakpoint_pairs, filtered_bpp = filter_on_overlap(breakpoint_pairs, masking)
    # filter by informative
    if uninformative_filter:
        log('filtering from', len(breakpoint_pairs), 'breakpoint pairs using informative filter')
        pass_clusters, uninformative_clusters = filter_uninformative(annotations, breakpoint_pairs)
        log(
            'filtered from', len(breakpoint_pairs),
            'down to', len(pass_clusters),
            '(removed {})'.format(len(uninformative_clusters))
        )
        breakpoint_pairs = pass_clusters
        output_tabbed_file(uninformative_clusters, UNINFORM_OUTPUT)
    else:
        log('did not apply uninformative filter')
    log('computing clusters')
    clusters = merge_breakpoint_pairs(
        breakpoint_pairs, cluster_radius=cluster_radius, cluster_initial_size_limit=cluster_initial_size_limit)

    hist = {}
    length_hist = {}
    for index, cluster in enumerate(clusters):
        input_pairs = clusters[cluster]
        hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
        c1 = round(len(cluster[0]), -2)
        c2 = round(len(cluster[1]), -2)
        length_hist[c1] = length_hist.get(c1, 0) + 1
        length_hist[c2] = length_hist.get(c2, 0) + 1
        cluster.data[COLUMNS.cluster_id] = str(uuid.uuid4())
        cluster.data[COLUMNS.cluster_size] = len(input_pairs)
        temp = set()
        data_items = set()
        for p in input_pairs:
            temp.update(p.data[COLUMNS.tools])
            data_items.update(p.data.keys())
        cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
        data_items -= {COLUMNS.tools}
        # retain all data where data is consistent between the input pairs
        for item in data_items:
            s = [p.data.get(item, None) for p in input_pairs]
            s = set(s)
            if len(s) == 1:
                cluster.data[item] = list(s)[0]
    log('computed', len(clusters), 'clusters', time_stamp=False)
    log('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
    log('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    mkdirp(output)

    rows = {}
    for cluster, input_pairs in clusters.items():
        for p in input_pairs:
            if p not in rows:
                rows[p] = p.flatten()
            rows[p][COLUMNS.tools].update(p.data[COLUMNS.tools])
            rows[p].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id])
    for row in rows.values():
        row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
        row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools])))
    output_tabbed_file(rows.values(), CLUSTER_ASSIGN_OUTPUT)

    output_files = []
    # filter clusters based on annotations
    # decide on the number of clusters to validate per job
    clusters = list(clusters.keys())
    JOB_SIZE = min_clusters_per_file
    if len(clusters) // min_clusters_per_file > max_files - 1:
        JOB_SIZE = int(round(len(clusters) / max_files, 0))
        assert(len(clusters) // JOB_SIZE <= max_files)

    bedfile = os.path.join(output, 'clusters.bed')
    write_bed_file(bedfile, itertools.chain.from_iterable([b.get_bed_repesentation() for b in clusters]))
    number_of_jobs = len(clusters) // min_clusters_per_file
    if number_of_jobs >= max_files:
        number_of_jobs = max_files
    elif number_of_jobs == 0:
        number_of_jobs = 1

    jobs = [[] for j in range(0, number_of_jobs)]
    clusters = sorted(clusters, key=lambda x: (x.break1.chr, x.break1.start, x.break2.chr, x.break2.start))

    if fetch_method_individual:
        # split up consecutive clusters
        for i, cluster in enumerate(clusters):
            jid = i % len(jobs)
            jobs[jid].append(cluster)
    else:  # group consecutive clusters
        extras = len(clusters) % number_of_jobs
        cluster_per_job = len(clusters) // number_of_jobs
        i = 0
        for job in jobs:
            job.extend(clusters[i:i + cluster_per_job + (1 if extras > 0 else 0)])
            extras -= 1
            i += len(job)
    assert(sum([len(j) for j in jobs]) == len(clusters))
    for i, job in enumerate(jobs):
        # generate an output file
        filename = split_file_name_func(i + 1)
        output_files.append(filename)
        output_tabbed_file(job, filename)

    generate_complete_stamp(output, log)
    return output_files
