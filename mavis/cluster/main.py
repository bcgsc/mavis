import os
import itertools
from .cluster import cluster_breakpoint_pairs
from ..constants import COLUMNS
from ..interval import Interval
from .constants import DEFAULTS
from ..util import read_inputs, output_tabbed_file, write_bed_file, generate_complete_stamp
from ..util import build_batch_id, filter_on_overlap, log, mkdirp


def main(
    inputs, output, stranded_bam, library, protocol, masking, annotations,
    cluster_clique_size=DEFAULTS.cluster_clique_size,
    cluster_radius=DEFAULTS.cluster_radius,
    uninformative_filter=DEFAULTS.uninformative_filter,
    max_proximity=DEFAULTS.max_proximity,
    min_clusters_per_file=DEFAULTS.min_clusters_per_file,
    max_files=DEFAULTS.max_files,
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

    # output files
    cluster_batch_id = build_batch_id(prefix='cluster-')
    UNINFORM_OUTPUT = os.path.join(output, 'uninformative_clusters.txt')
    CLUSTER_ASSIGN_OUTPUT = os.path.join(output, 'cluster_assignment.tab')
    CLUSTER_BED_OUTPUT = os.path.join(output, 'clusters.bed')
    COMPLETE_STAMP = os.path.join(output, 'CLUSTERING.COMPLETE')
    split_file_name_func = lambda x: os.path.join(output, '{}-{}.tab'.format(cluster_batch_id, x))
    # load the input files
    all_breakpoint_pairs = read_inputs(
        inputs,
        cast={COLUMNS.tools: lambda x: set(x.split(';')) if x else set()},
        add={COLUMNS.library: library, COLUMNS.protocol: protocol},
        expand_ns=True, explicit_strand=False
    )
    # ignore other library inputs
    other_libs = set()
    unfiltered_breakpoint_pairs = []
    for bpp in all_breakpoint_pairs:
        if bpp.library != library:
            other_libs.add(bpp.library)
        else:
            unfiltered_breakpoint_pairs.append(bpp)
    if len(other_libs) > 0:
        log('warning: ignoring breakpoints found for other libraries:', sorted([l for l in other_libs]))

    # filter by masking file
    breakpoint_pairs, filtered_bpp = filter_on_overlap(unfiltered_breakpoint_pairs, masking)

    log('computing clusters')
    clusters = cluster_breakpoint_pairs(breakpoint_pairs, r=cluster_radius, k=cluster_clique_size)

    hist = {}
    length_hist = {}
    for index, cluster in enumerate(clusters):
        input_pairs = clusters[cluster]
        hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
        c1 = round(len(cluster[0]), -2)
        c2 = round(len(cluster[1]), -2)
        length_hist[c1] = length_hist.get(c1, 0) + 1
        length_hist[c2] = length_hist.get(c2, 0) + 1
        cluster.data[COLUMNS.cluster_id] = '{}-{}'.format(cluster_batch_id, index + 1)
        cluster.data[COLUMNS.cluster_size] = len(input_pairs)
        temp = set()
        for p in input_pairs:
            temp.update(p.data[COLUMNS.tools])
        cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
    log('computed', len(clusters), 'clusters', time_stamp=False)
    log('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
    log('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    mkdirp(output)

    with open(CLUSTER_ASSIGN_OUTPUT, 'w') as fh:
        header = set()
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
            row[COLUMNS.library] = library
            row[COLUMNS.protocol] = protocol
        output_tabbed_file(rows.values(), CLUSTER_ASSIGN_OUTPUT)

    output_files = []
    # filter clusters based on annotations
    # decide on the number of clusters to validate per job
    pass_clusters = list(clusters)
    fail_clusters = []

    if uninformative_filter:
        pass_clusters = []
        for cluster in clusters:
            # loop over the annotations
            overlaps_gene = False
            w1 = Interval(cluster.break1.start - max_proximity, cluster.break1.end + max_proximity)
            w2 = Interval(cluster.break2.start - max_proximity, cluster.break2.end + max_proximity)
            for gene in annotations.get(cluster.break1.chr, []):
                if Interval.overlaps(gene, w1):
                    overlaps_gene = True
                    break
            for gene in annotations.get(cluster.break2.chr, []):
                if Interval.overlaps(gene, w2):
                    overlaps_gene = True
                    break
            if overlaps_gene:
                pass_clusters.append(cluster)
            else:
                fail_clusters.append(cluster)
    if len(fail_clusters) + len(pass_clusters) != len(clusters):
        raise AssertionError(
            'totals do not add up', len(fail_clusters), len(pass_clusters), 'does not total to', len(clusters))
    log('filtered', len(fail_clusters), 'clusters as not informative')
    output_tabbed_file(fail_clusters, UNINFORM_OUTPUT)

    JOB_SIZE = min_clusters_per_file
    if len(pass_clusters) // min_clusters_per_file > max_files - 1:
        JOB_SIZE = len(pass_clusters) // max_files
        assert(len(pass_clusters) // JOB_SIZE == max_files)

    bedfile = os.path.join(output, 'clusters.bed')
    write_bed_file(bedfile, itertools.chain.from_iterable([b.get_bed_repesentation() for b in pass_clusters]))
    job_ranges = list(range(0, len(pass_clusters), JOB_SIZE))
    if len(job_ranges) == 0 or job_ranges[-1] != len(pass_clusters):
        job_ranges.append(len(pass_clusters))
    job_ranges = zip(job_ranges, job_ranges[1::])

    for i, jrange in enumerate(job_ranges):
        # generate an output file
        filename = split_file_name_func(i + 1)
        output_files.append(filename)
        output_tabbed_file(pass_clusters[jrange[0]:jrange[1]], filename)
    
    generate_complete_stamp(output, log)
    return output_files
