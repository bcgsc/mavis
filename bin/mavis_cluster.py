"""
About
--------

This is the first step (other than preprocessing inputs) in the svmerge pipeline. Input files are taken in, separated by
library and protocol, and then clustered if they are similar types of event calls and are close together. The output is
a list of estimated calls based on the median of each cluster

This script is also responsible for setting up the directory structure of the outputs, which will be in the following
pattern

::

    <output_dir_name>/
    |-- clustering/
    |   `-- <library>_<protocol>/
    |       |-- uninformative_clusters.txt
    |       |-- clusters.bed
    |       |-- cluster_assignment.tab
    |       `-- clusterset-#.tab
    |-- validation/
    |   `-- <library>_<protocol>/
    |       |-- qsub.sh
    |       `-- log/
    |-- annotation/
    |   `--<library>_<protocol>/
    |-- pairing/
    `-- summary/


Filtering
------------

Clusters are optionally post-filtered based on a gene annotation file. This reduces the number of clusters that need to
go through validation which is the most expensive in terms of time
"""

import os
import sys
import argparse
import re
from datetime import datetime
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.constants import *
from mavis.error import *
from mavis.interval import Interval
from mavis.breakpoint import BreakpointPair, read_bpp_from_input_file
from mavis.cluster import cluster_breakpoint_pairs
from mavis.annotate import load_reference_genes, load_masking_regions
from mavis import __version__
from mavis.constants import log, build_batch_id

__prog__ = os.path.basename(os.path.realpath(__file__))
MIN_CLUSTERS_PER_FILE = 50
MAX_FILES = 100
CLUSTER_RADIUS = 20
CLUSTER_CLIQUE_SIZE = 15
MAX_PROXIMITY = 5000


def write_bed_file(filename, cluster_breakpoint_pairs):
    with open(filename, 'w') as fh:
        for bpp in cluster_breakpoint_pairs:
            if bpp.interchromosomal:
                fh.write('{}\t{}\t{}\tcluster={}\n'.format(
                    bpp.break1.chr, bpp.break1.start, bpp.break1.end, bpp.data['cluster_id']))
                fh.write('{}\t{}\t{}\tcluster={}\n'.format(
                    bpp.break2.chr, bpp.break2.start, bpp.break2.end, bpp.data['cluster_id']))
            else:
                fh.write('{}\t{}\t{}\tcluster={}\n'.format(
                    bpp.break1.chr, bpp.break1.start, bpp.break2.end, bpp.data['cluster_id']))


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number')
    parser.add_argument(
        '-f', '--force_overwrite', default=False, type=TSV.tsv_boolean,
        help='set flag to overwrite existing reviewed files')

    parser.add_argument(
        '-n', '--inputs', help='path to the input files', required=True, action='append')
    parser.add_argument('library', help='library name')
    parser.add_argument('protocol', help='the library protocol: genome or transcriptome', choices=PROTOCOL.values())
    parser.add_argument('output', help='path to the output directory')
    g = parser.add_argument_group('file splitting options')
    g.add_argument(
        '--max_files', '-j', default=MAX_FILES, type=int, dest='max_files',
        help='defines the maximum number of files that can be created')
    g.add_argument(
        '--min_clusters_per_file', '-e', default=MIN_CLUSTERS_PER_FILE, type=int,
        help='defines the minimum number of clusters per file')
    parser.add_argument(
        '-r', '--cluster_radius', help='radius to use in clustering', default=CLUSTER_RADIUS, type=int)
    parser.add_argument(
        '-k', '--cluster_clique_size',
        help='parameter used for computing cliques, smaller is faster, above 20 will be slow',
        default=CLUSTER_CLIQUE_SIZE, type=int)
    g = parser.add_argument_group('filter arguments')
    g.add_argument(
        '--uninformative_filter', default=True, help='If flag is False then the clusters will not be filtered '
        'based on lack of annotation', type=TSV.tsv_boolean)
    g.add_argument(
        '--max_proximity', '-p', type=int, default=MAX_PROXIMITY, help='maximum distance to look for annotations'
        'from evidence window')
    g.add_argument(
        '--annotations', '-a',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_annotations_20170203.json',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.')
    g.add_argument(
        '--masking', '-m',
        default='/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv'
    )
    args = parser.parse_args()

    if args.min_clusters_per_file < 1:
        print('\nerror: min_clusters_per_file cannot be less than 1')
        parser.print_help()
        exit(1)

    if args.max_files < 1:
        print('\nerror: max_files cannot be less than 1')
        parser.print_help()
        exit(1)

    if os.path.exists(args.output) and not args.force_overwrite:
        print(
            '\nerror: output directory {0} already exists. please use the --force_overwrite option'.format(args.output))
        parser.print_help()
        exit(1)

    for f in args.inputs:
        if not os.path.exists(f):
            print('\nerror: input file {0} does not exist'.format(f))
            parser.print_help()
            exit(1)

    args.output = os.path.abspath(args.output)

    return args


def main(args):
    # load the input files
    breakpoint_pairs = []
    masks = args.masking
    if isinstance(args.masking, str):
        log('loading:', args.masking)
        masks = load_masking_regions(args.masking)
        args.masking = masks
    mask_filtered = 0
    for f in args.inputs:
        log('loading:', f)
        bpps = read_bpp_from_input_file(
            f,
            validate={
                COLUMNS.tools: '^(\S+_v?\d+\.\d+\.\d+)(;\S+_v?\d+\.\d+\.\d+)*$',
                COLUMNS.library: '^[\w-]+$'
            },
            _in={
                COLUMNS.protocol: PROTOCOL
            },
            add={
                COLUMNS.protocol: args.protocol,
                COLUMNS.library: args.library
            },
            force_stranded=True if args.stranded_bam else False
        )
        for bpp in bpps:
            bpp.data[COLUMNS.tools] = set(';'.split(bpp.data[COLUMNS.tools]))
            bpp.data['files'] = {f}
            # filter if the breakpoint call is not in the library we are looking for
            # filter if overlaps a masking region
            if bpp.data[COLUMNS.library] == args.library and bpp.data[COLUMNS.protocol] == args.protocol:
                if not any([Interval.overlaps(m, bpp.break1) for m in masks.get(bpp.break1.chr, [])]) and \
                        not any([Interval.overlaps(m, bpp.break2) for m in masks.get(bpp.break2.chr, [])]):
                    breakpoint_pairs.append(bpp)
                else:
                    mask_filtered += 1

    log('loaded {} breakpoint pairs'.format(len(breakpoint_pairs)))
    log('stranded_bam', args.stranded_bam)
    log('filtered', mask_filtered, 'breakpoint pairs based on overlap with a masked region')

    # load the reference annotations for filtering uninformative clusters
    if args.uninformative_filter:
        REFERENCE_GENES = args.annotations
        if isinstance(args.annotations, str):
            log('loading:', args.annotations)
            REFERENCE_GENES = load_reference_genes(args.annotations, verbose=False)
            args.annotations = REFERENCE_GENES

    cluster_id_prefix = build_batch_id()
    cluster_id = 1

    # set up directories
    log('computing clusters')
    clusters = cluster_breakpoint_pairs(breakpoint_pairs, r=args.cluster_radius, k=args.cluster_clique_size)

    hist = {}
    length_hist = {}
    for cluster, input_pairs in clusters.items():
        hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
        c1 = round(len(cluster[0]), -2)
        c2 = round(len(cluster[1]), -2)
        length_hist[c1] = length_hist.get(c1, 0) + 1
        length_hist[c2] = length_hist.get(c2, 0) + 1
        cluster.data[COLUMNS.cluster_id] = 'cluster_{}-{}'.format(cluster_id_prefix, cluster_id)
        temp = set()
        for p in input_pairs:
            temp.update(p.data[COLUMNS.tools])
        cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
        cluster_id += 1
    log('computed', len(clusters), 'clusters', time_stamp=False)
    log('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
    log('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    if not os.path.exists(args.output):
        log("Creating output directory: '{}'".format(args.output))
        os.makedirs(args.output, mode=0o775)
    f = os.path.join(args.output, 'cluster_assignment.tab')
    with open(f, 'w') as fh:
        header = set()
        log('writing:', f)
        rows = {}

        for cluster, input_pairs in clusters.items():
            for p in input_pairs:
                if p in rows:
                    rows[p][COLUMNS.tools].update(p.data[COLUMNS.tools])
                else:
                    rows[p] = BreakpointPair.flatten(p)
                rows[p].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id])
        for row in rows.values():
            row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
            row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools])))
            row[COLUMNS.library] = args.library
            row[COLUMNS.protocol] = args.protocol
            header.update(row.keys())
        header = sort_columns(header)
        fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
        for row in rows.values():
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

    output_files = []
    # filter clusters based on annotations
    # decide on the number of clusters to validate per job
    pass_clusters = []
    fail_clusters = []

    for cluster in clusters:
        # don't need to generate transcriptome windows b/c will default to genome if not in a gene anyway
        if not args.uninformative_filter:
            pass_clusters.append(cluster)
        else:
            # loop over the annotations
            overlaps_gene = False
            w1 = Interval(cluster.break1.start - args.max_proximity, cluster.break1.end + args.max_proximity)
            w2 = Interval(cluster.break2.start - args.max_proximity, cluster.break2.end + args.max_proximity)
            if not cluster.interchromosomal:
                w1 = w1 | w2
            for gene in REFERENCE_GENES.get(cluster.break1.chr, []):
                if Interval.overlaps(gene, w1):
                    overlaps_gene = True
                    break
            if cluster.interchromosomal:
                for gene in REFERENCE_GENES.get(cluster.break2.chr, []):
                    if Interval.overlaps(gene, w2):
                        overlaps_gene = True
                        break
            if overlaps_gene:
                pass_clusters.append(cluster)
            else:
                fail_clusters.append(cluster)
    assert(len(fail_clusters) + len(pass_clusters) == len(clusters))

    log('filtered', len(fail_clusters), 'clusters as not informative')

    JOB_SIZE = args.min_clusters_per_file
    if len(pass_clusters) // args.min_clusters_per_file > args.max_files - 1:
        JOB_SIZE = len(pass_clusters) // args.max_files
        assert(len(pass_clusters) // JOB_SIZE == args.max_files)

    uninform = os.path.join(args.output, 'uninformative_clusters.txt')
    with open(uninform, 'w') as fh:
        log('writing:', uninform)
        for cluster in fail_clusters:
            fh.write('{}\n'.format(cluster.data[COLUMNS.cluster_id]))
    bedfile = os.path.join(args.output, 'clusters.bed')
    log('writing:', bedfile)
    write_bed_file(bedfile, clusters)

    clusterset_file_prefix = os.path.join(args.output, 'clusterset-')
    log('writing split outputs')
    jobs = []
    i = 0
    header = set()
    while i < len(pass_clusters):
        job = []
        for j in range(0, JOB_SIZE):
            if i >= len(pass_clusters):
                break
            curr = pass_clusters[i]
            row = BreakpointPair.flatten(curr)
            row[COLUMNS.cluster_size] = len(clusters[curr])
            row[COLUMNS.library] = args.library
            row[COLUMNS.protocol] = args.protocol
            job.append(row)
            header.update(row.keys())
            i += 1
        jobs.append(job)
        if len(jobs) == args.max_files:
            while i < len(pass_clusters):
                curr = pass_clusters[i]
                row = BreakpointPair.flatten(curr)
                row[COLUMNS.cluster_size] = len(clusters[curr])
                row[COLUMNS.library] = args.library
                row[COLUMNS.protocol] = args.protocol
                job.append(row)
                header.update(row.keys())
                i += 1
    assert(len(jobs) <= args.max_files)
    log('splitting {} clusters into {} files of size ~{}'.format(len(pass_clusters), len(jobs), JOB_SIZE))
    header = sort_columns(header)

    for i, job in enumerate(jobs):
        # generate an output file
        filename = '{}{}.tab'.format(clusterset_file_prefix, i + 1)
        output_files.append(filename)
        log('writing:', filename, time_stamp=False)
        with open(filename, 'w') as fh:
            fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
            for row in job:
                fh.write('\t'.join([str(row[c]) for c in header]) + '\n')

    return output_files

if __name__ == '__main__':
    args = parse_arguments()
    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)
    main(args)
