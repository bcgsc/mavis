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

import TSV
import os
import errno
import argparse
import warnings
import re
from datetime import datetime
from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.interval import Interval
from structural_variant.breakpoint import Breakpoint, BreakpointPair, read_bpp_from_input_file
from structural_variant.cluster import cluster_breakpoint_pairs
from structural_variant.annotate import load_reference_genes
from structural_variant.validate import EvidenceSettings, Evidence
from structural_variant import __version__
from sv_validate import add_evidence_args_to_parser

__prog__ = os.path.basename(os.path.realpath(__file__))


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def load_input_file(filename):
    """
    """
    def nullable_boolean(item):
        try:
            item = TSV.tsv_boolean(item)
        except TypeError:
            if item == '?':
                item = 'null'
            item = TSV.null(item)
        return item

    header, rows = TSV.read_file(
        filename,
        require=[
            'start_chromosome',
            'end_chromosome'
        ],
        rename={
            'tool_version': ['tools']
        },
        split={
            'start_position': '^(?P<start_pos1>\d+)-(?P<start_pos2>\d+)$',
            'end_position': '^(?P<end_pos1>\d+)-(?P<end_pos2>\d+)$',
        },
        cast={
            'start_pos1': int,
            'start_pos2': int,
            'end_pos1': int,
            'end_pos2': int,
            'stranded': TSV.tsv_boolean,
            'opposing_strands': nullable_boolean,
            'untemplated_sequence': TSV.null
        },
        _in={
            'start_strand': STRAND,
            'end_strand': STRAND,
            'start_orientation': ORIENT,
            'end_orientation': ORIENT,
            'protocol': PROTOCOL
        },
        validate={
            'tool_version': '^.+_v?\d+\.\d+\.\d+$',
            'libraries': '^[\w-]+$'
        },
        add={
            'opposing_strands': 'null',
            'start_strand': STRAND.NS,
            'end_strand': STRAND.NS,
            'stranded': False,
            'untemplated_sequence': 'null'
        },
        simplify=True
    )
    breakpoints = []

    for row in rows:
        row['files'] = set([filename])
        row[COLUMNS.tools.name] = set([row[COLUMNS.tools.name]])
        opposing_strands = [row['opposing_strands']]
        if row['opposing_strands'] is None and STRAND.NS in [row['start_strand'], row['end_strand']]:
            opposing_strands = [True, False]

        for opp in opposing_strands:
            b1 = Breakpoint(
                row['start_chromosome'],
                row['start_pos1'],
                row['start_pos2'],
                orient=row['start_orientation'],
                strand=row['start_strand'])
            b2 = Breakpoint(
                row['end_chromosome'],
                row['end_pos1'],
                row['end_pos2'],
                orient=row['end_orientation'],
                strand=row['end_strand'])
            try:
                bpp = BreakpointPair(
                    b1,
                    b2,
                    opposing_strands=opp,
                    untemplated_sequence=row['untemplated_sequence'],
                    stranded=row['stranded'],
                    data=row
                )
                if not row['stranded']:
                    bpp.break1.strand = STRAND.NS
                    bpp.break2.strand = STRAND.NS
                breakpoints.append(bpp)
            except InvalidRearrangement as e:
                warnings.warn(str(e) + '; reading: {}'.format(filename))
    return breakpoints


def mkdirp(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


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


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number')
    parser.add_argument(
        '-f', '--overwrite', action='store_true', default=False,
        help='set flag to overwrite existing reviewed files')
    parser.add_argument(
        '-o', '--output', help='path to the output directory', required=True)
    parser.add_argument(
        '-n', '--inputs', help='path to the input files', required=True, action='append')
    parser.add_argument(
        '-b', '--bamfile', metavar=('<library_name>', '</path/to/bam/file>'), nargs=2,
        help='specify a bam file for a given library', action='append', required=True)
    g = parser.add_argument_group('qsub job options')
    g.add_argument(
        '--max-jobs', '-j', default=100, type=int, dest='MAX_JOBS',
        help='defines the maximum number of jobs that can be created for the validation step')
    g.add_argument(
        '--min-events-per-job', '-e', default=50, type=int,
        help='defines the minimum number of clusters to validate per job', dest='MIN_EVENTS_PER_JOB')
    parser.add_argument(
        '-r', help='radius to use in clustering', default=20, type=int)
    g.add_argument(
        '-q', '--queue', default='transabyss.q',
        help='queue to submit validation jobs to, only affects the qsub script created')
    parser.add_argument(
        '-k', help='parameter used for computing cliques, smaller is faster, above 20 will be slow',
        default=15, type=int)
    g = parser.add_argument_group('filter arguments')
    g.add_argument(
        '--no_filter', default=False, help='If flag is given the clusters will not be filtered '
        'based on lack of annotation')
    g.add_argument(
        '--filter_proximity', '-p', type=int, default=5000, help='maximum distance to look for annotations'
        'from evidence window')
    g.add_argument(
        '--annotations', '-a',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_annotations_20170203.json',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.')
    g = parser.add_argument_group('evidence settings')
    add_evidence_args_to_parser(g)
    args = parser.parse_args()

    if args.MIN_EVENTS_PER_JOB < 1:
        print('\nerror: MIN_EVENTS_PER_JOB cannot be less than 1')
        parser.print_help()
        exit(1)

    if args.MAX_JOBS < 1:
        print('\nerror: MAX_JOBS cannot be less than 1')
        parser.print_help()
        exit(1)

    if os.path.exists(args.output) and not args.overwrite:
        print(
            '\nerror: output directory {0} already exists. please use the --overwrite option'.format(args.output))
        parser.print_help()
        exit(1)

    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)
    BAM_FILE_ARGS = {}

    for lib, bam in args.bamfile:
        if lib in BAM_FILE_ARGS:
            print('\nerror: library can only specify a single bam')
            parser.print_help()
            exit(1)
        BAM_FILE_ARGS[lib] = bam

    args.output = os.path.abspath(args.output)

    for f in args.inputs:
        if not os.path.exists(f):
            print('\nerror: input file {0} does not exist'.format(f))
            parser.print_help()
            exit(1)

    mkdirp(os.path.join(args.output, 'inputs'))

    breakpoint_pairs = []

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
            }
        )
        for bpp in bpps:
            bpp.data[COLUMNS.tools] = set(';'.split(bpp.data[COLUMNS.tools]))
        breakpoint_pairs.extend(bpps)

    log('loaded {} breakpoint pairs'.format(len(breakpoint_pairs)))

    if not args.no_filter:
        log('loading:', args.annotations)
        REFERENCE_GENES = load_reference_genes(args.annotations, verbose=False)

    # now split by library and protocol
    bpp_by_libprot = {}

    for bpp in breakpoint_pairs:
        lib = bpp.data[COLUMNS.library]
        protocol = bpp.data[COLUMNS.protocol.name]

        d = bpp_by_libprot.setdefault((lib, protocol), {})

        if bpp not in d:
            d[bpp] = bpp
        else:
            d[bpp].data['files'].update(bpp.data['files'])
            d[bpp].data[COLUMNS.tools].update(bpp.data[COLUMNS.tools])

    clusters_by_libprot = {}
    cluster_id_prefix = re.sub(' ', '_', str(datetime.now()))
    cluster_id = 1

    for lib, protocol in bpp_by_libprot:
        bpps = bpp_by_libprot[(lib, protocol)].values()
        if lib not in BAM_FILE_ARGS:
            log(
                'warning: found breakpoints for library',
                lib, ', but bam was not given therefore breakpoints will be ignored')
            continue
        # set up directories
        mkdirp(os.path.join(args.output, 'clustering/{}_{}'.format(lib, protocol)))
        mkdirp(os.path.join(args.output, 'validation/{}_{}/log'.format(lib, protocol)))
        log('computing clusters for', lib, protocol)
        c = cluster_breakpoint_pairs(bpps, r=args.r, k=args.k)
        clusters_by_libprot[(lib, protocol)] = c

        hist = {}
        for cluster, input_pairs in c.items():
            hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
            cluster.data[COLUMNS.cluster_id.name] = 'cluster_{}-{}'.format(cluster_id_prefix, cluster_id)
            temp = set()
            for p in input_pairs:
                temp.update(p.data[COLUMNS.tools.name])
            cluster.data[COLUMNS.tools.name] = ';'.join(sorted(list(temp)))
            cluster_id += 1
        log('input', len(bpps), 'breakpoint pairs', time_stamp=False)
        log('computed', len(c), 'clusters', time_stamp=False)
        log('cluster distribution', sorted(hist.items()), time_stamp=False)

    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    for lib, protocol in clusters_by_libprot:
        clusters = clusters_by_libprot[(lib, protocol)]
        f = os.path.join(args.output, 'clustering/{}_{}/cluster_assignment.tab'.format(lib, protocol))
        with open(f, 'w') as fh:
            header = set()
            log('writing:', f)
            rows = {}

            for cluster, input_pairs in clusters.items():
                for p in input_pairs:
                    if p in rows:
                        rows[p][COLUMNS.tools.name].update(p.data[COLUMNS.tools.name])
                    else:
                        rows[p] = BreakpointPair.flatten(p)
                    rows[p].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id.name])
            for row in rows.values():
                row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
                row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools.name])))
                row[COLUMNS.library] = lib
                row[COLUMNS.protocol] = protocol
                header.update(row.keys())
            header = sort_columns(header)
            fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
            for row in rows.values():
                fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

    settings = EvidenceSettings.parse_args(args)

    for lib, protocol in clusters_by_libprot:
        clusters = clusters_by_libprot[(lib, protocol)]
        # decide on the number of clusters to validate per job
        pass_clusters = []
        fail_clusters = []

        for cluster in clusters:
            # don't need to generate transcriptome windows b/c will default to genome if not in a gene anyway
            w1 = Evidence.generate_window(
                cluster.break1,
                read_length=settings.read_length,
                median_insert_size=settings.median_insert_size,
                call_error=settings.call_error,
                stdev_isize=settings.stdev_isize,
                stdev_count_abnormal=settings.stdev_count_abnormal
            )
            w2 = Evidence.generate_window(
                cluster.break2,
                read_length=settings.read_length,
                median_insert_size=settings.median_insert_size,
                call_error=settings.call_error,
                stdev_isize=settings.stdev_isize,
                stdev_count_abnormal=settings.stdev_count_abnormal
            )
            if args.no_filter:
                pass_clusters.append(cluster)
            else:
                # loop over the annotations
                overlaps_gene = False
                w1 = Interval(w1.start - args.filter_proximity, w1.end + args.filter_proximity)
                w2 = Interval(w2.start - args.filter_proximity, w2.end + args.filter_proximity)
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

        log('filtered', len(fail_clusters), 'clusters as not informative')

        JOB_SIZE = args.MIN_EVENTS_PER_JOB
        if len(pass_clusters) // args.MIN_EVENTS_PER_JOB > args.MAX_JOBS - 1:
            JOB_SIZE = len(pass_clusters) // args.MAX_JOBS

        log('info: splitting {} clusters into {} jobs of size {}'.format(
            len(pass_clusters), len(pass_clusters) // JOB_SIZE, JOB_SIZE))

        index = 0
        fileno_start = 1
        fileno = fileno_start
        parent_dir = os.path.join(args.output, 'clustering/{}_{}'.format(lib, protocol))
        uninform = os.path.join(parent_dir, 'uninformative_clusters.txt')

        with open(uninform, 'w') as fh:
            log('writing:', uninform)
            for cluster in fail_clusters:
                fh.write('{}\n'.format(cluster.data[COLUMNS.cluster_id]))
        bedfile = filename = '{}/clusters.bed'.format(parent_dir)
        log('writing:', bedfile)
        write_bed_file(bedfile, clusters)

        rows = [c for c in clusters]

        clusterset_file_prefix = parent_dir + '/clusterset-'
        log('writing split outputs')
        while index < len(rows):
            # generate an output file
            filename = '{}{}.tab'.format(clusterset_file_prefix, fileno)
            log('writing:', filename, time_stamp=False)
            with open(filename, 'w') as fh:
                limit = index + JOB_SIZE
                header = None
                if len(rows) - limit < args.MIN_EVENTS_PER_JOB or fileno == args.MAX_JOBS:
                    limit = len(rows)

                while index < len(rows) and index < limit:
                    row = BreakpointPair.flatten(rows[index])
                    row[COLUMNS.cluster_size.name] = len(clusters[rows[index]])
                    row[COLUMNS.library.name] = lib
                    row[COLUMNS.protocol.name] = protocol
                    if not header:
                        header = sort_columns(row.keys())
                        fh.write('#' + '\t'.join(header) + '\n')
                    fh.write('\t'.join([str(row[c]) for c in header]) + '\n')
                    index += 1
            fileno += 1

        # create the qsub script as well
        output_folder = '{}/validation/{}_{}'.format(args.output, lib, protocol)
        qsub_file = '{}/qsub.sh'.format(output_folder)
        with open(qsub_file, 'w') as fh:
            log('writing:', qsub_file)
            fh.write('#!/bin/sh\n')
            fh.write('#$ -t {}-{}\n'.format(fileno_start, fileno - 1))  # array job
            fh.write('#$ -V\n')  # copy environment variables
            fh.write('#$ -N svmV_{}\n'.format(lib[-5:]))
            fh.write('#$ -j y\n')
            fh.write('#$ -q {}\n'.format(args.queue))
            fh.write('#$ -o {}/log/\n'.format(output_folder))
            fh.write('#$ -l mem_free=12G,mem_token=12G,h_vmem=12G\n')
            fh.write('echo "Starting job: $SGE_TASK_ID"\n')
            fh.write('python sv_validate.py -n {}$SGE_TASK_ID.tab \\\n\t-o {} \\\n\t-b {} \\\n\t-l {} \\\n'.format(
                clusterset_file_prefix,
                output_folder,
                BAM_FILE_ARGS[lib],
                lib
            ))
            temp = EvidenceSettings()
            opt = []
            for param, val in args.__dict__.items():
                if hasattr(temp, param):
                    opt.append('--{} {}'.format(param, val))
            fh.write('\t{}\n'.format(' \\\n\t'.join(opt)))
            fh.write('echo "Job complete: $SGE_TASK_ID"\n')

if __name__ == '__main__':
    main()
