"""
this is the script used for merging a set of input structural variant calls
into clusters

Input File Expected Format
--------------------------
::

    | column name       | expected value                    | description                                    |
    |-------------------|-----------------------------------|------------------------------------------------|
    | start_position    | <int>-<int>                       |                                                |
    | start_strand      | <+,-,?>                           | the reference strand aligned to                |
    | start_orientation | <L,R,?>                           |                                                |
    | end_chromosome    | <1-22,X,Y,MT>                     |                                                |
    | end_position      | <int>-<int>                       |                                                |
    | end_strand        | <+,-,?>                           | the reference strand aligned to                |
    | end_orientation   | <L,R,?>                           |                                                |
    | protocol          | <genome or transcriptome>         |                                                |
    | library           | library id                        |                                                |
    | tool_version      | <tool name>_<tool version number> |                                                |
    | opposing_strand   | <True,False,?>                    |                                                |

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
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.cluster import cluster_breakpoint_pairs
from structural_variant import __version__

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
            item = TSV.bool(item)
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
            'stranded': TSV.bool,
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
            'tool_version': '^.+_v\d+\.\d+\.\d+$',
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
    args = parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                        help='Outputs the version number')
    parser.add_argument('-f', '--overwrite', action='store_true', default=False,
                        help='set flag to overwrite existing reviewed files')
    parser.add_argument(
        '-o', '--output', help='path to the output directory', required=True)
    parser.add_argument(
        '-n', '--inputs', help='path to the input files', required=True, action='append')
    parser.add_argument('-b', '--bamfile', metavar=('<library_name>', '</path/to/bam/file>'), nargs=2,
                        help='specify a bam file for a given library', action='append', required=True)
    parser.add_argument('--max-jobs', '-j', default=100, type=int, dest='MAX_JOBS',
                        help='defines the maximum number of jobs that can be created for the validation step')
    parser.add_argument('--min-events-per-job', '-e', default=50, type=int,
                        help='defines the minimum number of clusters to validate per job', dest='MIN_EVENTS_PER_JOB')
    parser.add_argument('-r', help='radius to use in clustering', default=50, type=int)
    parser.add_argument('-q', '--queue', default='transabyss.q',
                        help='queue to submit validation jobs to, only affects the qsub script created')
    parser.add_argument(
        '-k', help='parameter used for computing cliques, smaller is faster, above 20 will be slow',
        default=15, type=int)
    args = parser.parse_args()

    if args.MIN_EVENTS_PER_JOB < 1:
        print('\nerror: MIN_EVENTS_PER_JOB cannot be less than 1')
        parser.print_help()
        exit(1)

    if args.MAX_JOBS < 1:
        print('\nerror: MAX_JOBS cannot be less than 1')
        parser.print_help()
        exit(1)

    BAM_FILE_ARGS = {}

    for lib, bam in args.bamfile:
        if lib in BAM_FILE_ARGS:
            print('\nerror: library can only specify a single bam')
            parser.print_help()
            exit(1)
        BAM_FILE_ARGS[lib] = bam

    if os.path.exists(args.output) and not args.overwrite:
        print(
            '\nerror: output directory {0} already exists. please use the --overwrite option'.format(args.output))
        parser.print_help()
        exit(1)

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
        temp = load_input_file(f)
        breakpoint_pairs.extend(temp)

    log('loaded {} breakpoint pairs'.format(len(breakpoint_pairs)))

    # now split by library and protocol
    bpp_by_libprot = {}

    for bpp in breakpoint_pairs:
        lib = bpp.data['libraries']
        protocol = bpp.data[COLUMNS.protocol.name]

        d = bpp_by_libprot.setdefault((lib, protocol), {})

        if bpp not in d:
            d[bpp] = bpp
        else:
            d[bpp].data['files'].update(bpp.data['files'])
            d[bpp].data[COLUMNS.tools.name].update(bpp.data[COLUMNS.tools.name])

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

        log('cluster distribution', sorted(hist.items()))

    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    for lib, protocol in clusters_by_libprot:
        clusters = clusters_by_libprot[(lib, protocol)]
        f = os.path.join(args.output, 'clustering/{}_{}/cluster_assignment.tsv'.format(lib, protocol))
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
                row[COLUMNS.tools.name] = ';'.join(sorted(list(row[COLUMNS.tools.name])))
                row[COLUMNS.library.name] = lib
                row[COLUMNS.protocol.name] = protocol
                header.update(row.keys())
            header = sort_columns(header)
            fh.write('#' + '\t'.join(header) + '\n')
            for row in rows.values():
                fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')

    for lib, protocol in clusters_by_libprot:
        clusters = clusters_by_libprot[(lib, protocol)]
        # decide on the number of clusters to validate per job
        JOB_SIZE = args.MIN_EVENTS_PER_JOB
        if len(clusters) // args.MIN_EVENTS_PER_JOB > args.MAX_JOBS - 1:
            JOB_SIZE = len(clusters) // args.MAX_JOBS

        log('info: splitting {} clusters into {} jobs of size {}'.format(
            len(clusters), len(clusters) // JOB_SIZE, JOB_SIZE))

        index = 0
        fileno_start = 1
        fileno = fileno_start
        parent_dir = os.path.join(args.output, 'clustering/{}_{}'.format(lib, protocol))

        bedfile = filename = '{}/clusters.bed'.format(parent_dir)
        log('writing bed file:', bedfile)
        write_bed_file(bedfile, clusters)

        rows = [c for c in clusters]

        clusterset_file_prefix = parent_dir + '/clusterset-'
        while index < len(rows):
            # generate an output file
            filename = '{}{}'.format(clusterset_file_prefix, fileno)
            log('writing:', filename)
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
        qsub_file = '{}/qsub.sh'.format(parent_dir)
        with open(qsub_file, 'w') as fh:
            log('writing:', qsub_file)
            fh.write('#!/bin/sh\n')
            fh.write('#$ -t {}-{}\n'.format(fileno_start, fileno - 1))  # array job
            fh.write('#$ -V\n')  # copy environment variables
            fh.write('#$ -N svmV_{}\n'.format(lib[-5:]))
            fh.write('#$ -j y\n')
            fh.write('#$ -q {}\n'.format(args.queue))
            output_folder = '{}/validation/{}_{}'.format(args.output, lib, protocol)
            fh.write('#$ -o {}/log/\n'.format(output_folder))
            fh.write('#$ -l mem_free=12G,mem_token=12G,h_vmem=12G\n')
            fh.write('echo "Starting job: $SGE_TASK_ID"\n')
            fh.write('python sv_validate.py -n {}$SGE_TASK_ID -o {} -b {} -l {}\n'.format(
                clusterset_file_prefix,
                output_folder,
                BAM_FILE_ARGS[lib],
                lib
            ))
            fh.write('echo "Job complete: $SGE_TASK_ID"\n')

if __name__ == '__main__':
    main()
