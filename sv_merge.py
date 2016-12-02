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

Output File Format
-----------------------------
::

    | column name       | expected value                    | description                                    |
    |-------------------|-----------------------------------|------------------------------------------------|
    | cluster_id        | int                               |                                                |
    | cluster_size      | int > 1                           | the number of individual breakpoint pairs that |
    |                   |                                   | participate in this cluster                    |
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

MAX_JOBS = 20
MIN_EVENTS_PER_JOB = 50

TSV._verbose = False


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
        }
    )
    breakpoints = []

    for row in rows:
        row['tools'] = set([row['tool_version']])
        row['files'] = set([filename])
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
    global MAX_JOBS, MIN_EVENTS_PER_JOB
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
    parser.add_argument('--max-jobs', '-j', default=MAX_JOBS, type=int, dest='MAX_JOBS',
                        help='defines the maximum number of jobs that can be created for the validation step')
    parser.add_argument('--min-events-per-job', '-e', default=MIN_EVENTS_PER_JOB, type=int,
                        help='defines the minimum number of clusters to validate per job', dest='MIN_EVENTS_PER_JOB')
    parser.add_argument('-r', help='radius to use in clustering', default=50, type=int)
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

    MIN_EVENTS_PER_JOB = args.MIN_EVENTS_PER_JOB
    MAX_JOBS = args.MAX_JOBS

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

    for f in args.inputs:
        if not os.path.exists(f):
            print('\nerror: input file {0} does not exist'.format(f))
            parser.print_help()
            exit(1)
    
    
    mkdirp(os.path.join(args.output, 'inputs'))

    breakpoint_pairs = []

    for f in args.inputs:
        print('loading:', f)
        temp = load_input_file(f)
        breakpoint_pairs.extend(temp)

    print('loaded {} breakpoint pairs'.format(len(breakpoint_pairs)))

    # now split by library and protocol
    bpp_by_libprot = {}

    for bpp in breakpoint_pairs:
        lib = bpp.data['libraries']
        protocol = bpp.data['protocol']
        
        d = bpp_by_libprot.setdefault((lib, protocol), {})
        
        if bpp.key not in d:
            d[bpp.key] = bpp
        else:
            d[bpp.key].data['files'].update(bpp.data['files'])
            d[bpp.key].data['tools'].update(bpp.data['tools'])
    
    clusters_by_libprot = {}
    cluster_id_prefix = re.sub(' ', '_', str(datetime.now()))
    cluster_id = 1

    for lib, protocol in bpp_by_libprot:
        bpps = bpp_by_libprot[(lib, protocol)].values()
        if lib not in BAM_FILE_ARGS:
            print(
                'warning: found breakpoints for library', lib,
                ', but bam was not given therefore breakpoints will be ignored'
            )
            continue
        # set up directories
        mkdirp(os.path.join(args.output, 'clustering/{}_{}'.format(lib, protocol)))
        mkdirp(os.path.join(args.output, 'validation/{}_{}'.format(lib, protocol)))
        print('computing clusters for', lib, protocol)
        c = cluster_breakpoint_pairs(bpps, r=args.r, k=args.k)
        clusters_by_libprot[(lib, protocol)] = c
        
        hist = {}
        for cluster, input_pairs in c.items():
            hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
            cluster.data['cluster_id'] = '{}-{}'.format(cluster_id_prefix, cluster_id)
            temp = set()
            for p in input_pairs:
                temp.update(p.data['tools'])
            cluster.data['tools'] = ';'.join(sorted(list(temp)))
            cluster_id += 1

        print('cluster distribution', sorted(hist.items()))
    
    


    header = [
        'cluster_id',
        'cluster_size',
        'break1_chromosome',
        'break1_position_start',
        'break1_position_end',
        'break1_orientation',
        'break1_strand',
        'break2_chromosome',
        'break2_position_start',
        'break2_position_end',
        'break2_orientation',
        'break2_strand',
        'opposing_strands',
        'protocol',
        'library',
        'untemplated_sequence',
        'stranded',
        'tools'
    ]

    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    for lib, protocol in clusters_by_libprot:
        clusters = clusters_by_libprot[(lib, protocol)]
        f = os.path.join(args.output, 'clustering/{}_{}_cluster_assignment.tsv'.format(lib, protocol))
        with open(f, 'w') as fh:
            print('writing:', f)
            h = ['clusters'] + [c for c in header if c not in ['cluster_id', 'cluster_size']] + ['files']
            rows = {}
            fh.write('#' + '\t'.join(h) + '\n')
            for cluster, input_pairs in clusters.items():
                for p in input_pairs:
                    if p in rows:
                        rows[p]['tools'].update(p.data['tools'])
                    else:
                        rows[p] = BreakpointPair.flatten(p)
                    rows[p].setdefault('clusters', set()).add(cluster.data['cluster_id'])

            for row in rows.values():
                row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
                row['tools'] = ';'.join(sorted(list(row['tools'])))
                row['library'] = lib
                row['protocol'] = protocol
                fh.write('\t'.join([str(row[c]) for c in h]) + '\n')



    with open(os.path.join(args.output, 'qsub_all.sh'), 'w') as qsub:
        qsub.write(
            "# script to submit jobs to the cluster\n\n"
            "# array to hold the job ids so we can create the dependency\n"
            "VALIDATION_JOBS=()\n"
        )

        for lib, protocol in clusters_by_libprot:
            clusters = clusters_by_libprot[(lib, protocol)]
            # decide on the number of clusters to validate per job
            JOB_SIZE = MIN_EVENTS_PER_JOB
            if len(clusters) // MIN_EVENTS_PER_JOB > MAX_JOBS - 1:
                JOB_SIZE = len(clusters) // MAX_JOBS

            print('info: splitting', len(clusters), 'clusters into', len(clusters) // JOB_SIZE, 'jobs of size', JOB_SIZE)
            
            index = 0
            fileno = 1
            files = []

            file_prefix = os.path.join(args.output, 'clustering/{}_{}/clusterset'.format(lib, protocol))
            
            bedfile = filename = file_prefix + '.bed'
            print('writing bed file', bedfile)
            write_bed_file(bedfile, clusters)
            
            rows = [c for c in clusters]

            while index < len(rows):
                # generate an output file
                filename = file_prefix + '-{}.tsv'.format(fileno)
                print('writing:', filename)
                with open(filename, 'w') as fh:
                    limit = index + JOB_SIZE
                    fh.write('#' + '\t'.join(header) + '\n')
                    if len(rows) - limit < MIN_EVENTS_PER_JOB or fileno == MAX_JOBS:
                        limit = len(rows)

                    while index < len(rows) and index < limit:
                        row = BreakpointPair.flatten(rows[index])
                        row['cluster_size'] = len(clusters[rows[index]])
                        row['library'] = lib
                        row['protocol'] = protocol
                        fh.write('\t'.join([str(row[c]) for c in header]) + '\n')
                        index += 1
                fileno += 1

                qsub.write(
                    'jid=$(qsub -b sv_validate.py -n {0} -o {1} -b {2} -l {3} -N sv_validate -terse)\n'.format(
                        filename, 
                        os.path.join(args.output, 'validation', '{}_{}'.format(lib, protocol)), 
                        BAM_FILE_ARGS[lib], lib) + 'VALIDATION_JOBS+=("$jid")\n'
                )

if __name__ == '__main__':
    main()
