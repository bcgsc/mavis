import re
import TSV
import os
import errno
import argparse
import warnings
import datetime
from structural_variant.constants import *
from structural_variant.align import *
from structural_variant.error import *
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.cluster import cluster_breakpoints
from structural_variant import __version__

__prog__ = os.path.basename(os.path.realpath(__file__))

MAX_JOBS = 50
MIN_EVENTS_PER_JOB = 10

TSV._verbose = False


def load_input_file(filename):
    """
    returns sets of breakpoint pairs keyed by (library, protocol, tool)
    """
    header, rows = TSV.read_file(
        filename,
        retain=[
            'start_chromosome',
            'end_chromosome',
            'start_orientation',
            'end_orientation',
            'start_strand',
            'end_strand',
            'protocol',
            'tool_version'],
        split={
            'start_position': ('^(\d+)-(\d+)$', ['start_pos1', 'start_pos2']),
            'end_position': ('^(\d+)-(\d+)$', ['end_pos1', 'end_pos2']),
        },
        cast={'start_pos1': 'int', 'start_pos2': 'int',
              'end_pos1': 'int', 'end_pos2': 'int'},
        validate={
            'start_orientation': '^{0}$'.format('|'.join([re.escape(x) for x in ORIENT.values()])),
            'end_orientation': '^{0}$'.format('|'.join([re.escape(x) for x in ORIENT.values()])),
            'start_strand': '^{0}$'.format('|'.join([re.escape(x) for x in STRAND.values()])),
            'end_strand': '^{0}$'.format('|'.join([re.escape(x) for x in STRAND.values()])),
            'tool_version': '^.+_v\d+\.\d+\.\d+$',
            'protocol': '^(genome|transcriptome)$',
            'libraries': '^[\w-]+(;[\w-]+)*$'
        }
    )
    breakpoints = {}

    for row in rows:
        for lib, prot, tool in itertools.product(
                row['libraries'].split(';'), row['protocol'].split(';'), row['tool_version'].split(';')):
            label = {'library': lib, 'protocol': prot,
                     'tool_version': tool, 'input_file': filename}
            key = tuple([k[1] for k in sorted(label.items())])
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
                if key not in breakpoints:
                    breakpoints[key] = set()
                bpp = BreakpointPair(b1, b2, flags=[] if FLAGS.LQ not in row.get(
                    'filters', '').split(';') else [FLAGS.LQ])
                bpp.label = label
                breakpoints[key].add(bpp)
            except InvalidRearrangement as e:
                warnings.warn(str(e))
    return breakpoints


def mkdirp(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


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

    for l, b in args.bamfile:
        mkdirp(os.path.join(args.output, 'clustering', l))
        mkdirp(os.path.join(args.output, 'validation', l))
    mkdirp(os.path.join(args.output, 'log'))

    pairs = {}

    for f in args.inputs:
        print('loading:', f)
        temp = load_input_file(f)
        pairs.update(temp)

    bpp_by_library = {}

    for key, bpp_set in pairs.items():
        filename, libname, protocol, tool = key
        if libname not in BAM_FILE_ARGS:
            raise UserWarning('error: input file library specified does not have a corresponding bam file. '
                              'must be specified', libname)
        print('loaded', len(bpp_set), 'breakpoint pairs for', key)
        if (libname, protocol) not in bpp_by_library:
            bpp_by_library[(libname, protocol)] = []
        bpp_by_library[(libname, protocol)].extend(list(bpp_set))

    cluster_id = 0
    cluster_rows_by_lib = {}

    for l in BAM_FILE_ARGS:
        cluster_rows_by_lib[l] = []

    temp = os.path.join(args.output, 'log/cluster_assignment.tsv')
    with open(temp, 'w') as fh:
        assignment_header = [
            'cluster_ids',
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
            'protocol',
            'library',
            'tool_version',
            'input_file',
        ]
        fh.write('## {0} v{1} {2}\n'.format(
            __prog__, __version__, datetime.datetime.now()))
        fh.write(
            '## this file details the inputs and the cluster ids they were assigned to\n')
        fh.write('#' + '\t'.join(assignment_header) + '\n')
        for key, bpp_list in bpp_by_library.items():
            libname, protocol = key
            print('for', key, 'there are', len(
                bpp_list), 'input breakpoint pairs')

            clusters = cluster_breakpoints(bpp_list, r=20, k=15)
            initial_count = len(clusters)

            # filter out the low quality clusters
            temp = list(clusters.keys())
            for cluster in temp:  # set the ids of the clusters
                cluster.label = cluster_id
                cluster_id += 1
                if len(clusters[cluster]) == 1:
                    if FLAGS.LQ in list(clusters[cluster])[0].flags:
                        del clusters[cluster]
            print('after clustering there are', len(clusters),
                  'quality breakpoint pairs and', initial_count, 'total pairs')

            # track the inputs to their clusters
            for bpp in bpp_list:
                temp = []
                for c, cset in clusters.items():
                    if bpp in cset:
                        temp.append(c.label)
                row = {
                    'cluster_ids': ';'.join(sorted([str(k) for k in temp])),
                    'break1_chromosome': bpp.break1.chr,
                    'break1_position_start': bpp.break1.start,
                    'break1_position_end': bpp.break1.end,
                    'break1_orientation': bpp.break1.orient,
                    'break1_strand': bpp.break1.strand,
                    'break2_chromosome': bpp.break2.chr,
                    'break2_position_start': bpp.break2.start,
                    'break2_position_end': bpp.break2.end,
                    'break2_orientation': bpp.break2.orient,
                    'break2_strand': bpp.break2.strand
                }
                for tag, value in bpp.label.items():
                    if tag in row:
                        warnings.warn(
                            'label tag {0}={1} is dropped b/c conflicts with header'.format(tag. value))
                    else:
                        row[tag] = value
                fh.write('\t'.join(str(row[k])
                                   for k in assignment_header) + '\n')

            for c in clusters:
                tools = ';'.join(
                    sorted(list(set([k.label['tool_version'] for k in clusters[c]]))))
                row = {
                    'cluster_id': c.label,
                    'cluster_size': len(clusters[c]),
                    'break1_chromosome': c.break1.chr,
                    'break1_position_start': c.break1.start,
                    'break1_position_end': c.break1.end,
                    'break1_orientation': c.break1.orient,
                    'break1_strand': c.break1.strand,
                    'break2_chromosome': c.break2.chr,
                    'break2_position_start': c.break2.start,
                    'break2_position_end': c.break2.end,
                    'break2_orientation': c.break2.orient,
                    'break2_strand': c.break2.strand,
                    'opposing_strands': c.opposing_strands,
                    'protocol': protocol,
                    'library': libname,
                    'tools': tools
                }
                cluster_rows_by_lib[libname].append(row)
                # add the reciprocal for translocations and inversions
                if len(set([SVTYPE.INV, SVTYPE.ITRANS, SVTYPE.TRANS]) & set(BreakpointPair.classify(c))) > 0:
                    reciprocal = {}
                    reciprocal.update(row)
                    if reciprocal['break1_orientation'] == ORIENT.LEFT:
                        reciprocal['break1_orientation'] = ORIENT.RIGHT
                    elif reciprocal['break1_orientation'] == ORIENT.RIGHT:
                        reciprocal['break1_orientation'] = ORIENT.LEFT

                    if reciprocal['break2_orientation'] == ORIENT.LEFT:
                        reciprocal['break2_orientation'] = ORIENT.RIGHT
                    elif reciprocal['break2_orientation'] == ORIENT.RIGHT:
                        reciprocal['break2_orientation'] = ORIENT.LEFT
                    cluster_rows_by_lib[libname].append(reciprocal)

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
        'tools'
    ]

    qsub = open(os.path.join(args.output, 'qsub_all.sh'), 'w')
    qsub.write(
        "# script to submit jobs to the cluster\n\n"
        "# array to hold the job ids so we can create the dependency\n"
        "VALIDATION_JOBS=()\n"
    )

    for lib, cluster_rows in cluster_rows_by_lib.items():
        # decide on the number of clusters to validate per job
        JOB_SIZE = MIN_EVENTS_PER_JOB
        if len(cluster_rows) // MIN_EVENTS_PER_JOB > MAX_JOBS - 1:
            JOB_SIZE = len(cluster_rows) // MAX_JOBS

        print('splitting', len(cluster_rows), 'clusters into', len(
            cluster_rows) // JOB_SIZE, 'jobs of size', JOB_SIZE)
        # split the clusters by chromosomes
        cluster_rows.sort(key=lambda x: (
            x['break1_chromosome'], x['break2_chromosome']))

        index = 0
        fileno = 1
        files = []
        while index < len(cluster_rows):
            # generate an output file
            filename = os.path.abspath(os.path.join(
                args.output, 'clustering/{0}/{0}-clusterset-{1}.tsv'.format(lib, fileno)))
            files.append(filename)
            row_subset = []
            print('writing:', filename)
            with open(filename, 'w') as fh:
                temp = index
                limit = index + JOB_SIZE
                fh.write('#' + '\t'.join(header) + '\n')
                if len(cluster_rows) - limit < MIN_EVENTS_PER_JOB or fileno == MAX_JOBS:
                    limit = len(cluster_rows)
                while temp < len(cluster_rows) and temp < limit:
                    row = [str(cluster_rows[temp][k]) for k in header]
                    row_subset.append(cluster_rows[temp])
                    fh.write('\t'.join(row) + '\n')
                    temp += 1
            bedfile = os.path.abspath(os.path.join(
                args.output, 'clustering/{0}/{0}-clusterset-{1}.bed'.format(lib, fileno)))
            with open(bedfile, 'w') as fh:
                print('writing:', bedfile)
                for row in row_subset:
                    if row['break1_chromosome'] == row['break2_chromosome']:
                        fh.write('{0}\t{1}\t{2}\t{0}:{1}{3}{4}-{2}{5}{6}\n'.format(
                            row['break2_chromosome'],
                            row['break1_position_start'],
                            row['break2_position_end'],
                            row['break1_orientation'],
                            row['break1_strand'],
                            row['break2_orientation'],
                            row['break2_strand']
                        ))
                    else:
                        fh.write('{0}\t{1}\t{2}\t{0}:{1}{3}{4}\n'.format(
                            row['break1_chromosome'],
                            row['break1_position_start'],
                            row['break1_position_end'],
                            row['break1_orientation'],
                            row['break1_strand']
                        ))
                        fh.write('{0}\t{1}\t{2}\t{0}:{1}{3}{4}\n'.format(
                            row['break2_chromosome'],
                            row['break2_position_start'],
                            row['break2_position_end'],
                            row['break2_orientation'],
                            row['break2_strand']
                        ))
            index = limit
            fileno += 1

            output_file = os.path.abspath(os.path.join(args.output, 'validation', lib))
            qsub.write(
                'jid=$(qsub -b sv_validate.py -n {0} -o {1} -b {2} -l {3} -N sv_validate -terse)\n'.format(
                    filename, output_file, BAM_FILE_ARGS[lib], lib) + 'VALIDATION_JOBS+=("$jid")\n'
            )
        # now create the qsub scripts
    # merge the bam files
    qsub.close()

main()
