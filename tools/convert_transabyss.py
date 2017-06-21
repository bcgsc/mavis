"""
Script for converting Trans-ABySS output file into the MAVIS accepted input format
"""

import argparse
import os
import sys
import time
import TSV
import warnings

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, STRAND, PROTOCOL
from mavis.error import InvalidRearrangement
from mavis import __version__

__prog__ = os.path.basename(os.path.realpath(__file__))

TSV._verbose = True

SVTYPES = {'ins': SVTYPE.INS,
           'del': SVTYPE.DEL,
           'dup': SVTYPE.DUP,
           'ITD': SVTYPE.DUP,
           'insertion': SVTYPE.INS,
           'deletion': SVTYPE.DEL,
           'duplication': SVTYPE.DUP,
           'translocation': SVTYPE.TRANS,
           'inversion': SVTYPE.INV
           }


def main():
    parser = argparse.ArgumentParser(
        description='Convert a TA output file to the MAVIS input format, takes both indels and fusions',
        add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
    required.add_argument('-o', '--output', help='path to the output file', required=True)
    required.add_argument('-p', '--protocol', choices=[PROTOCOL.GENOME, PROTOCOL.TRANS], required=True)
    required.add_argument('-l', '--library', default=None,
                          help="the library id of that was used as input")

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='outputs the version number')
    optional.add_argument('--stranded', action='store_true', default=False)
    optional.add_argument('--tool-version', help='the version of trans abyss that was used in the analysis',
                          default='1.4.10')
    optional.add_argument('--no-filter', action='store_true', default=False,
                          help='turn off filtering of events that are in the "MT" and "GL" chromosomes')

    # /projects/POG/POG_data/POG098/wgs/GV2/POG098_POG098-OCT-1-unique-14-filters/POG098-OCT-1_genome_fusions_concat.tsv
    warnings.warn('currently assuming that trans-abyss is calling the strand exactly opposite and swapping them')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print('error: input file {0} does not exist'.format(args.input))
        sys.exit(1)
    print('reading:', args.input)
    file_type = None
    try:
        header, rows = TSV.read_file(
            args.input,
            require=['id'],
            rename={'rearrangement': [COLUMNS.event_type]},
            split={
                'breakpoint': '^(?P<chr1>[^:]+):(?P<pos1>\d+)\|(?P<chr2>[^:]+):(?P<pos2>\d+)$',
                'orientations': '^(?P<or1>[RL]),(?P<or2>[RL])$',
                'strands': '^(?P<strand1>[\+-]),(?P<strand2>[\+-])$'
            },
            cast={
                'pos1': int,
                'pos2': int,
                'strand1': lambda x: STRAND.NEG if x == STRAND.POS else STRAND.POS,
                'strand2': lambda x: STRAND.NEG if x == STRAND.POS else STRAND.POS
            },
            strict=False,
            in_={
                'strand1': STRAND,
                'strand2': STRAND,
                'or1': ORIENT,
                'or2': ORIENT
            }
        )
        file_type = 'fusions'
    except KeyError:
        pass
    if not file_type:
        try:
            header, rows = TSV.read_file(
                args.input,
                require=['id', 'type'],
                cast={
                    'chr_start': int,
                    'chr_end': int
                }
            )
            file_type = 'indels'
        except KeyError:
            print('error: input file {0} does not have the expected columns'.format(args.input))
            sys.exit(1)

    bpps = set()
    filter_count = 0
    for row in rows:
        if file_type == 'indels':
            strand1 = strand2 = row['ctg_strand']
            chr1 = chr2 = row['chr']
            event_type = SVTYPES[row['type']]

            if event_type == SVTYPE.DUP:
                orient1, orient2 = (ORIENT.RIGHT, ORIENT.LEFT)
            elif event_type == SVTYPE.INS or event_type == SVTYPE.DEL:
                orient1, orient2 = (ORIENT.LEFT, ORIENT.RIGHT)
                if row['chr_start'] == row['chr_end']:
                    row['chr_end'] += 1
            else:
                print("ERROR: Found an unexpected event type in the indel file {}".format(event_type))
                sys.exit(1)

            pos1, pos2 = (row['chr_start'], row['chr_end'])

        else:
            strand1, strand2 = (row['strand1'], row['strand2'])
            orient1, orient2 = (row['or1'], row['or2'])
            pos1, pos2 = (row['pos1'], row['pos2'])
            chr1, chr2 = (row['chr1'], row['chr2'])
            event_type = SVTYPES[row['rearrangement']]

        if not args.no_filter and ("GL" in chr1 + chr2 or 'M' in chr1 + chr2):
            filter_count += 1
            continue

        bpp = BreakpointPair(
            Breakpoint(chr1, pos1, strand=strand1, orient=orient1),
            Breakpoint(chr2, pos2, strand=strand2, orient=orient2),
            data={
                COLUMNS.library: args.library,
                COLUMNS.protocol: args.protocol,
                COLUMNS.tools: 'TransABySS_v{0}'.format(args.tool_version),
                COLUMNS.event_type: event_type
            },
            stranded=args.stranded
        )
        if len(set([bpp.break1.orient, bpp.break2.orient, ORIENT.NS])) == 2:
            if not bpp.interchromosomal:
                bpp.data[COLUMNS.event_type] = SVTYPE.INV
            else:
                bpp.data[COLUMNS.event_type] = SVTYPE.ITRANS
        if bpp.data[COLUMNS.event_type] not in BreakpointPair.classify(bpp):
            print(bpp.break1, bpp.break2)
            print(row)
            raise InvalidRearrangement('expected', BreakpointPair.classify(bpp), 'found',
                                       bpp.data[COLUMNS.event_type])
        bpps.add(bpp)

    header = set()
    rows = []
    for bpp in bpps:
        f = bpp.flatten()
        header.update(f.keys())
        rows.append(f)

    header = sort_columns(header)

    if filter_count > 0:
        print('{0} events have been filtered'.format(filter_count))
    with open(args.output, 'w') as fh:
        print('writing:', args.output)
        fh.write('## {1} v{0}\n'.format(__version__, __prog__))
        fh.write('## inputs {0}\n'.format(" ".join(sys.argv)))
        # fh.write('## input: {0}\n'.format(args.input))
        # fh.write('## output: {0}\n'.format(args.output))
        # fh.write('## library: {0}\n'.format(args.library))
        # fh.write('## protocol: {0}\n'.format(args.protocol))
        fh.write('## file generated on {0}\n'.format(time.strftime('%B %d, %Y')))
        fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row[c]) for c in header]) + '\n')
        print("Wrote {0} gene {1} events to {2}".format(len(rows), file_type, args.output))

if __name__ == '__main__':
    main()
