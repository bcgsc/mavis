"""
script for converting Trans-ABySS indel output files into the MAVIS accepted input format
"""
import TSV
import argparse
import warnings
import os
import sys
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.error import InvalidRearrangement

__version__ = '0.0.1'
__prog__ = os.path.basename(os.path.realpath(__file__))

TSV._verbose = True

# chr     chr_start       chr_end      len      ref     alt type
#types ins, del, dup, ITD


SVTYPES = {'ins': SVTYPE.INS,
           'del': SVTYPE.DEL,
           'dup': SVTYPE.DUP,
           'ITD': SVTYPE.DUP
           }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number')
    parser.add_argument(
        '-f', '--overwrite', action='store_true', default=False,
        help='set flag to overwrite existing reviewed files')
    parser.add_argument('-o', '--output', help='path to the output file', required=True)
    parser.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
    parser.add_argument('-p', '--protocol', choices=['genome', 'transcriptome'], required=True)
    parser.add_argument('-l', '--library_id', required=True)
    parser.add_argument('--stranded', action='store_true', default=False)
#    warnings.warn('currently assuming that trans-abyss is calling the strand exactly opposite and swapping them')
    args = parser.parse_args()

    if os.path.exists(args.output) and not args.overwrite:
        print('error: output file {0} already exists. please use the --overwrite option'.format(args.output))
        parser.print_help()
        exit()
    if not os.path.exists(args.input):
        print('error: input file {0} does not exist'.format(args.input))
        parser.print_help()
        exit()
    print('reading:', args.input)

    # /projects/POG/POG_data/POG098/wgs/GV2/POG098_POG098-OCT-1-unique-14-filters/POG098-OCT-1_genome_fusions_concat.tsv

    header, rows = TSV.read_file(
        args.input,
        require=['id'],
        cast={
            'chr_start': int,
            'chr_end': int
        }
    )

    bpps = set()
    for row in rows:
        strand1 =  strand2 = row['ctg_strand']
        event_type = SVTYPES[row['type']]

        if event_type == SVTYPE.DUP:
            orient1, orient2 = (ORIENT.RIGHT, ORIENT.LEFT)
        elif event_type == SVTYPE.INS or event_type == SVTYPE.DEL:
            orient1, orient2 = (ORIENT.LEFT, ORIENT.RIGHT)
        else:
            sys.exit("Found an unexpected event type {}".format(event_type))

        if row['chr_start'] == row['chr_end']:
            row['chr_end'] += 1

        if "GL" in row['chr']:
            continue

        bpp = BreakpointPair(
            Breakpoint(row['chr'], row['chr_start'], strand=strand1, orient=orient1),
            Breakpoint(row['chr'], row['chr_end'], strand=strand2, orient=orient2),
            data={
                COLUMNS.library: args.library_id,
                COLUMNS.protocol: args.protocol,
                COLUMNS.tools: '{}_v{}'.format(__prog__,__version__),
                COLUMNS.event_type: event_type
            },
            stranded=args.stranded
        )
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

    with open(args.output, 'w') as fh:
        print('writing:', args.output)
        fh.write('## {} v{}\n'.format(__prog__, __version__))
        fh.write('## output: {0}\n'.format(args.output))
        fh.write('## overwrite: {0}\n'.format(args.overwrite))
        fh.write('## library: {0}\n'.format(args.library_id))
        fh.write('## protocol: {0}\n'.format(args.protocol))
        fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
        for row in rows:
            fh.write('\t'.join([row[c] for c in header]) + '\n')
        print("Wrote {} gene fusion events to {}".format(len(rows), args.output))
if __name__ == '__main__':
    main()
