"""
script for converting Trans-ABySS output file into the SVMerge accepted input format
"""
import TSV
from structural_variant.constants import ORIENT
import argparse
import os

__version__ = '0.0.1'
__prog__ = os.path.basename( os.path.realpath(__file__) )

TSV._verbose = True

def main():
    args = parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                            help='Outputs the version number')
    parser.add_argument('-f', '--overwrite', action='store_true', default=False,
                            help='set flag to overwrite existing reviewed files')
    parser.add_argument('-o', '--output', help='path to the output file', required=True)
    parser.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
    parser.add_argument('-p', '--protocol', choices=['genome', 'transcriptome'], required=True)
    parser.add_argument('-l', '--library_id', required=True)
    
    # /projects/POG/POG_data/POG098/wgs/GV2/POG098_POG098-OCT-1-unique-14-filters/POG098-OCT-1_genome_fusions_concat.tsv
    
    args = parser.parse_args()
    
    if os.path.exists(args.output) and not args.overwrite:
        print('error: output file {0} already exists. please use the --overwrite option'.format(args.out))
        parser.print_help()
        exit()
    if not os.path.exists(args.input):
        print('error: input file {0} does not exist'.format(args.input))
        parser.print_help()
        exit()
    
    
    header, rows = TSV.read_file(
            args.input,
            require = ['id'],
            split = {
                'breakpoint': '^(?P<chr1>[^:]+):(?P<pos1>\d+)\|(?P<chr2>[^:]+):(?P<pos2>\d+)$',
                'orientations': '^(?P<or1>[RL]),(?P<or2>[RL])$',
                'strands': '^(?P<strand1>[\+-]),(?P<strand2>[\+-])$'
                },
            cast = {'pos1': int, 'pos2': int},
            strict=False
            )
    
    output_header = [
            'start_chromosome', 
            'start_position',
            'end_chromosome',
            'end_position',
            'start_orientation', 
            'end_orientation', 
            'start_strand', 
            'end_strand',
            'protocol',
            'tool_version',
            'libraries'
            ]
    
    # split files into protocol and library separated files (NA here)
    
    with open(args.output, 'w') as fh:
        fh.write('## {1} v{0}\n'.format(__version__, __prog__))
        fh.write('## input: {0}\n'.format(args.input))
        fh.write('## output: {0}\n'.format(args.output))
        fh.write('## overwrite: {0}\n'.format(args.overwrite))
        fh.write('## library: {0}\n'.format(args.library_id))
        fh.write('## protocol: {0}\n'.format(args.protocol))
        fh.write('#' + '\t'.join(output_header) + '\n')
        output_rows = set()
        for row in rows:
            output_row = [
                    row['chr1'],
                    '{0}-{0}'.format(row['pos1']),
                    row['chr2'],
                    '{0}-{0}'.format(row['pos2']),
                    row['or1'],
                    row['or2'],
                    row['strand1'],
                    row['strand2'],
                    args.protocol,
                    '{1}_v{0}'.format(__version__, __prog__),
                    args.library_id
                    ]
            output_rows.add(tuple(output_row))
        for row in output_rows:
            fh.write('\t'.join([x for x in row]) + '\n')

if __name__ == '__main__':
    main()
