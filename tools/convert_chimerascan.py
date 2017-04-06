#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
"""
Script for converting Chimerascan output into the MAVIS accepted input format
"""

from __future__ import print_function
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, PROTOCOL
from mavis.breakpoint import Breakpoint, BreakpointPair
import argparse
import TSV
import sys
import os
import time

__version__ = '0.0.1'
__prog__ = os.path.basename(os.path.realpath(__file__))


def parse_arguments():
    """
    Function to parse the arguments.
    """
    parser = argparse.ArgumentParser(
        description='Convert a ChimeraScan bedpe file to the MAVIS input format.',
        add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input_file', required=True,
                          help="The path to the input ChimeraScan bedpe file")
    required.add_argument('-l', '--library', required=True,
                          help="The library id of that was used as input")

    optional = parser.add_argument_group('Optional arguemts')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-o', '--output', help="The output file name", default="mavis_chimerascan.tsv")
    optional.add_argument('-v', '--version',
                          help='the version of ChimeraScan that was used in the analysis', default='0.4.5')
    optional.add_argument('--no-filter', action='store_true', default=False,
                          help='turn off filtering of events that are in the "MT" and "GL" chromosomes')
    args = parser.parse_args()
    return args


def chromosome_str(chr_repr):
    """
    Adjust the chromosome names of from the ChimeraScan output
    """
    mapping = {'23': 'X', 'M': 'MT', '24': 'Y', '25': 'MT'}
    ret_val = str(chr_repr).strip().upper().replace('CHR', '')
    if ret_val in mapping:
        ret_val = mapping[ret_val]
    return ret_val


def load_bedpe(input_bedpe, library_name, version, filter_event=True):
    """
    Function to parse the bedpe file.
    """
    events = []
    filter_count = 0
    header, rows = TSV.read_file(input_bedpe, require=['chrom5p', 'start5p', 'end5p',
                                                       'chrom3p', 'start3p', 'end3p',
                                                       'strand5p', 'strand3p'])
    for row in rows:
        output = {}
        output[COLUMNS.break1_chromosome] = chromosome_str(row['chrom5p'])
        output[COLUMNS.break2_chromosome] = chromosome_str(row['chrom3p'])

        if filter_event and ("GL" in output[COLUMNS.break1_chromosome] + output[COLUMNS.break2_chromosome] or
                             "MT" in output[COLUMNS.break1_chromosome] + output[COLUMNS.break2_chromosome]):
            filter_count += 1
            continue

    # Chimerascan's breakpoint is based on the strand of the gene, if it is + + then it is 5pend -> 3pstart
        if row['strand5p'] == '+':
            output[COLUMNS.break1_position_start] = output[COLUMNS.break1_position_end] = row['end5p']
            output[COLUMNS.break1_orientation] = ORIENT.LEFT
        else:
            output[COLUMNS.break1_position_start] = output[COLUMNS.break1_position_end] = row['start5p']
            output[COLUMNS.break1_orientation] = ORIENT.RIGHT
        if row['strand3p'] == '+':
            output[COLUMNS.break2_position_start] = output[COLUMNS.break2_position_end] = row['start3p']
            output[COLUMNS.break2_orientation] = ORIENT.LEFT
        else:
            output[COLUMNS.break2_position_start] = output[COLUMNS.break2_position_end] = row['end3p']
            output[COLUMNS.break2_orientation] = ORIENT.RIGHT

        output[COLUMNS.opposing_strands] = True if output[COLUMNS.break1_orientation] == \
            output[COLUMNS.break2_orientation] else False

        output[COLUMNS.protocol] = PROTOCOL.TRANS  # Chimerascan is assumed to only be run on transcriptomes
        output[COLUMNS.library] = library_name
        output[COLUMNS.tools] = "ChimeraScan_v"+version
        evidence = "total_spanning_frags:{}".format(row['spanning_frags'])
        output['chimerascan_evidence'] = evidence
        output[COLUMNS.stranded] = False
        bpp = BreakpointPair(
            Breakpoint(row['chrom5p'], output[COLUMNS.break1_position_start],
                       orient=output[COLUMNS.break1_orientation]),
            Breakpoint(row['chrom3p'], output[COLUMNS.break2_position_start],
                       orient=output[COLUMNS.break2_orientation]),
            data={},
            opposing_strands=output[COLUMNS.opposing_strands]
            )
        event_types = BreakpointPair.classify(bpp)
        if len(event_types) == 1 or event_types == [SVTYPE.DEL, SVTYPE.INS]:
            output[COLUMNS.event_type] = event_types[0]
        else:
            print("ERROR: event_type generated was not one of the expected event types")
            sys.exit(2)
        events.append(output)
    if filter_count > 0:
        print('{0} events have been filtered'.format(filter_count))

    return events


def write_output(events, output_file_name):
    """
    Function to write the output file.
    """
    elements = sort_columns(events[0].keys())
    header = "\t".join(elements)
    with open(output_file_name, 'w') as fh:
        fh.write('## {} v{}\n'.format(__prog__, __version__))
        fh.write('## inputs: {}\n'.format(" ".join(sys.argv)))
        fh.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
        fh.write(header + "\n")
        for event in events:
            line = []
            for element in elements:
                line.append(str(event[element]))
            fh.write("\t".join(line) + "\n")
    print("Wrote {} gene fusion events to {}".format(len(events), output_file_name))


def main():
    args = parse_arguments()

    if os.path.isfile(args.input_file):
        output = load_bedpe(args.input_file, args.library, args.version, not args.no_filter)
        write_output(output, args.output)
    else:
        print("ERROR: Cannot find file: " + args.input_file)

if __name__ == "__main__":
    main()
