"""
Script for converting pindel output into the MAVIS accepted input format
"""

import sys
import os
import argparse
import TSV
from mavis.constants import COLUMNS, sort_columns, ORIENT, STRAND, SVTYPE, PROTOCOL

from mavis.util import get_version

__version__ = get_version()
__prog__ = os.path.basename(os.path.realpath(__file__))
default_version = '0.2.5b9'

def make_duplicate(output, events):
    """
    Function to create duplicate of an event but with opposing orientation.
    """
    duplicate = output.copy()

    if (duplicate[COLUMNS.break1_orientation] == ORIENT.RIGHT):
        duplicate[COLUMNS.break1_orientation] = ORIENT.LEFT
    else:
        duplicate[COLUMNS.break1_orientation] = ORIENT.RIGHT
    
    if (duplicate[COLUMNS.break2_orientation] == ORIENT.RIGHT):
        duplicate[COLUMNS.break2_orientation] = ORIENT.LEFT
    else:
        duplicate[COLUMNS.break2_orientation] = ORIENT.RIGHT
    
    events.append(duplicate)

def get_info_map(info_string):
    """
    Function to parse the INFO column of pindel's vcf file. 
    Returns dictionary mapping the ID and value of the INFO column contents
    """
    result = {}
    split_info = info_string.split(';')

    for info in split_info:
        temp = info.split('=')
        result[temp[0]] = temp[1]

    return result


def make_tsv(pd_files, library_name, output_file, tool_version):
    """
    Function to parse the pindel vcf file and output the MAVIS tsv file
    """

    SVTYPES = {'DEL': SVTYPE.DEL,
               'INV': SVTYPE.INV,
               'DUP': SVTYPE.DUP,
               'CNV': SVTYPE.DUP,
               'DUP:TANDEM': SVTYPE.DUP,
               'INS': SVTYPE.INS,
               'RPL': SVTYPE.INS
              }

    tsv_header = None
    with open(output_file, 'w') as fh:
        for pd_file in pd_files:
            print("Now processing: " + pd_file)
            events = []
            header, rows = TSV.read_file(pd_file)
    
            for row in rows:
                output = {}
                info = get_info_map(row["INFO"])
                output[COLUMNS.break1_chromosome] = output[COLUMNS.break2_chromosome] = row["CHROM"]
                output[COLUMNS.break1_position_start] = output[COLUMNS.break1_position_end] = row["POS"]
                output[COLUMNS.break2_position_start] = output[COLUMNS.break2_position_end] = info["END"]
                output[COLUMNS.break1_strand] = output[COLUMNS.break2_strand] = STRAND.NS
                output[COLUMNS.library] = library_name
                output[COLUMNS.protocol] = PROTOCOL.GENOME
                output[COLUMNS.tools] = 'pindel_v' + tool_version
                output[COLUMNS.stranded] = False
                output[COLUMNS.event_type] = event_type = SVTYPES[info["SVTYPE"]]
                
                if (event_type == SVTYPE.DEL or event_type == SVTYPE.INS):
                    output[COLUMNS.break1_orientation] = ORIENT.LEFT
                    output[COLUMNS.break2_orientation] = ORIENT.RIGHT
                    output[COLUMNS.opposing_strands] = False
                
                elif (event_type == SVTYPE.DUP):
                    output[COLUMNS.break1_orientation] = ORIENT.RIGHT
                    output[COLUMNS.break2_orientation] = ORIENT.LEFT
                    output[COLUMNS.opposing_strands] = False
        
                elif (event_type == SVTYPE.INV):
                    output[COLUMNS.break1_orientation] = ORIENT.RIGHT
                    output[COLUMNS.break2_orientation] = ORIENT.RIGHT
                    output[COLUMNS.opposing_strands] = True
                    make_duplicate(output, events)
                
                events.append(output)
        
            elements = sort_columns(events[0].keys())
            if (tsv_header is None):
                tsv_header = "\t".join(elements)
                fh.write(tsv_header + "\n")

            for event in events:
                line = []
                for element in elements:
                    line.append(str(event[element]))
                fh.write("\t".join(line) + "\n")

            print("Number of events: " + str(len(events)))

def main():
    parser = argparse.ArgumentParser(description='Converts pindel vcf to mavis-compatible tab file.',
            add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', required=True, help='Pindel *.vcf output', nargs='+')
    required.add_argument('-o', '--output', required=True, help='Full name of the converted output file')
    required.add_argument('-l', '--library', help='The library id for file.', required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
            help='outputs version number')
    optional.add_argument('--tool-version', help='the version of pindel that was used in the analysis',
            default=default_version)

    args = parser.parse_args()

    for f in args.input:
        if not os.path.isfile(f):
            print("ERROR: Cannot find file: " + f)
            print("Exiting.")
            exit(1)

    make_tsv(args.input, args.library, args.output, args.tool_version)

if __name__ == '__main__':
    main()
