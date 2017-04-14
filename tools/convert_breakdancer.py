import sys
import os
import argparse
import TSV
from mavis.constants import COLUMNS, sort_columns, ORIENT, STRAND, SVTYPE, PROTOCOL
from mavis.util import get_version

__version__ = get_version()
__prog__ = os.path.basename(os.path.realpath(__file__))
default_version = '1.3.7'

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

def make_tsv(bd_list, library_name, output_file, tool_version):
    """
    Function that parses the breakdancer files and outputs a single
    mavis-compatible file. 
    """
    events = []

    for bd_file in bd_list:
        print("Now processing: " + bd_file)
        header, rows = TSV.read_file(bd_file)

        for row in rows:
            output = {}
            output[COLUMNS.break1_chromosome] = row['Chr1']
            output[COLUMNS.break2_chromosome] = row['Chr2']
            output[COLUMNS.break1_position_start] = row['Pos1']
            output[COLUMNS.break2_position_start] = row['Pos2']
            output[COLUMNS.break1_position_end] = row['Pos1']
            output[COLUMNS.break2_position_end] = row['Pos2']
            output[COLUMNS.break1_strand] = STRAND.NS
            output[COLUMNS.break2_strand] = STRAND.NS
            output[COLUMNS.library] = library_name
            output[COLUMNS.protocol] = PROTOCOL.GENOME
            output[COLUMNS.tools] = 'breakdancer_v' + tool_version
            output[COLUMNS.stranded] = False
    
            bd_event = row['Type']
    
            if (bd_event == 'DEL' or bd_event == 'INS'):
                output[COLUMNS.break1_orientation] = ORIENT.LEFT
                output[COLUMNS.break2_orientation] = ORIENT.RIGHT
                output[COLUMNS.opposing_strands] = False
    
            elif (bd_event == "ITX"): # duplication
                output[COLUMNS.break1_orientation] = ORIENT.RIGHT
                output[COLUMNS.break2_orientation] = ORIENT.LEFT
                output[COLUMNS.opposing_strands] = False
    
            elif (bd_event == "CTX"): # translocation
                output[COLUMNS.break1_orientation] = ORIENT.RIGHT
                output[COLUMNS.break2_orientation] = ORIENT.LEFT
                output[COLUMNS.opposing_strands] = False
                make_duplicate(output, events)
    
            elif (bd_event == "INV"):
                output[COLUMNS.break1_orientation] = ORIENT.RIGHT
                output[COLUMNS.break2_orientation] = ORIENT.RIGHT
                output[COLUMNS.opposing_strands] = True
                make_duplicate(output, events)
    
            else:
                output[COLUMNS.break1_orientation] = ORIENT.NS
                output[COLUMNS.break2_orientation] = ORIENT.NS
                output[COLUMNS.opposing_strands] = False # not sure if this is correct
    
            events.append(output)

    elements = sort_columns(events[0].keys())
    header = "\t".join(elements)
    with open(output_file, 'w') as fh:
        fh.write(header + "\n")
        for event in events:
            line = []
            for element in elements:
                line.append(str(event[element]))
            fh.write("\t".join(line) + "\n")

def main():
    parser = argparse.ArgumentParser(
            description='Converts breakdancer file(s) to mavis-compatible tab file.',
            add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', required=True, help='Breakdancer *.max output', nargs='+')
    required.add_argument('-o', '--output', required=True, help='Full name of the converted output file')
    required.add_argument('-l', '--library', help='The library id for the input bam file.',
            required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
            help='outputs version number')
    optional.add_argument('--tool-version', default=default_version,
            help='the version of breakdander that was used in the analysis')

    args = parser.parse_args()

    for f in args.input:
        if not os.path.isfile(f):
            print("ERROR: Cannot find file: " + f)
            print("Exiting.")
            exit(1)

    make_tsv(args.input, args.library, args.output, args.tool_version)

if __name__ == '__main__':
    main()
