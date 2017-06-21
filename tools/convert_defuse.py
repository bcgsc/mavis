"""
Script for converting defuse output into the MAVIS accepted input format
"""

import argparse
import os
import sys
import time
import TSV
import warnings

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, PROTOCOL
from mavis import __version__

__prog__ = os.path.basename(os.path.realpath(__file__))


def make_tsv(tsv, library_name, version, output_file_name, filter_event=True):
    """
    Function to parse the defuse tsv file and output the MAVIS tsv file
    """
    events = []
    filter_count = 0
    header, rows = TSV.read_file(tsv)

    for row in rows:
        output = {}
        o1 = ORIENT.LEFT if row['genomic_strand1'] == '+' else ORIENT.RIGHT
        o2 = ORIENT.LEFT if row['genomic_strand2'] == '+' else ORIENT.RIGHT
        output[COLUMNS.break1_orientation] = o1
        output[COLUMNS.break2_orientation] = o2
        output[COLUMNS.break1_chromosome] = row['gene_chromosome1']
        output[COLUMNS.break1_position_start] = output[COLUMNS.break1_position_end] = row['genomic_break_pos1']
        output[COLUMNS.break2_chromosome] = row['gene_chromosome2']
        output[COLUMNS.break2_position_start] = output[COLUMNS.break2_position_end] = row['genomic_break_pos2']

        if filter_event and ("GL" in output[COLUMNS.break1_chromosome] + output[COLUMNS.break2_chromosome] or
                             "MT" in output[COLUMNS.break1_chromosome] + output[COLUMNS.break2_chromosome]):
            filter_count += 1
            continue

        # use the opposite here since it is based on genes (which act similar to reads) rather than contigs
        output[COLUMNS.opposing_strands] = True if o1 == o2 else False
        output[COLUMNS.library] = library_name
        output[COLUMNS.tools] = 'deFuse_v{0}'.format(version)
        output[COLUMNS.stranded] = False
        output[COLUMNS.protocol] = PROTOCOL.TRANS
        output['defuse_spanning_read_count'] = row['span_count']
        output['defuse_split_read_count'] = row['splitr_count']
        output['defuse_cluster_id'] = row['cluster_id']
        output['defuse_probability'] = row['probability']
        bpp = BreakpointPair(
            Breakpoint(row['gene_chromosome1'], row['genomic_break_pos1'], orient=o1),
            Breakpoint(row['gene_chromosome2'], row['genomic_break_pos2'], orient=o2),
            data={},
            opposing_strands=output[COLUMNS.opposing_strands]
            )
        event_type = False
        event = {SVTYPE.DEL: row['deletion'],
                 SVTYPE.TRANS: row['interchromosomal'],
                 SVTYPE.INV: row['inversion'],
                 SVTYPE.DUP: row['eversion']}
        for key in event.keys():
            if event[key] == 'Y':
                if event_type:
                    warnings.warn("WARNING: deFuse has classified an event as more than one type")
                event_type = key
        if not event_type:
            warnings.warn("WARNING: deFuse has not given an event a classification")

        event_types = BreakpointPair.classify(bpp)
        output[COLUMNS.event_type] = event_type
        if event_type not in event_types:
            for found_event_type in event_types:
                # ignore insertions as defuse is not expected to report this type of event.
                if found_event_type == SVTYPE.INS:
                    continue
                warnings.warn("WARNING: Expected {0}, found {1}. Will add \"{2}\" for the event type".format(
                    event_types, event_type, found_event_type))
                output[COLUMNS.event_type] = found_event_type
                events.append(output)
        else:
            events.append(output)

    if filter_count > 0:
        print('{0} events have been filtered'.format(filter_count))

    elements = sort_columns(events[0].keys())
    header = "\t".join(elements)
    with open(output_file_name, 'w') as fh:
        fh.write('## {0} v{1}\n'.format(__prog__, __version__))
        fh.write('## inputs: {0}\n'.format(" ".join(sys.argv)))
        fh.write('## file generated on {0}\n'.format(time.strftime('%B %d, %Y')))
        fh.write(header+"\n")
        for event in events:
            line = []
            for element in elements:
                line.append(str(event[element]))
            fh.write("\t".join(line) + "\n")
    print("Wrote {0} gene fusion events to {1}".format(len(events), output_file_name))


def main():
    parser = argparse.ArgumentParser(
        description="Pulls gene breakpoint coordinates, strands, and orientations from deFuse result \
tsv and generates the MAVIS input file",
        add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
    required.add_argument('-o', '--output', help='path to the output file', required=True)
    required.add_argument('-l', '--library', help='The library id for the input bam file.', required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='outputs the version number')
    optional.add_argument('--tool-version', help='the version of defuse that was used in the analysis',
                          default='0.6.2')
    optional.add_argument('--no-filter', action='store_true', default=False,
                          help='turn off filtering of events that are in the "MT" and "GL" chromosomes')

    args = parser.parse_args()

    if os.path.isfile(args.input):
        make_tsv(args.input, args.library, args.tool_version, args.output, not args.no_filter)
    else:
        print("ERROR: Cannot find file: " + args.input)
        sys.exit(1)

if __name__ == "__main__":
    main()
