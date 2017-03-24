#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
"""
Script for converting defuse output into the MAVIS accepted input format
"""

import argparse
import TSV
import sys
import os
import time
import warnings
from mavis.constants import COLUMNS, sort_columns, ORIENT, STRAND, SVTYPE
from mavis.breakpoint import Breakpoint, BreakpointPair

__version__ = '0.0.1'
__prog__ = os.path.basename(os.path.realpath(__file__))


def make_tsv(patient_id, tsv, library_name, version=None, output_dir=""):
    """
    Function to parse the defuse tsv file and output the MAVIS tsv file
    """
    events = []
    header, rows = TSV.read_file(tsv)
    output_file_name = output_dir + patient_id + '.defuse_mavis.tsv'

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

        # use the opposite here since it is based on genes (which act similar to reads) rather than contigs
        output[COLUMNS.opposing_strands] = True if o1 == o2 else False
        output[COLUMNS.library] = library_name
        output[COLUMNS.tools] = 'deFuse_v{}'.format(version)
        output[COLUMNS.stranded] = False
        output[COLUMNS.protocol] = 'transcriptome'
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
        event_type = ''
        event = {SVTYPE.DEL: row['deletion'],
                 SVTYPE.TRANS: row['interchromosomal'],
                 SVTYPE.INV: row['inversion'],
                 SVTYPE.DUP: row['eversion']}
        for key in event.keys():
            if event[key] == 'Y':
                event_type = key
        event_types = BreakpointPair.classify(bpp)
        output[COLUMNS.event_type] = event_type
        if event_type not in event_types:
            warnings.warn("WARNING: Expected {}, found {}. Will add \"{}\" for the event type".format(
                    event_types, event_type, event_types[0]))
            # grabs the first event_type found, might not be the right choice, expects deletion
            output[COLUMNS.event_type] = event_types[0]
        events.append(output)

    elements = sort_columns(events[0].keys())
    header = "\t".join(elements)
    with open(output_file_name, 'w') as fh:
        fh.write('## {} v{}\n'.format(__prog__, __version__))
        fh.write('## inputs: {}\n'.format(" ".join(sys.argv)))
        fh.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
        fh.write(header+"\n")
        for event in events:
            line = []
            for element in elements:
                line.append(str(event[element]))
            fh.write("\t".join(line) + "\n")
    print("Wrote {} gene fusion events to {}".format(len(events), output_file_name))


def __main__():
    parser = argparse.ArgumentParser(
        description="Pulls gene breakpoint coordinates, strands, and orientations from deFuse result \
tsv and generates the MAVIS input file",
        add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--id', help='The id used in the deFuse run.', type=str, required=True)
    required.add_argument('-l', '--library_id', help='The library id for the input bam file.', required=True)
    required.add_argument('path_to_result', type=str,
                          help='absolute path to results.filtered.tsv from deFuse')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', help='the version of defuse that was used in the analysis',
                          default='0.6.2')
    optional.add_argument('-o', '--outdir', help='the directory where the output will be placed', default='')

    args = parser.parse_args()

    result = args.path_to_result
    pogid = args.id
    libraryid = args.library_id
    outputDir = args.outdir
    version = args.version

    if os.path.isfile(result):
        make_tsv(pogid, result, libraryid, version, outputDir)
    else:
        print("ERROR: Cannot find file: " + result)

if __name__ == "__main__":
    __main__()
