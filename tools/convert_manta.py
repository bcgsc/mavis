#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v4.3.0/envs/python3.4/bin/python
"""
Script for converting Manta output into the MAVIS accepted input format
"""

import vcf
import argparse
import sys
import os
import time
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, STRAND, PROTOCOL

__version__ = '0.0.1'
__prog__ = os.path.basename(os.path.realpath(__file__))

SVTYPES = {'DEL': SVTYPE.DEL,
           'INV': SVTYPE.INV,
           'DUP': SVTYPE.DUP,
           'BND': SVTYPE.TRA,
           'INS': SVTYPE.INS
           }


def load_vcf(vcf_filename, library, version, filter=None):
    """
    Function to parse the manta vcf file.
    """
    events = []
    with open(vcf_filename) as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)
        for record in vcf_reader:
            event = {}
            if filter:
                if 'SR' in record.samples[0].data._fields and 'PR' in record.samples[0].data._fields:
                    paired_read = record.samples[0].data.PR[1]
                    split_read = record.samples[0].data.SR[1]
                    if int(paired_read) < 2 or int(split_read) < 2:
                        continue

            record_id = record.ID
            chrom_a = record.CHROM
            position_a = record.POS
            event_type = record.INFO['SVTYPE']

            if event_type in ['DEL', 'DUP', 'DEL', 'INS', 'INV']:
                chrom_b = str(record.CHROM)
                position_b = record.INFO['END']
            elif event_type == 'BND':  # Manta specific
                if record_id[-1] == "1":
                    continue
                chrom_b = str(record.ALT[0].chr)
                position_b = record.ALT[0].pos

            if "GL" in chrom_a or "GL" in chrom_b:  # filter GL chromosomes
                continue

            if 'CIPOS' in record.INFO:
                start_a = position_a + record.INFO['CIPOS'][0]
                end_a = position_a + record.INFO['CIPOS'][1]
            else:
                start_a = position_a
                end_a = position_a

            if 'CIEND' in record.INFO:
                start_b = position_b + record.INFO['CIEND'][0]
                end_b = position_b + record.INFO['CIEND'][1]
            else:
                start_b = position_b
                end_b = position_b

            key1 = (chrom_a, start_a, end_a)
            key2 = (chrom_b, start_b, end_b)

            if key1 > key2:
                event[COLUMNS.break1_chromosome], event[COLUMNS.break1_position_start], \
                    event[COLUMNS.break1_position_end] = key2
                event[COLUMNS.break2_chromosome], event[COLUMNS.break2_position_start], \
                    event[COLUMNS.break2_position_end] = key1
            else:
                event[COLUMNS.break1_chromosome], event[COLUMNS.break1_position_start], \
                    event[COLUMNS.break1_position_end] = key1
                event[COLUMNS.break2_chromosome], event[COLUMNS.break2_position_start], \
                    event[COLUMNS.break2_position_end] = key2

            event[COLUMNS.protocol] = PROTOCOL.GENOME
            event[COLUMNS.event_type] = SVTYPES[event_type]
            event['manta_evidence'] = str(record.ID) + " " + str(record.INFO) + " " + str(record.samples)
            event[COLUMNS.stranded] = 'False'
            event[COLUMNS.library] = library
            event[COLUMNS.tools] = "Manta_v{}".format(version)

            if event_type == 'DEL' or event_type == 'INS':
                event[COLUMNS.break1_orientation], event[COLUMNS.break2_orientation] = (ORIENT.LEFT, ORIENT.RIGHT)
                event[COLUMNS.opposing_strands] = False
                events.append(event)
            elif event_type == 'DUP':
                event[COLUMNS.break1_orientation], event[COLUMNS.break2_orientation] = (ORIENT.RIGHT, ORIENT.LEFT)
                event[COLUMNS.opposing_strands] = False
                events.append(event)
            elif event_type == 'INV':
                event[COLUMNS.break1_orientation], event[COLUMNS.break2_orientation] = (ORIENT.LEFT, ORIENT.LEFT)
                event[COLUMNS.opposing_strands] = True
                events.append(event)
                event[COLUMNS.break1_orientation], event[COLUMNS.break2_orientation] = (ORIENT.RIGHT, ORIENT.RIGHT)
                event[COLUMNS.opposing_strands] = True
                events.append(event)
            elif event_type == 'BND':
                event[COLUMNS.break1_orientation], event[COLUMNS.break2_orientation] = (ORIENT.NS, ORIENT.NS)
                event[COLUMNS.opposing_strands] = STRAND.NS
                events.append(event)

    return events

parser = argparse.ArgumentParser(
    description='Convert a Manta vcf file to the MAVIS pre processed output. \
Note this has the option to filter the results for diploid.tsv based on if an event has 2 split reads \
 and 2 flanking reads',
    add_help=False)
required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', required=True, help='Manta/DeFuse vcf file to process')
required.add_argument('-o', '--output', required=True, help='output file name')
required.add_argument('-l', '--library', required=True, help='libary name of the tumor bam')
optional = parser.add_argument_group('Optional arguments')
# optional.add_argument('-b', '--bam', help = 'path to the evidence bam file for the tumour library')
optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
optional.add_argument('-f', '--filter', action='store_true',
                      help='Filter the events based on flanking and split read evidence')
optional.add_argument('-v', '--version', help='the version of Manta that was used in the analysis', default='1.0.0')
if len(sys.argv) == 1:
    parser.print_help()
    exit(2)

args = parser.parse_args()
vcf_filename = args.input
output_filename = args.output

events = load_vcf(vcf_filename, args.library, args.version, args.filter)
elements = sort_columns(events[0].keys())
header = "\t".join(elements)

with open(output_filename, 'w') as fh:
    fh.write('## {} v{}\n'.format(__prog__, __version__))
    fh.write('## inputs: {}\n'.format(" ".join(sys.argv)))
    fh.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
    fh.write('#{}\n'.format(header))
    for event in events:
        line = []
        for element in elements:
            line.append(str(event[element]))
        fh.write("{}\n".format("\t".join(line)))
    print("Wrote {} gene fusion events {}".format(len(events), output_filename))
