# Requires the PyVcf (vcf) non-standard module.
"""
Script for converting DELLY VCF output into the MAVIS accepted input format
"""

import argparse
import os
import sys
import time
import vcf

from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, PROTOCOL
from mavis.util import get_version

__version__ = get_version()
__prog__ = os.path.basename(os.path.realpath(__file__))

SVTYPES = {'DEL': SVTYPE.DEL,
           'INV': SVTYPE.INV,
           'DUP': SVTYPE.DUP,
           'TRA': SVTYPE.TRANS,
           'INS': SVTYPE.INS
           }


def chromosome_str(chrom):
    """Convert and of the chromosome notations to a single standard.
    """
    mapping = {'23': 'X', 'M': 'MT', '24': 'Y', '25': 'MT'}
    if not isinstance(chrom, str):
        chrom = str(chrom)
    # eliminate any 'chr' type notations or spaces
    chrom = chrom.strip().upper().replace('CHR', '')
    # Use letters X, Y instead of chromosome numbers.
    # Use 'MT' for mitochondrial
    if chrom in mapping:
        chrom = mapping[chrom]
    return chrom


def delly_vcf_to_tsv(delly_vcf_list, output_filename=None, filter_event=True):
    """
    Converts from the DELLY VCF format to the SV_Merging TSV format
    """
    events = []
    filter_count = 0
    for vcf_fn in delly_vcf_list:
        print('reading:', vcf_fn)
        vcf_reader = vcf.Reader(filename=vcf_fn)
        for record in vcf_reader:
            extra_info = {}
            extra_info['filename'] = vcf_fn
            extra_info['record_id'] = record.ID
            extra_info['mapping_qualitiy'] = record.INFO['MAPQ']
            extra_info['SVTYPE'] = record.INFO['SVTYPE']
            extra_info['CT'] = record.INFO['CT']

            call = {}
            # Should be a semi-colon delimited list of <tool name>_<tool version>
            call[COLUMNS.tools] = record.INFO['SVMETHOD'].replace('EMBL.DELLY', 'DELLY_')

            if filter_event and ("GL" in chromosome_str(record.CHROM) + chromosome_str(record.INFO['CHR2']) or
                                 "MT" in chromosome_str(record.CHROM) + chromosome_str(record.INFO['CHR2'])):
                filter_count += 1
                continue

            position1 = (chromosome_str(record.CHROM), max(1, record.POS + record.INFO['CIPOS'][0]),
                         record.POS + record.INFO['CIPOS'][1])
            position2 = (chromosome_str(record.INFO['CHR2']),  record.INFO['END'] + record.INFO['CIEND'][0],
                         record.INFO['END'] + record.INFO['CIEND'][1])

            # sort the positions so the orientations match
            if position1 > position2:
                call[COLUMNS.break1_chromosome], call[COLUMNS.break1_position_start], call[COLUMNS.break1_position_end] = position2
                call[COLUMNS.break2_chromosome], call[COLUMNS.break2_position_start], call[COLUMNS.break2_position_end] = position1

            else:
                call[COLUMNS.break1_chromosome], call[COLUMNS.break1_position_start], call[COLUMNS.break1_position_end] = position1
                call[COLUMNS.break2_chromosome], call[COLUMNS.break2_position_start], call[COLUMNS.break2_position_end] = position2

            call[COLUMNS.protocol] = PROTOCOL.GENOME  # Just hardcoded by DELLY usage
            call['delly_comments'] = repr(extra_info)
            call[COLUMNS.stranded] = False  # Never stranded for genomes

            # Orientations on the genome are somewhat ambiguous
            # only returning half the possible orientations for now as
            # only one will be used in clustering.
            if record.INFO['CT'] == "3to5":  # Deletions
                # call[COLUMNS.break1_strand] = call[COLUMNS.break2_strand] = '+'
                call[COLUMNS.break1_orientation], call[COLUMNS.break2_orientation] = (ORIENT.LEFT, ORIENT.RIGHT)
                call[COLUMNS.opposing_strands] = False
            elif record.INFO['CT'] == "5to3":  # Tandem Duplication
                # call['break1_strand'] = call['break2_strand'] = '+'
                call[COLUMNS.break1_orientation], call[COLUMNS.break2_orientation] = (ORIENT.RIGHT, ORIENT.LEFT)
                call[COLUMNS.opposing_strands] = False
            elif record.INFO['CT'] == "3to3":  # Inversion
                # call[COLUMNS.break1_strand], call[COLUMNS.break2_strand] = ('+', '-')
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = ORIENT.LEFT
                call[COLUMNS.opposing_strands] = True
            elif record.INFO['CT'] == "5to5":  # Inversion
                # call[COLUMNS.break1_strand], call[COLUMNS.break2_strand] = ('-', '+')
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = ORIENT.RIGHT
                call[COLUMNS.opposing_strands] = True
            elif record.INFO['CT'] == "NtoN":  # Blank No-data
                # Insertions in v0.7.3 are NtoN
                # call['break1_strand'] = call['break2_strand'] = '?'
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = ORIENT.NS
                call[COLUMNS.opposing_strands] = False
            else:
                raise ValueError("Unrecognized record['CT'] value of '{0}'".format(record.INFO['CT']))
            call[COLUMNS.event_type] = SVTYPES[record.INFO['SVTYPE']]
            if record.INFO['SVTYPE'] == 'TRA' and call['opposing_strands']:
                call[COLUMNS.event_type] = 'inverted translocation'

            for sample in record.samples:
                lib = sample.sample.split('_')[0]  # by naming convention
                flanking_pairs_reference = sample['DR']  # TODO need to determine the normal library, right now using the tumor evidence for both?
                flanking_pairs_variant = sample['DV']
                split_read_reference = sample['RR']
                split_read_variants = sample['RV']
                filters = []
                if record.FILTER:
                    if isinstance(record.FILTER, str):
                        filters.append(record.FILTER)
                    elif isinstance(record.FILTER, list):
                        filters.extend(record.FILTER)
                filters.extend(sample['FT'].split(';'))
                filters = sorted(list(set(filters)))
                # Must be some evidence for the sample
                if flanking_pairs_variant or split_read_variants:
                    new_row = {}
                    new_row.update(call)
                    new_row[COLUMNS.library] = lib
                    # TODO Check if start and end evidence can be different.
                    new_row['delly_split_reads'] = split_read_variants
                    new_row['delly_flanking_reads'] = flanking_pairs_variant
                    new_row['delly_mapping_quality'] = extra_info['mapping_qualitiy']
                    new_row['delly_filters'] = ';'.join(filters)
                    events.append(new_row)
    if filter_count > 0:
        print('{0} events have been filtered'.format(filter_count))

    elements = sort_columns(events[0].keys())
    header = "\t".join(elements)
    with open(output_filename, 'w') as fh:
        fh.write('## {0} v{1}\n'.format(__prog__, __version__))
        fh.write('## inputs: {0}\n'.format(" ".join(sys.argv)))
        fh.write('## file generated on {0}\n'.format(time.strftime('%B %d, %Y')))
        fh.write('#{0}\n'.format(header))
        for event in events:
            line = []
            for element in elements:
                line.append(str(event[element]))
            fh.write("{0}\n".format("\t".join(line)))
        print("Wrote {0} gene fusion events to {1}".format(len(events), output_filename))


def main():
    parser = argparse.ArgumentParser(
        description="Convert a DELLY VCF output file to the MAVIS input format",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #library and tool version are determined from the vcf file
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', required=True, help='the path to all the DELLY *.vcf output', nargs='+')
    required.add_argument('-o', '--output', help='path to the output file', required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                          help='outputs the version number')
    optional.add_argument('--no-filter', action='store_true', default=False,
                          help='turn off filtering of events that are in the "MT" and "GL" chromosomes')
    args = parser.parse_args()
    delly_vcf_to_tsv(delly_vcf_list=args.input, output_filename=args.output,
                     filter_event=not args.no_filter)


if __name__ == "__main__":
    main()
