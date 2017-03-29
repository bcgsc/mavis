#! /projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v4.3.0/envs/python3.4/bin/python
# Requires the PyVcf (vcf) non-standard module.
"""
Script for converting DELLY VCF output into the MAVIS accepted input format
"""

import argparse
import datetime
import os
import sys
import pprint
import time
import vcf
from mavis.constants import COLUMNS, sort_columns, ORIENT, SVTYPE, STRAND, PROTOCOL

__version__ = '0.0.1'
__prog__ = os.path.abspath(os.path.realpath(__file__))


SVTYPES = {'DEL': 'deletion',
           'INV': 'inversion',
           'DUP': 'duplication',
           'TRA': 'translocation',
           'INS': 'insertion',
           }


def chromosome_standard_str(chrom):
    """Convert and of the chromosome notations to a single standard.
    """
    if not isinstance(chrom, str):
        chrom = str(chrom)
    # eliminate any 'chr' type notations or spaces
    chrom = chrom.strip().upper().replace('CHR', '')
    # Use letters X, Y instead of chromosome numbers.
    # Use 'MT' for mitochondrial
    if chrom == '23':
        chrom = 'X'
    elif chrom == '24':
        chrom = 'Y'
    elif chrom == '25':
        chrom = 'MT'
    elif chrom == 'MT':
        chrom = 'MT'
    elif chrom == 'M':
        chrom = 'MT'
    return chrom


def delly_vcf_to_tsv(delly_vcf_list, output_filename=None):
    """
    Converts from the DELLY VCF format to the SV_Merging TSV format
    """
    events = []
    for vcf_fn in delly_vcf_list:
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

            position1 = (chromosome_standard_str(record.CHROM), max(1, record.POS + record.INFO['CIPOS'][0]), record.POS + record.INFO['CIPOS'][1])
            position2 = (chromosome_standard_str(record.INFO['CHR2']),  record.INFO['END'] + record.INFO['CIEND'][0], record.INFO['END'] + record.INFO['CIEND'][1])

            if position1 > position2:
                call[COLUMNS.break1_chromosome], call[COLUMNS.break1_position_start],call[COLUMNS.break1_position_end] = position2
                call[COLUMNS.break2_chromosome], call[COLUMNS.break2_position_start],call[COLUMNS.break2_position_end] = position1

            else:
                call[COLUMNS.break1_chromosome], call[COLUMNS.break1_position_start],call[COLUMNS.break1_position_end] = position1
                call[COLUMNS.break2_chromosome], call[COLUMNS.break2_position_start],call[COLUMNS.break2_position_end] = position2

            call[COLUMNS.protocol] = PROTOCOL.GENOME  # Just hardcoded by DELLY usage
            call['delly_comments'] = repr(extra_info)
            call[COLUMNS.stranded] = False  # Never stranded for genomes

            # Orientations on the genome are somewhat ambiguous
            # only returning half the possible orientations for now as
            # only one will be used in clustering.
            if record.INFO['CT'] == "3to5":  # Deletions
                # call[COLUMNS.break1_strand] = call[COLUMNS.break2_strand] = '+'
                call[COLUMNS.break1_orientation], call[COLUMNS.break2_orientation] = ('L', 'R')
                call[COLUMNS.opposing_strands] = False
            elif record.INFO['CT'] == "5to3":  # Tandem Duplication
                # call['break1_strand'] = call['break2_strand'] = '+'
                call[COLUMNS.break1_orientation], call[COLUMNS.break2_orientation] = ('R', 'L')
                call[COLUMNS.opposing_strands] = False
            elif record.INFO['CT'] == "3to3":  # Inversion
                # call[COLUMNS.break1_strand], call[COLUMNS.break2_strand] = ('+', '-')
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = 'L'
                call[COLUMNS.opposing_strands] = True
            elif record.INFO['CT'] == "5to5":  # Inversion
                # call[COLUMNS.break1_strand], call[COLUMNS.break2_strand] = ('-', '+')
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = 'R'
                call[COLUMNS.opposing_strands] = True
            elif record.INFO['CT'] == "NtoN":  # Blank No-data
                # Insertions in v0.7.3 are NtoN
                # call['break1_strand'] = call['break2_strand'] = '?'
                call[COLUMNS.break1_orientation] = call[COLUMNS.break2_orientation] = '?'
                call[COLUMNS.opposing_strands] = False
            else:
                raise ValueError("Unrecognized record['CT'] value of '{}'".format(record.INFO['CT']))
            call[COLUMNS.event_type] = SVTYPES[record.INFO['SVTYPE']]
            if record.INFO['SVTYPE'] == 'TRA' and call['opposing_strands']:
                call[COLUMNS.event_type] = 'inverted translocation'

            for sample in record.samples:
                lib = sample.sample.split('_')[0]  # by naming convention
                flanking_pairs_reference = sample['DR']
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
                    # tool_evidence = {}
                    # tool_evidence['split_reads'] = new_row['delly_split_reads']
                    # tool_evidence['flanking_reads'] = new_row['delly_flanking_reads']
                    # tool_evidence['mapping_quality'] = new_row['delly_mapping_quality']
                    # tool_evidence['filters'] = new_row['delly_filters']
                    # new_row['delly_evidence'] = pprint.pformat(tool_evidence).replace('\n', '')
                    events.append(new_row)

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
        print("Wrote {} gene fusion events to {}".format(len(events), output_filename))


def main():
    parser = argparse.ArgumentParser(description="DELLY VCF to TSV conversion for merging of values",
                                     fromfile_prefix_chars='@',
                                     epilog="As an alternative to the commandline, params can be placed \
in a file, one per line, and specified on the commandline like '%(prog)s @filename'",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--vcf-files', required=True, help='DELLY *.vcf output', nargs='+')
    parser.add_argument('-o', '--output-filename', help='DELLY *.tsv output', default='mavis_delly.tsv')
    args = parser.parse_args()
    delly_vcf_to_tsv(delly_vcf_list=args.vcf_files, output_filename=args.output_filename)


if __name__ == "__main__":
    main()
