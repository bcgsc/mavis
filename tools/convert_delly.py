#! /projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v4.3.0/envs/python3.4/bin/python
# Requires the PyVcf (vcf) non-standard module.
##! /projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python3
"""Script to convert DELLY result to VCF for SV Merging.

APA-626
"""
import argparse
import datetime
import glob
import os
import pandas
import pdb
import pprint
import sys

import vcf

import __init__
from ProjectInfo.chromosome_naming import chromosome_standard_str
from delly_sv_report import SVTYPES

__version__ = __init__.version
__program__ = '{} v{}'.format(os.path.abspath(os.path.realpath(__file__)), __version__)

pandas.set_option('display.width', None)
pandas.set_option('display.max_colwidth', 100)

REQUIRED_COLUMNS = [
    'break1_chromosome',
    'break1_position_start',
    'break1_position_end',
    'break1_strand',
    'break1_orientation',
    'break2_chromosome',
    'break2_position_start',
    'break2_position_end',
    'break2_strand',
    'break2_orientation',
    'opposing_strands',
    'stranded',
    'library',
    'protocol',
    'tools',
    'event_type',
]


def delly_vcf_to_tsv(delly_vcf_list, output_filename=None):
    """
    Converts from the DELLY VCF format to the SV_Merging TSV format

    # Test POG181
    >>> test_vcf_files = glob.glob('/projects/POG/POG_data/POG181/sv/delly/delly-0.6.1/P00481_P00512/original_delly_vcfs/*.vcf')
    >>> table = delly_vcf_to_tsv(test_vcf_files)
    >>> table[REQUIRED_COLUMNS].head() # doctest: +NORMALIZE_WHITESPACE
         break1_chromosome       start_position break1_strand break1_orientation break2_chromosome         end_position break2_strand break2_orientation protocol library  tool_version                                                                              tool_evidence
    2415                1  100259670-100260814            +                 L              1  100260303-100261447          +               R   genome    P00512  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    2416                1  100259670-100260814            +                 L              1  100260303-100261447          +               R   genome    P00481  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    2417                1  101388177-101388873            +                 L              1  101388583-101389279          +               R   genome    P00512  DELLY_v0.6.1          {'filters': 'PASS', 'flanking_reads': 6, 'mapping_quality': 60, 'split_reads': 0}
    2418                1  102283435-102284043            +                 L              1  102283842-102284450          +               R   genome    P00512  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    2419                1  102283435-102284043            +                 L              1  102283842-102284450          +               R   genome    P00481  DELLY_v0.6.1       {'filters': 'LowQual', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}

    >>> delly073 = delly_vcf_to_tsv(['/projects/POG/POG_data/POG129/sv/delly/delly-0.7.3/P00279_P00291/Somatic_Germline_Quality_tagged.vcf',])
    >>> delly073[REQUIRED_COLUMNS].head()
        break1_chromosome       start_position break1_strand break1_orientation break2_chromosome         end_position break2_strand break2_orientation protocol library  tool_version                                                                                      tool_evidence
    704                1  100092388-100093408            +                 L              1  100093216-100094236          +               R   genome    P00279  DELLY_v0.7.3          {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    705                1  100092388-100093408            +                 L              1  100093216-100094236          +               R   genome    P00291  DELLY_v0.7.3          {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    706                1  100222492-100222598            +                 L              1  100441310-100441416          +               R   genome    P00279  DELLY_v0.7.3  {'filters': 'PASS;Quality;Somatic', 'flanking_reads': 7, 'mapping_quality': 37, 'split_reads': 0}
    707                1  100359771-100360709            +                 L              1  100360500-100361438          +               R   genome    P00279  DELLY_v0.7.3          {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    708                1  100359771-100360709            +                 L              1  100360500-100361438          +               R   genome    P00291  DELLY_v0.7.3          {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    >>> delly061 = delly_vcf_to_tsv(['/projects/POG/POG_data/POG129/sv/delly/delly-0.6.1/P00279_P00291/Somatic_Germline_Quality_tagged.vcf',])
    >>> delly061[REQUIRED_COLUMNS].head()
        break1_chromosome       start_position break1_strand break1_orientation break2_chromosome         end_position break2_strand break2_orientation protocol library  tool_version                                                                              tool_evidence
    427                1  100092287-100093509            +                 L              1  100093115-100094337          +               R   genome    P00279  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    428                1  100092287-100093509            +                 L              1  100093115-100094337          +               R   genome    P00291  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    429                1  100222253-100222879            +                 L              1  100441047-100441673          +               R   genome    P00279  DELLY_v0.6.1  {'filters': 'PASS;Somatic', 'flanking_reads': 5, 'mapping_quality': 37, 'split_reads': 2}
    430                1  100359540-100360940            +                 L              1  100360269-100361669          +               R   genome    P00291  DELLY_v0.6.1       {'filters': 'LowQual', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    431                1  100375425-100376537            +                 L              1  100376170-100377282          +               R   genome    P00279  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 3, 'mapping_quality': 60, 'split_reads': 0}
    >>> delly_table = delly_vcf_to_tsv(['/projects/POG/POG_data/POG129/sv/delly/delly-0.6.1/P00279_P00291/Somatic_Germline_Quality_tagged.vcf', '/projects/POG/POG_data/POG129/sv/delly/delly-0.7.3/P00279_P00291/Somatic_Germline_Quality_tagged.vcf'])
    >>> delly_table[REQUIRED_COLUMNS].head()
          break1_chromosome       start_position break1_strand break1_orientation break2_chromosome         end_position break2_strand break2_orientation protocol library  tool_version                                                                              tool_evidence
    427                  1  100092287-100093509            +                 L              1  100093115-100094337          +               R   genome    P00279  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    428                  1  100092287-100093509            +                 L              1  100093115-100094337          +               R   genome    P00291  DELLY_v0.6.1  {'filters': 'LowQual;PASS', 'flanking_reads': 2, 'mapping_quality': 60, 'split_reads': 0}
    19437                1  100092388-100093408            +                 L              1  100093216-100094236          +               R   genome    P00279  DELLY_v0.7.3  {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    19438                1  100092388-100093408            +                 L              1  100093216-100094236          +               R   genome    P00291  DELLY_v0.7.3  {'filters': 'LowQual;PASS', 'flanking_reads': 1, 'mapping_quality': 60, 'split_reads': 0}
    429                  1  100222253-100222879            +                 L              1  100441047-100441673          +               R   genome    P00279  DELLY_v0.6.1  {'filters': 'PASS;Somatic', 'flanking_reads': 5, 'mapping_quality': 37, 'split_reads': 2}
    """
    rows = []
    info = []
    timestamp = datetime.datetime.now().strftime('%H:%M of %Y%m%d')
    info.append("Created by: {}".format(__program__))
    info.append("Created on: {}".format(timestamp))
    input_files = [os.path.abspath(os.path.realpath(fn)) for fn in delly_vcf_list]
    info.append("Input files:")
    for fn in input_files:
        info.append("\t{}".format(fn))
    for vcf_fn in delly_vcf_list:

        vcf_reader = vcf.Reader(filename=vcf_fn)
        for record in vcf_reader:
            extra_info = {}
            extra_info['filename'] = vcf_fn
            extra_info['record_id'] = record.ID
            extra_info['mapping_qualitiy'] = record.INFO['MAPQ']
            extra_info['SVTYPE'] =  record.INFO['SVTYPE']
            extra_info['CT'] = record.INFO['CT']

            call = {}
            # Should be a semi-colon delimited list of <tool name>_<tool version>
            call['tools'] = record.INFO['SVMETHOD'].replace('EMBL.DELLY', 'DELLY_')
            call['break1_chromosome'] = chromosome_standard_str(record.CHROM)
            # Prevent minimal start position from being negative
            call['break1_position_start'] = max(1, record.POS + record.INFO['CIPOS'][0])
            call['break1_position_end'] = record.POS + record.INFO['CIPOS'][1]
            call['protocol'] = 'genome' # Just hardcoded by DELLY usage
            call['break2_chromosome'] = chromosome_standard_str(record.INFO['CHR2'])
            end_position = (record.INFO['END'] + record.INFO['CIEND'][0], record.INFO['END'] + record.INFO['CIEND'][1])
            call['break2_position_start'] = record.INFO['END'] + record.INFO['CIEND'][0]
            call['break2_position_end'] = record.INFO['END'] + record.INFO['CIEND'][1]
            call['comments'] = repr(extra_info)
            call['stranded'] = False  # Never stranded for genomes
            # Orientations on the genome are somewhat ambiguous
            # only returning half the possible orientations for now as
            # only one will be used in clustering.
            if record.INFO['CT'] == "3to5":  # Deletions
                call['break1_strand'] = call['break2_strand'] = '+'
                call['break1_orientation']   = 'L'
                call['break2_orientation']     = 'R'
                call['opposing_strands'] = False
                # or
                #call['break1_strand'] = call['break2_strand'] = '-'
                #call['break1_orientation']   = 'R'
                #call['break2_orientation']     = 'L'
            elif record.INFO['CT'] == "5to3":  # Tandem Duplication
                call['break1_strand'] = call['break2_strand'] = '+'
                call['break1_orientation']   = 'R'
                call['break2_orientation']     = 'L'
                call['opposing_strands'] = False
                # or
                #call['break1_strand'] = call['break2_strand'] = '-'
                #call['break1_orientation']   = 'L'
                #call['break2_orientation']     = 'R'
            elif record.INFO['CT'] == "3to3":  # Inversion
                call['break1_strand']        = '+'
                call['break1_orientation'] = call['break2_orientation']  = 'L'
                call['break2_strand']          = '-'
                call['opposing_strands'] = True
                # or
                #call['break1_strand']        = '-'
                #call['break1_orientation'] = call['break2_orientation']  = 'R'
                #call['break2_strand']          = '+'
            elif record.INFO['CT'] == "5to5": # Inversion
                call['break1_strand']        = '-'
                call['break1_orientation'] = call['break2_orientation']  = 'R'
                call['break2_strand']          = '+'
                call['opposing_strands'] = True
                # or
                #call['break1_strand']        = '+'
                #call['break1_orientation'] = call['break2_orientation']  = 'L'
                #call['break2_strand']          = '-'
            elif record.INFO['CT'] == "NtoN":  # Blank No-data
                # Insertions in v0.7.3 are NtoN
                call['break1_strand'] = call['break2_strand'] = '?'
                call['break1_orientation'] = call['break2_orientation']  = '?'
                call['opposing_strands'] = False
            else:
                # Can the following case and opposite not happen?
                #call['break1_strand']        = '-'
                #call['break1_orientation']   = 'R'
                #call['break2_strand']          = '+'
                #call['break2_orientation']     = 'L'
                raise ValueError("Unrecognized record['CT'] value of '{}'".format(record.INFO['CT']))
            call['event_type'] = SVTYPES[record.INFO['SVTYPE']]
            if record.INFO['SVTYPE'] == 'TRA' and call['opposing_strands']:
                call['event_type'] = 'inverted translocation'

            for sample in record.samples:
                lib = sample.sample.split('_')[0] # by naming convention
                flanking_pairs_reference= sample['DR']
                flanking_pairs_variant  = sample['DV']
                split_read_reference    = sample['RR']
                split_read_variants     = sample['RV']
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
                    new_row['library'] = lib
                    # TODO Check if start and end evidence can be different.
                    new_row['split_reads'] = split_read_variants
                    new_row['flanking_reads'] = flanking_pairs_variant
                    new_row['mapping_quality'] = extra_info['mapping_qualitiy']
                    new_row['filters'] = ';'.join(filters)
                    tool_evidence = {}
                    tool_evidence['split_reads'] = new_row['split_reads']
                    tool_evidence['flanking_reads'] = new_row['flanking_reads']
                    tool_evidence['mapping_quality'] = new_row['mapping_quality']
                    tool_evidence['filters'] = new_row['filters']
                    new_row['tool_evidence'] = pprint.pformat(tool_evidence).replace('\n', '')  # pprint is consistent for the doctest
                    rows.append(new_row)
    # Create the final table - add any missing table columns
    table = pandas.DataFrame(rows)
    header = REQUIRED_COLUMNS + [col for col in table.columns if col not in REQUIRED_COLUMNS]
    for col in header:
        if col not in table.columns:
            table[col] = '?'  # Not Specified
    # Some sorting for convenice
    table.sort_values(by=['break1_chromosome', 'break1_position_start', 'break2_chromosome', 'break2_position_start', 'flanking_reads', 'split_reads'], inplace=True)
    table = table[header].copy()
    if output_filename:
        with open(output_filename, 'w') as out_file:
            for line in info:
                out_file.write("##{}\n".format(line))
            table.to_csv(out_file, sep='\t', header=True, index=False)
        print("Wrote: '{}'".format(os.path.abspath(output_filename)))
    return table


def main():
    parser = argparse.ArgumentParser(description="DELLY VCF to TSV conversion for merging of values",
                                     fromfile_prefix_chars='@',
                                     epilog = "As an alternative to the commandline, params can be placed in a file, one per line, and specified on the commandline like '%(prog)s @filename'", 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--vcf-files', required=True, help='DELLY *.vcf output', nargs='+')
    parser.add_argument('-o', '--output-filename', help='DELLY *.tsv output', default='delly_results.tsv')
    args = parser.parse_args()
    delly_vcf_to_tsv(delly_vcf_list=args.vcf_files, output_filename=args.output_filename)


if __name__=="__main__":
    main()
