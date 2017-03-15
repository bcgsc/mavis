#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python

"""
Date: 2016-Aug-8
written by: cchoo@bcgsc.ca
for python 2/3
"""

from __future__ import print_function

import argparse,TSV,sys,os,time

# Script to take in a chimerascan input file and output a svmerge output file

# Expected format is
# break1_chromosome      <1-22,X,Y,MT>
# break1_position_start        <int>
# break1_position_end          <int>
# start_strand  <+,-,NS>        the reference strand aligned to #not required
# start_orientation     <L,R,NS>
# break2_chromosome      <1-22,X,Y,MT>
# break2_position_start         <int>
# break2_position_end         <int>
# end_strand    <+,-,NS>        the reference strand aligned to #not required
# end_orientation       <L,R,NS>
# protocol      <genome or transcriptome>
# library     library id's
# tools  <tool name>_<tool version number>
# tool_evidence         free-form

#ChimeraScan bedpe has the following as output
#chrom5p        start5p end5p   chrom3p start3p end3p   chimera_cluster_id      score   strand5p        strand3p        transcript_ids_5p       transcript_ids_3p       genes5p genes3p type    distance        total_frags     spanning_frags  unique_alignment_positions      isoform_fraction_5p     isoform_fraction_3p     breakpoint_spanning_reads       chimera_ids

versions = ['0.4.5']

def parse_arguments():
    # read the command-line arguments
    parser = argparse.ArgumentParser(description = 'Convert a ChimeraScan bedpe file to the SVmerge pre processed output.',
                                 add_help=False,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i','--input_file', required=True, help="The path to the input ChimeraScan bedpe file")
    required.add_argument('-l','--library', required=True, help="The library id of that was used as input to the ChimeraScan")

    optional = parser.add_argument_group('Optional arguemts')
    optional.add_argument('-h','--help', action='help', help='Show this help message and exit')
    optional.add_argument('-o','--output_file', help="The output file name", default="svmerge_chimerascan.tsv")
    optional.add_argument('-v', '--version', help='the version of defuse that was used in the analysis', default='0.4.5')
    args = parser.parse_args()
    return args

def chromosome_str(chr_repr):
    ret_val = str(chr_repr).strip()
    ret_val = chr_repr.replace('chr', '').replace('Chr', '').replace('CHR', '')
    ret_val.replace('23', 'X').replace('24', 'Y').replace('25', 'MT')
    ret_val.replace('M', 'MT')
    return ret_val


def load_bedpe(input_bedpe, library_name, version):
    """
    TODO
    """
    events=[]
    header, rows = TSV.read_file(input_bedpe, require=['chrom5p', 'start5p', 'end5p',
                                                       'chrom3p', 'start3p', 'end3p',
                                                       'strand5p', 'strand3p'])
    for row in rows:
        output = {}
        output['break1_chromosome'] = chromosome_str(row['chrom5p'])
        output['break2_chromosome'] = chromosome_str(row['chrom3p'])

        ## Chimerascan's breakpoint is based on the strand of the gene, if it is + + then it is 5pend -> 3pstart
        if row['strand5p'] == '+':
            output['break1_position_start'] = output['break1_position_end']= row['end5p']
            output['start_orientation'] = 'L'
        else:
            output['break1_position_start'] = output['break1_position_end'] =row['start5p']
            output['start_orientation'] = 'R'
        if row['strand3p'] == '+':
            output['break2_position_start'] = output['break2_position_end'] = row['start3p']
            output['end_orientation'] = 'L'
        else:
            output['break2_position_start'] = output['break2_position_end'] = row['end3p']
            output['end_orientation'] = 'R'

        output['opposing_strands'] = True if output['start_orientation'] == output['end_orientation'] else False
        output['protocol'] = "transcriptome" #Chimerascan is assumed to only be run on transcriptomes
        output['library'] = library_name
        output['tools'] = "ChimeraScan_v"+version
        evidence = "total_frags:{}".format(row['spanning_frags'])
        output['tool_evidence'] = evidence
        output['stranded'] = False
        events.append(output)
    return events

def write_output(events, output_file_name):
    """
    TODO
    """
    elements=["break1_chromosome", "break1_position_start", "break1_position_end", "start_orientation",
              "break2_chromosome", "break2_position_start", "break2_position_end", "end_orientation",
              "opposing_strands", "protocol", "library", "stranded", "tools", "tool_evidence"]
    header="\t".join(elements)
    with open(output_file_name, 'w') as fh:
        fh.write('## inputs: {}\n'.format(" ".join(sys.argv)))
        fh.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
        fh.write(header+"\n")
        for event in events:
            line=[]
            for element in elements:
                line.append(str(event[element]))
            fh.write("\t".join(line) + "\n")
    print("Wrote {} gene fusion events to {}".format(len(events),output_file_name))
def __main__():
    args = parse_arguments()
    if args.version not in versions:
        print("ERROR: The chimerascan version given is not valid")
        sys.exit()

    if os.path.isfile(args.input_file):
        output = load_bedpe(args.input_file, args.library, args.version)
        write_output(output,args.output_file)
    else:
        print("ERROR: Cannot find file: " + args.input_file)

if __name__ == "__main__":
    __main__()
