#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v4.3.0/envs/python3.4/bin/python
"""
#Basic Process

#Load the Manta vcf file, parse out the information provided


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
# opposing_strands
# protocol      <genome or transcriptome>
# library     library id's
# tools  <tool name>_<tool version number>
# tool_evidence         free-form

"""

import vcf
import argparse
import pysam
import sys
import time

# def getOrientationFromBam(bam_file,event,chrom,pos):
#     """
#     Takes an evidence bam file and searches for the reads that support the event. Then determines the strand and orientation of the evidence reads to determine strand and orientation for the event.

#     """
#     evidence_reads = []
#     strand = orientation = None
#     #print(event,chrom,pos)
#     for read in bam_file.fetch(chrom,pos-100,pos+100):
#         if event + "|" in dict(read.get_tags())['ZM']:
#             evidence_reads.append(read)
#     types = {}
#     for read in evidence_reads:
#         if read.is_reverse:
#             strand = "-"
#         else:
#             strand = "+"

#         if read.reference_start > pos:
#             orientation = "R"
#         elif read.reference_start <= pos and read.reference_start + read.query_length >= pos: #split reads
#             orientation = "NA"
#         elif read.reference_start < pos:
#             orientation = "L"
#         else:
#             orientation = "ok"

#         if orientation != "NA" and 'SR' not in dict(read.get_tags())['ZM']:
#             if (strand,orientation) not in types:
#                 types[(strand,orientation)] = True
#     return list(types.keys())

def load_vcf(vcf_filename, library, version, filter=None):
    """
    """
    events=[]
    with open(vcf_filename) as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)
        for record in vcf_reader:
            event = {}
            if filter:
                if 'SR'in  record.samples[0].data._fields and 'PR' in record.samples[0].data._fields:
                    paired_read = record.samples[0].data.PR[1]
                    split_read = record.samples[0].data.SR[1]
                    if int(paired_read) < 2 or  int(split_read) < 2:
                        continue

            record_id = record.ID
            chrom_a = record.CHROM
            position_a = record.POS
            event_type = record.INFO['SVTYPE']

            if event_type in ['DEL','DUP','DEL','INS','INV']:
                chrom_b = str(record.CHROM)
                position_b = record.INFO['END']
            elif event_type  == 'BND': #Manta specific
                if record_id[-1] == "1":
                    continue
                chrom_b = str(record.ALT[0].chr)
                position_b = record.ALT[0].pos

            if "GL" in chrom_a or "GL" in chrom_b: #filter GL chromosomes
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

            event['break1_chromosome'] = chrom_a
            event['break2_chromosome'] = chrom_b
            event['break1_position_start'] = start_a
            event['break1_position_end'] = end_a
            event['break2_position_start'] = start_b
            event['break2_position_end'] = end_b
            event['protocol'] = 'genome'
            event['tool_evidence'] = str(record.ID)+str(record.INFO)+str(record.samples)
            event['stranded'] = 'False'
            event['library'] = library
            event['tools'] = "Manta_v{}".format(version)

            if event_type == 'DEL' or event_type == 'INS':
                event['start_orientation'],event['end_orientation'] = ('L','R')
                event['opposing_strands'] = False
                events.append(event)
            elif event_type == 'DUP':
                event['start_orientation'],event['end_orientation'] = ('R','L')
                event['opposing_strands'] = False
                events.append(event)
            elif event_type == 'INV':
                event['start_orientation'],event['end_orientation'] = ('L','L')
                event['opposing_strands'] = False
                events.append(event)
                event['start_orientation'],event['end_orientation'] = ('R','R')
                event['opposing_strands'] = False
                events.append(event)
            elif event_type == 'BND':
                event['start_orientation'],event['end_orientation'] = ('?','?')
                event['opposing_strands'] = 'null'
                events.append(event)

    return events

parser = argparse.ArgumentParser(description = 'Convert a Manta vcf file to the SVmerge pre processed output. Note this has the option to filter the results for diploid.tsv based on if an event has 2 split reads and 2 flanking reads',
                                 add_help=False)
required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', required=True, help='Manta/DeFuse vcf file to process')
required.add_argument('-o', '--output', required=True, help='output file name')
required.add_argument('-l','--library',required=True, help='libary name of the tumor bam')
optional = parser.add_argument_group('Optional arguments')
#optional.add_argument('-b', '--bam', help = 'path to the evidence bam file for the tumour library')
optional.add_argument('-h','--help', action='help', help='Show this help message and exit')
optional.add_argument('-f', '--filter', action='store_true', help='Filter the events based on flanking and split read evidence')
optional.add_argument('-v', '--version',help='the version of Manta that was used in the analysis', default='1.0.0')
if len(sys.argv) == 1:
    parser.print_help()
    exit(2)

args = parser.parse_args()
vcf_filename = args.input
output_filename = args.output

events = load_vcf(vcf_filename, args.library, args.version, args.filter)
#if args.bam:
#    bam_file = pysam.AlignmentFile(arguments.bam,'rb')

with open(output_filename, 'w') as output_file:
    output_file.write('## inputs: {}\n'.format(" ".join(sys.argv)))
    output_file.write('## file generated on {}\n'.format(time.strftime('%B %d, %Y')))
    header=['break1_chromosome','break1_position_start','break1_position_end','start_orientation',
            'break2_chromosome','break2_position_start','break2_position_end','end_orientation',
            'opposing_strands','protocol','library','stranded','tools', 'tool_evidence']
    output_file.write('#{}\n'.format('\t'.join(header)))

    for event in events:
        line=[]
        for element in header:
            line.append(str(event[element]))
        output_file.write("{}\n".format("\t".join(line)))
    print("Wrote {} gene fusion events {}".format(len(events), output_filename))

