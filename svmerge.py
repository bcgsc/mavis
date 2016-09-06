import re
import TSV
import structural_variant.annotate as ann
from structural_variant.constants import *
from structural_variant.align import *
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.validate import Evidence
from structural_variant.annotate import load_reference_genome, load_reference_genes 

def load_input_file(filename):
    header, rows = TSV.read_file(filename,
            retain=[
                'start_chromosome', 
                'end_chromosome', 
                'start_orientation', 
                'end_orientation', 
                'start_strand', 
                'end_strand',
                'protocol',
                'tool_version'],
            split={
                'start_position': ('^(\d+)-(\d+)$', ['start_pos1', 'start_pos2']),
                'end_position': ('^(\d+)-(\d+)$', ['end_pos1', 'end_pos2']),
                },
            cast={'start_pos1': 'int', 'start_pos2': 'int', 'end_pos1': 'int', 'end_pos2': 'int'},
            validate={
                'start_orientation': '^{0}$'.format('|'.join([ re.escape(x) for x in ORIENT.values()])),
                'end_orientation': '^{0}$'.format('|'.join([ re.escape(x) for x in ORIENT.values()])),
                'start_strand': '^{0}$'.format('|'.join([ re.escape(x) for x in STRAND.values()])),
                'end_strand': '^{0}$'.format('|'.join([ re.escape(x) for x in STRAND.values()])),
                'tool_version': '^.+_v\d+\.\d+\.\d+$',
                'protocol': '^(genome|transcriptome)$',
                'libraries': '^[\w-]+(;[\w-]+)*$'
                },
            row_index = 'row_index'
            )
    breakpoints = []
    for row in rows:
        label = row['tool_version'] + '_' + str(row['row_index'])
        b1 = Breakpoint(
                row['start_chromosome'], 
                row['start_pos1'], 
                row['start_pos2'], 
                orient = row['start_orientation'], 
                strand = row['start_strand'],
                label = label)
        b2 = Breakpoint(
                row['end_chromosome'], 
                row['end_pos1'], 
                row['end_pos2'], 
                orient = row['end_orientation'], 
                strand = row['end_strand'],
                label = label)
        breakpoints.append(BreakpointPair(b1, b2))
    return breakpoints

print('loading the reference genes')
#ref = load_reference('/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv')
ref = load_reference_genes('/home/creisle/svn/sv_compile/trunk/reference_genes_CCDC6.tsv')
print('finished loading')
# '10:61665878' '10:43612031'
pos = 61665878
b = Breakpoint('10', pos, orient = ORIENT.LEFT, strand=STRAND.NEG)
#pos = 43612031
for gene in ref['10']:
    for t in gene.transcripts:
        if pos < t.genomic_start or pos > t.genomic_end:
            continue
        print('t', t, t.gene.aliases, t.length)
        for e in t.get_exons():
            try:
                st = t.convert_genomic_to_cdna(e.start) - t.cds_start + 1
                en = t.convert_genomic_to_cdna(e.end) - t.cds_start + 1
                print(e, st, en, st%3, en%3)
            except DiscontiuousMappingError as err:
                print(e, err)

# open up the reference sequence
#f = '/projects/seqref/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'
f = 'chr11_chr22.fa'
print('loading the human reference genome', f)
load_reference_genome(f)
print('finished loading:', f)

print('reference_genome', ann.HUMAN_REFERENCE_GENOME.keys())

bp = BreakpointPair(
            Breakpoint('11', 128664209, 128664209, orient = ORIENT.RIGHT),
            Breakpoint('22', 29684365, 29684365, orient = ORIENT.LEFT),
            stranded = False, opposing_strands = False)

bf = '/projects/analysis/analysis14/P00159/merge_bwa/125nt/hg19a/P00159_5_lanes_dupsFlagged.bam'
e = Evidence(bp, SVTYPE.TRANS, bf)
e.load_evidence()
for start, end, support, start_reads, end_reads in e.resolve_breakpoints():
    print(start, end, support, len(start_reads), len(end_reads))

for read in sorted(e.split_reads[e.break1], key=lambda x: breakpoint_pos(x) ):
    #print(read.query_name, read.cigar, read.reference_name,  read.reference_start, read.reference_end, read.reference_id)
    #print(score_cigar(read.cigar))
    print(str_cigar(read), read.query_name, breakpoint_pos(read))
print()
for read in sorted(e.split_reads[e.break2], key=lambda x: breakpoint_pos(x) ):
    #print(read.query_name, read.cigar, read.reference_name,  read.reference_start, read.reference_end, read.reference_id)
    #print(score_cigar(read.cigar))
    print(str_cigar(read), read.query_name, breakpoint_pos(read))
