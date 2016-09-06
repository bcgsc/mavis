import TSV
from structural_variant.constants import ORIENT
import argparse
import os

__version__ = '0.0.1'
__prog__ = os.path.basename( os.path.realpath(__file__) )

TSV._verbose = True

#f = '/projects/seqref/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'
#print('loading the reference file', f)
#validate.load_reference(f)
args = parser = argparse.ArgumentParser()
parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                        help='Outputs the version number')
parser.add_argument('-f', '--overwrite', action='store_true', default=False,
                        help='set flag to overwrite existing reviewed files')
parser.add_argument('-o', '--output', help='path to the output file', required=True)
parser.add_argument('-n', '--input', help='path to the input file to be converted', required=True)
parser.add_argument('-p', '--protocol', choices=['genome', 'transcriptome'], required=True)

# /projects/POG/POG_data/POG098/wgs/GV2/POG098_POG098-OCT-1-unique-14-filters/POG098-OCT-1_genome_fusions_concat.tsv

args = parser.parse_args()

if os.path.exists(args.output) and not args.overwrite:
    print('error: output file {0} already exists. please use the --overwrite option'.format(args.out))
    parser.print_help()
    exit()
if not os.path.exists(args.input):
    print('error: input file {0} does not exist'.format(args.input))
    parser.print_help()
    exit()


header, rows = TSV.read_file(
        '/projects/POG/POG_data/POG098/wgs/GV2/POG098_POG098-OCT-1-unique-14-filters/POG098-OCT-1_genome_fusions_concat.tsv',
        retain = ['id'],
        split = {
            'breakpoint': ('^([^:]+):(\d+)\|([^:]+):(\d+)$', ['chr1', 'pos1', 'chr2', 'pos2']),
            'orientations': ('^([RL]),([RL])$', ['or1', 'or2']),
            'strands': ('^([\+-]),([\+-])$', ['strand1', 'strand2'])
            },
        cast = {'pos1': 'int', 'pos2': 'int'},
        strict = False
        )

output_header = [
        'start_chromosome', 
        'end_chromosome', 
        'start_orientation', 
        'end_orientation', 
        'start_strand', 
        'end_strand',
        'protocol',
        'tool_version'
        ]

with open(args.output, 'w') as fh:
    fh.write('## {1} v{0}\n'.format(__version__, __prog__))
    fh.write('## input: {0}\n'.format(args.input))
    fh.write('## output: {0}\n'.format(args.output))
    fh.write('## overwrite: {0}\n'.format(args.overwrite))
    fh.write('#' + '\t'.join(output_header) + '\n')
    for row in rows:
        output_row = [
                row['chr1'],
                row['chr2'],
                row['or1'],
                row['or2'],
                row['strand1'],
                row['strand2'],
                args.protocol,
                '{1}_v{0}'.format(__version__, __prog__)
                ]
        fh.write('\t'.join(output_row) + '\n')

exit()
"""
UNC = 10
print('loaded', len(rows), 'rows')
breakpoints = []
for row in rows:
    b1 = Breakpoint(row['chr1'], row['pos1'] - UNC, row['pos1'] + UNC, orient = row['or1'], strand = row['strand1'], label=row['id'])
    b2 = Breakpoint(row['chr2'], row['pos2'] - UNC, row['pos2'] + UNC, orient = row['or2'], strand = row['strand2'], label=row['id'])
    breakpoints.append(BreakpointPair(b1, b2))
print()
clusters = sv.cluster_breakpoints(breakpoints, r=20, k=15)

more_breakpoints = svmerge.load_input_file('/projects/trans_scratch/validations/DELLY/POG/POG098/delly-0.6.1/P00157_P00159/delly_svmerge.tsv')

print('loaded', len(breakpoints), 'breakpoints from TA and', len(more_breakpoints), 'from delly')

more_clusters = sv.cluster_breakpoints(breakpoints + more_breakpoints, r=20, k=15)

print('initially found', len(clusters), 'clusters. with delly found', len(more_clusters), 'clusters')

high_supported_clusters = sum([1 for k, c in more_clusters.items() if len(c) > 1])
print('found', high_supported_clusters, 'high_supported_clusters')
print()

bedfh = open('result.bed', 'w')

with open('result.tsv', 'w') as fh:
    fh.write('centroid_breakpoint_pair\tunion\tintersection\tsupport\tstart_dist\t_end_dist\tcumu_dist\toriginal_pairs\n')
    last_pair = None
    for pair, support in sorted(more_clusters.items(), key=lambda x: x[0].key):
        if len(support) <= 0:
            continue
        print(pair)
        d1 = -1
        d2 = -1
        if last_pair is not None \
                and last_pair.break1.chr == pair.break1.chr \
                and last_pair.break2.chr == pair.break2.chr:
            d1 = abs(last_pair.break1.pos - pair.break1.pos)
            d2 = abs(last_pair.break2.pos - pair.break2.pos)
        u1 = Interval.union([s.break1 for s in support])
        u2 = Interval.union([s.break2 for s in support])
        i1 = Interval.intersection([s.break1 for s in support])
        i2 = Interval.intersection([s.break2 for s in support])
        
        # try adding breakpoint resolution
        bf = '/projects/analysis/analysis14/P00159/merge_bwa/125nt/hg19a/P00159_5_lanes_dupsFlagged.bam'
        e = Evidence(pair, bf)
        e.load_evidence()
        for start, end, support, ss, es in e.resolve_breakpoints():
            print('putative breakpoint: {0}:{1}-{2}:{3} #{4} {5} {6}'.format(pair.break1.chr, start, pair.break2.chr, end, support, ss , es))
        e.classify()
        continue
        for bp, reads in e.resolve_breakpoint(True):
            print(bp)
            for read in reads:
                print(validate.str_cigar(read), read.query_name, read.cigar)
        print()

        for bp, reads in e.resolve_breakpoint(False):
            print(bp)
            for read in reads:
                print(validate.str_cigar(read), read.query_name, read.cigar)
        print()
        w = e._window(e.break1)
        bedfh.write('{0}\t{1}\t{2}\t{3}'.format(e.break1.chr, w[0], w[1], str(e.breakpoint_pair)))
        print()
        fh.write( '{type}\t{pair}\t{u1}==>{u2}\t{i1}==>{i2}'.format(
            type=event_type, pair=pair, u1=len(u1) , u2=len(u2), i1=len(i1 if i1 is not None else []), i2=len(i2 if i2 is not None else [])) 
                +  '\t{3}\t{0}\t{1}\t{2}\t'.format(d1, d2, d1 + d2, len(support))
                + ';'.join([str(k) for k in support]) + '\n')
        last_pair = pair
bedfh.close()
"""
