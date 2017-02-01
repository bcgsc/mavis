from structural_variant.draw import Diagram, ScatterPlot
from structural_variant.annotate.base import BioInterval
from structural_variant.annotate.file_io import load_reference_genes
from structural_variant import __version__
import os
from structural_variant.interval import Interval
import argparse
import TSV
import re


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        'gene',
        help='hugo gene name'
    )
    parser.add_argument(
        'output',
        help='path to the output folder'
    )
    parser.add_argument(
        '-b', '--buffer', help='genomic length to plot on either side of the target gene', default=0, type=int
    )
    parser.add_argument(
        '-c', '--cna_file', help='path to the cna file to plot'
    )
    
    parser.add_argument('--markers', '-m', nargs=3, metavar=('start', 'end', 'name'), action='append', default=[])
    g = parser.add_argument_group('reference files')
    g.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )
    g.add_argument(
        '-r', '--reference_genome',
        default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
        help='path to the human reference genome in fa format'
    )
    g.add_argument(
        '--template_metadata', default='/home/creisle/svn/svmerge/trunk/cytoBand.txt',
        help='file containing the cytoband template information')
    
    args = parser.parse_args()
    
    return args

args = parse_arguments()

print('load the reference genes')
GENES = load_reference_genes(args.annotations)

genes_to_draw = []

for chr in GENES:
    for gene in GENES[chr]:
        if args.gene in gene.aliases:
            print('found gene', gene)
            genes_to_draw.append(gene)

cna_by_chr = {}

if args.cna_file:
    print('loading:', args.cna_file)
    header, rows = TSV.read_file(
        args.cna_file,
        header=['chr', 'start', 'end', 'cna', 'hmm'],
        cast={
            'start': int, 
            'end': int, 
            'cna': float, 
            'chr': lambda x: re.sub('^chr', '', x)
        })

    for row in rows:
        cna_by_chr.setdefault(row['chr'], []).append(row)
        row['pos'] = int(round((row['start'] + row['end']) / 2, 0))
    
    for chr in cna_by_chr:
        cna_by_chr[chr] = sorted(cna_by_chr[chr], key=lambda x: x['start'])

d = Diagram(WIDTH=1000)

for g in genes_to_draw:
    svg = os.path.join(args.output, '{}_{}_overlay.svg'.format(g.name, args.gene))
    markers = []

    st = max([g.start - args.buffer, 1])
    end = g.end + args.buffer
    print('window:', st, end)

    for m in args.markers:
        m = BioInterval(g.chr, int(m[0]), int(m[1]), name=m[2])
        markers.append(m)

    plots = []
    if args.cna_file:
        points = []
        for row in cna_by_chr.get(g.chr, []):
            if row['start'] >= st and row['end'] <= end:
                points.append((Interval(row['start'], row['end']), Interval(row['cna'])))
        if g.chr not in cna_by_chr:
            print('chromosome not given in data', g.chr, cna_by_chr.keys())
        elif len(points) == 0:
            print('warning: did not find any cna points for this gene')
        else:
            s = ScatterPlot(points, 'cna', ymin=-1, ymax=1, hmarkers=[-1, 0, 1])
            plots.append(s)
    canvas = d.draw_ustranscripts_overlay(g, vmarkers=markers, plots=plots)
    print('writing:', svg)
    canvas.saveas(svg)

