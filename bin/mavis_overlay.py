from mavis.illustrate.scatter import ScatterPlot
from mavis.illustrate.diagram import draw_multi_transcript_overlay
from mavis.illustrate.constants import DiagramSettings
from mavis.annotate.base import BioInterval
from mavis.annotate.file_io import load_reference_genes
from mavis import __version__
from mavis.error import DrawingFitError
import os
from mavis.interval import Interval
import argparse
import TSV
import pysam
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
        '--buffer', help='genomic length to plot on either side of the target gene', default=0, type=int
    )
    parser.add_argument(
        '-c', '--cna_file', help='path to the cna file to plot'
    )
    parser.add_argument(
        '-b', '--bam_file', help='read to add a coverage plot', default=[], action='append', nargs=3,
        metavar=('<yaxis name>', '</path/to/bam/file>', '<bin size>'))

    parser.add_argument('--markers', '-m', nargs=3, metavar=('start', 'end', 'name'), action='append', default=[])
    g = parser.add_argument_group('reference files')
    if os.environ.get('MAVIS_ANNOTATIONS', None):
        g.add_argument(
            '-a', '--annotations',
            default=os.environ.get('MAVIS_ANNOTATIONS', None),
            help='path to the reference annotations of genes, transcript, exons, domains, etc.'
        )
    else:
        g.add_argument(
            '-a', '--annotations', required=True,
            help='path to the reference annotations of genes, transcript, exons, domains, etc.'
        )
    args = parser.parse_args()
    return args


args = parse_arguments()

print('load the reference genes')
GENES = load_reference_genes(args.annotations)

genes_to_draw = []

for chr in GENES:
    for gene in GENES[chr]:
        if args.gene in gene.aliases:
            print('found gene', gene, gene.aliases)
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

d = DiagramSettings(width=1000)
d.domain_name_regex_filter = '.*'

for g in genes_to_draw:
    plots = []
    for bname, bfile, bin_size in args.bam_file:
        bin_size = int(bin_size)
        # one plot per bam
        samfile = pysam.AlignmentFile(bfile, "rb")
        points = []
        for pileupcolumn in samfile.pileup(g.chr, g.start, g.end):
            points.append((pileupcolumn.pos, pileupcolumn.n))

        temp = [x for x in range(0, len(points), bin_size)]
        temp.append(None)
        avg_points = []
        for st, end in zip(temp[0::], temp[1::]):
            pos = [x for x, y in points[st:end]]
            pos = Interval(min(pos), max(pos))
            cov = [y for x, y in points[st:end]]
            cov = Interval(sum(cov) / len(cov))
            avg_points.append((pos, cov))
        ymax = max([y.start for x, y in avg_points] + [100])
        s = ScatterPlot(avg_points, bname, ymin=0, ymax=ymax)
        print('scatter plot has', len(avg_points), 'points')
        plots.append(s)

    svg = os.path.join(args.output, '{}_{}_overlay.svg'.format(g.name, args.gene))
    markers = []

    st = max([g.start - args.buffer, 1])
    end = g.end + args.buffer
    print('window:', st, end)

    for m in args.markers:
        m = BioInterval(g.chr, int(m[0]), int(m[1]), name=m[2])
        markers.append(m)

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
    initial_width = d.width
    while d.width < initial_width + 1000:
        try:
            canvas = draw_multi_transcript_overlay(d, g, vmarkers=markers, plots=plots)
            break
        except DrawingFitError as err:
            print(repr(err), 'extending window')
            d.width += 200
    if d.width >= initial_width + 1000:
        raise DrawingFitError('could not draw', initial_width, d.width)
    d.width = initial_width
    print('writing:', svg)

    canvas.saveas(svg)
