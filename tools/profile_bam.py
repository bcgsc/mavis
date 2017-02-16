#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import argparse
import os
from structural_variant import __version__
from structural_variant.annotate.file_io import load_reference_genes
from structural_variant.error import DiscontinuousMappingError
from structural_variant.interval import Interval
from structural_variant.constants import STRAND
import multiprocessing
import pysam
import math
import numpy as np


__prog__ = os.path.basename(os.path.realpath(__file__))
POOL_SIZE = 5
MIN_MAPPING_QUALITY = 20
BAM = None
BIN_SIZE = 1000
INTER_BIN_SIZE = 100000
SAMPLE_CAP = 3000
READ_LENGTH = None
ANNOTATIONS = None
STRANDED = True
MIN_GENES_PER_POOL = 50


def histogram_average(hist):
    sum_v = 0
    sum_f = 0
    for v, f in hist.items():
        sum_v += v * f
        sum_f += f
    return sum_v / sum_f


def histogram_median(hist):
    values = []
    for v, f in hist.items():
        for i in range(0, f):
            values.append(v)
    values.sort()
    if len(values) % 2 == 0:
        m = len(values) // 2
        n = len(values) // 2 + 1
        return (values[m - 1] + values[n - 1]) / 2
    else:
        m = len(values) // 2 + 1
        return values[m - 1]


def histogram_stderr(hist, median):
    return histogram_distrib_stderr(hist, median, 1)


def histogram_distrib_stderr(hist, median, fraction, error_function=lambda x, y: math.pow(x - y, 2)):
    values = []
    for v, f in hist.items():
        for i in range(0, f):
            values.append(error_function(v, median))
    values.sort()
    
    end = int(len(values) * fraction)
    return sum(values[0:end]) / end


def histogram_abserr(hist, avg=None):
    if avg is None:
        avg = histogram_average(hist)
    running_total = 0
    N = 0
    for v, f in hist.items():
        running_total += abs(v - avg) * f
        N += f
    return running_total / N


def sum_histograms(hist):
    sum_hist = {}
    for h in hist:
        for v, f in h.items():
            sum_hist[v] = sum_hist.get(v, 0) + f
    return sum_hist


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        'bamfile',
        help='path to the bam file you want to profile'
    )
    g = parser.add_argument_group('genome')
    g.add_argument(
        '--chromosomes', '-c', nargs='*', default=None,
        help='pick chromosomes to use as samples for profiling the bam')
    g.add_argument(
        '--bin_size', type=int, default=BIN_SIZE,
        help='number of reads to consume per bin'
    )
    g.add_argument(
        '--inter_bin_size', type=int, default=INTER_BIN_SIZE,
        help='space to skip between bins'
    )
    parser.add_argument(
        '--sample_cap', type=int, default=SAMPLE_CAP,
        help='max number of reads to parse per bin')
    g = parser.add_argument_group('transcriptome')
    g.add_argument(
        '--annotations', '-a',
        default='/home/creisle/svn/svmerge/trunk/tools/ensembl69_annotations_20170207.json',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )
    g.add_argument(
        '--transcriptome', '-t', action='store_true', default=False,
        help='flag to indicate that this is a transcriptome bam (requires annotations)'
    )
    g.add_argument(
        '--genes_set_size', default=1000, type=int,
        help='randomly samples this many genes per chromosome'
    )
    g.add_argument(
        '--unstranded', default=False, action='store_true', help='this is not a stranded protocol')
    parser.add_argument(
        '-p', '--pool_size', default=2, type=int,
        help='number of processing to use'
    )
    g.add_argument(
        '-r', '--read_length', default=75, type=int,
        help='read length')
    args = parser.parse_args()
    return args


def profile_chr(chr):
    print('profiling chr', chr)
    hist = {}
    bin_start = 1
    total_reads = 0
    bin_end = bin_start
    for ref, length in zip(BAM.references, BAM.lengths):
        if ref == chr:
            bin_end = length - BIN_SIZE
            break
    if bin_start == bin_end:
        raise KeyError('did not find reference template', chr)
    while True:
        if bin_start >= bin_end:
            break
        try:
            sample_count = 0
            for read in BAM.fetch(chr, bin_start, bin_start + BIN_SIZE, multiple_iterators=True):
                if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                        or read.next_reference_id != read.reference_id or read.is_secondary \
                        or not read.is_proper_pair:
                    continue
                sample_count += 1
                if sample_count > SAMPLE_CAP:
                    break
                i = abs(read.template_length)
                hist[i] = hist.get(i, 0) + 1
            bin_start += BIN_SIZE + INTER_BIN_SIZE
        except ValueError:
            break
    return chr, hist


def profile_genes(set_size):
    genes = []
    while len(genes) < set_size:
        i = np.random.choice(range(0, len(ANNOTATIONS)))
        genes.append(ANNOTATIONS[i])
    hist = dict()
    for gene in genes:
        inserts = []
        read_count = 0
        for read in BAM.fetch(gene.chr, gene.start, gene.end, multiple_iterators=True):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                    or read.next_reference_id != read.reference_id or read.is_secondary \
                    or not read.is_proper_pair:
                continue
            if STRANDED:
                if read.is_read1:
                    if read.is_reverse:
                        if gene.get_strand() == STRAND.POS:
                            continue
                    elif gene.get_strand() == STRAND.NEG:
                        continue
                else:
                    if read.is_reverse:
                        if gene.get_strand() == STRAND.NEG:
                            continue
                    elif gene.get_strand() == STRAND.POS:
                        continue

            if read.reference_end > read.next_reference_start:
                continue

            for t in gene.spliced_transcripts:
                try:
                    c1 = t.convert_genomic_to_cdna(read.reference_start)
                    c2 = t.convert_genomic_to_cdna(read.next_reference_start)
                    inserts.append(abs(c1 - c2) - 2 + READ_LENGTH)
                except IndexError:
                    pass
            read_count += 1
            if read_count >= SAMPLE_CAP:
                break
        for ins in inserts:
            hist[ins] = hist.get(ins, 0) + 1
    return hist


def profile_transcriptome(BAM, GENES_SET_SIZE, best_transcripts_only=True):
    hist = {}
    bin_start = 1
    total_reads = 0

    print('selecting', GENES_SET_SIZE, 'random genes from', len(ANNOTATIONS), 'total genes')
    sets = []
    temp = GENES_SET_SIZE
    each = max([GENES_SET_SIZE // POOL_SIZE, MIN_GENES_PER_POOL])
    while temp >= 2 * each:
        sets.append(each)
        temp -= each
    if temp > 0:
        sets.append(temp)
    
    with multiprocessing.Pool(POOL_SIZE) as p:
        fragment_sizes = p.map(profile_genes, sets)
    return sum_histograms(fragment_sizes)


def main(args):
    global BIN_SIZE, INTER_BIN_SIZE, SAMPLE_CAP, BAM, POOL_SIZE, READ_LENGTH, ANNOTATIONS, STRANDED
    STRANDED = not args.unstranded
    POOL_SIZE = args.pool_size
    BIN_SIZE = args.bin_size
    INTER_BIN_SIZE = args.inter_bin_size
    SAMPLE_CAP = args.sample_cap
    BAM = pysam.AlignmentFile(args.bamfile)
    READ_LENGTH = args.read_length

    # run different protocols for transcriptom vs genome
    hist = dict()
    if args.transcriptome:
        print('loading:', args.annotations)
        ANNOTATIONS = load_reference_genes(args.annotations, best_transcripts_only=True)
        temp = []
        for chr, l in ANNOTATIONS.items():
            if chr not in BAM.references:
                continue
            temp.extend(l)
        ANNOTATIONS = temp
        hist = profile_transcriptome(BAM, args.genes_set_size)
    else:
        references = args.chromosomes if args.chromosomes is not None else [r for r in BAM.references]
        # map processing pools, one per chromosome
        with multiprocessing.Pool(POOL_SIZE) as p:
            fragment_sizes = p.map(profile_chr, references)
        hist = sum_histograms([v for k, v in fragment_sizes])
    
    # output the stats
    print('\nFINAL')
    print('sample size', sum(hist.values()), 'reads')
    a = histogram_average(hist)
    print('average                   \t{:.2f}'.format(a))
    s = histogram_stderr(hist, a)
    print('average stdev             \t{:.2f}'.format(math.sqrt(s)))
    m = histogram_median(hist)
    print('median                    \t{}'.format(m))
    
    e80 = histogram_distrib_stderr(hist, m, 0.8)
    print('median_distrib[0.80] stdev\t{:.2f}'.format(math.sqrt(e80)))
    e90 = histogram_distrib_stderr(hist, m, 0.9)
    print('median distrib[0.90] stdev\t{:.2f}'.format(math.sqrt(e90)))
    e95 = histogram_distrib_stderr(hist, m, 0.95)
    print('median distrib[0.95] stdev\t{:.2f}'.format(math.sqrt(e95)))
    e97 = histogram_distrib_stderr(hist, m, 0.97)
    print('median distrib[0.97] stdev\t{:.2f}'.format(math.sqrt(e97)))
    e98 = histogram_distrib_stderr(hist, m, 0.98)
    print('median distrib[0.98] stdev\t{:.2f}'.format(math.sqrt(e98)))
    e99 = histogram_distrib_stderr(hist, m, 0.99)
    print('median distrib[0.99] stdev\t{:.2f}'.format(math.sqrt(e99)))
    e100 = histogram_stderr(hist, m)
    print('median distrib[1.00] stdev\t{:.2f}'.format(math.sqrt(e100)))

    # now make a chart?
    try:
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab
        import numpy

        simple_hist = []
        for ins, freq in hist.items():
            for f in range(0, freq):
                simple_hist.append(ins)
        hist = simple_hist
        fig, ax = plt.subplots()

        x = numpy.array(sorted(hist))
        # try:
        #     from scipy.stats.mstats import normaltest
        #     k2, p = normaltest(x)
        #     print('normal test p-value', p)
        #     # take the inner fraction only

        # except ImportError:
        #     print('error: no scipy support. cannot test normal fit')
        binwidth = 10
        ax.hist(x, normed=True, bins=range(min(hist), max(hist) + binwidth, binwidth), color='b')
        ax.set_xlabel('abs fragment size')
        ax.set_ylabel('frequency')
        ax.set_xticks([0, m, max(x)])
        ax.set_title('histogram of absolute fragment sizes')
        plt.grid(True)
        x = numpy.linspace(min(hist), max(hist), 100)
        ax.plot(x, mlab.normpdf(x, a, math.sqrt(s)), color='r', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e100)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e95)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e90)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e80)), color='b', linewidth=2)
        plt.savefig('fragment_sizes_histogram_chr22.svg')
        print('wrote figure: fragment_sizes_histogram_chr22.svg')
    except ImportError:
        print('cannot import matplotlib, will not generate figure')

if __name__ == '__main__':
    main(parse_arguments())
