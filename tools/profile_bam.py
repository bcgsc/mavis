#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import argparse
import os
from structural_variant import __version__
import multiprocessing
import pysam
import math


__prog__ = os.path.basename(os.path.realpath(__file__))
MAX_POOL_SIZE = 5
MIN_MAPPING_QUALITY = 20
BAM = None


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

if __name__ == '__main__':
    def parse_arguments():
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '-v', '--version', action='version', version='%(prog)s version ' + __version__,
            help='Outputs the version number'
        )
        parser.add_argument(
            '-b', '--bamfile', required=True,
            help='path to the bam file you want to profile'
        )
        parser.add_argument(
            '--chromosomes', '-c', nargs='*', default=None,
            help='pick chromosomes to use as samples for profiling the bam')
        parser.add_argument(
            '--bin_size', type=int, default=1000,
            help='number of reads to consume per bin'
        )
        parser.add_argument(
            '--inter_bin_size', type=int, default=100000,
            help='space to skip between bins'
        )
        parser.add_argument(
            '--sample_cap', type=int, default=1000,
            help='max number of reads to parse per bin')
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
                bin_end = length - args.bin_size
                break
        if bin_start == bin_end:
            raise KeyError('did not find reference template', chr)
        while True:
            if bin_start >= bin_end:
                break
            try:
                sample_count = 0
                for read in BAM.fetch(chr, bin_start, bin_start + args.bin_size, multiple_iterators=True):
                    if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                            or read.next_reference_id != read.reference_id or read.is_secondary \
                            or not read.is_proper_pair:
                        continue
                    sample_count += 1
                    if sample_count > args.sample_cap:
                        break
                    i = abs(read.template_length)
                    hist[i] = hist.get(i, 0) + 1
                bin_start += args.bin_size + args.inter_bin_size
            except ValueError:
                break
        return chr, hist
    insert_size_by_chr = {}
    args = parse_arguments()
    BAM = pysam.AlignmentFile(args.bamfile)
    insert_sizes = None
    references = args.chromosomes if args.chromosomes is not None else [r for r in BAM.references]
    with multiprocessing.Pool(MAX_POOL_SIZE) as p:
        insert_sizes = p.map(profile_chr, references)
    hist = sum_histograms([v for k, v in insert_sizes])
    print('\nFINAL')
    a = histogram_average(hist)
    print('average                    {:.2f}'.format(a))
    s = histogram_stderr(hist, a)
    print('average stdev              {:.2f}'.format(math.sqrt(s)))
    m = histogram_median(hist)
    print('median                     {}'.format(m))
    
    e80 = histogram_distrib_stderr(hist, m, 0.8)
    print('median distrib[0.80] stdev {:.2f}'.format(math.sqrt(e80)))
    e90 = histogram_distrib_stderr(hist, m, 0.9)
    print('median distrib[0.90] stdev {:.2f}'.format(math.sqrt(e90)))
    e95 = histogram_distrib_stderr(hist, m, 0.95)
    print('median distrib[0.95] stdev {:.2f}'.format(math.sqrt(e95)))
    e99 = histogram_distrib_stderr(hist, m, 0.99)
    print('median distrib[0.99] stdev {:.2f}'.format(math.sqrt(e99)))
    e100 = histogram_stderr(hist, m)
    print('median distrib[1.00] stdev {:.2f}'.format(math.sqrt(e100)))

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
        ax.set_xlabel('abs insert size')
        ax.set_ylabel('frequency')
        ax.set_xticks([0, m, max(x)])
        ax.set_title('histogram of absolute insert sizes')
        plt.grid(True)
        x = numpy.linspace(min(hist), max(hist), 100)
        ax.plot(x, mlab.normpdf(x, a, math.sqrt(s)), color='r', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e100)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e95)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e90)), color='b', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e80)), color='b', linewidth=2)
        plt.savefig('insert_sizes_histogram_chr22.svg')
        print('wrote figure: insert_sizes_histogram_chr22.svg')
    except ImportError:
        print('cannot import matplotlib, will not generate figure')
