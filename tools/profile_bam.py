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

def histogram_distrib_stderr_correct_zero_bias(hist, median, fraction):
    arr = []
    for v, f in sorted(hist.items()):
        if v > 2 * median:
            break
        for i in range(0, f):
            arr.append(v)
    # m = subset[len(subset) // 2] if len(subset) % 2 == 0 else subset[(len(subset) - 1) // 2 + 1]
    err = 0
    total = 0
    for v in arr:
        err += math.pow(median - v, 2)
        total += 1
    return err / total


def histogram_distrib_stderr(hist, median, fraction):
    values = []
    for v, f in hist.items():
        for i in range(0, f):
            values.append(v)
    values.sort()

    count = int(len(values) * fraction)
    start = (len(values) - count) // 2

    err = 0
    total = 0
    for i in range(start, start + count):
        err += math.pow(values[i] - median, 2)
        total += 1
    return err / total


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

        args = parser.parse_args()
        return args
    
    def profile_chr(chr):
        print('profiling chr', chr)
        hist = {}
        for read in BAM.fetch(chr, multiple_iterators=True):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                    or read.next_reference_id != read.reference_id or read.is_secondary or not read.is_proper_pair:
                continue
            i = abs(read.template_length)
            hist[i] = hist.get(i, 0) + 1
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
    print('average', a)
    s = histogram_stderr(hist, a)
    print('stderr', s)
    print('stdev', math.sqrt(s))
    m = histogram_median(hist)
    print('median', m)
    e100 = histogram_stderr(hist, m)
    print('median stderr', e100)
    print('median stdev', math.sqrt(e100))
    e80 = histogram_distrib_stderr(hist, m, 0.8)
    print('median distrib[0.8] stderr', e80)
    print('median distrib[0.8] stdev', math.sqrt(e80))
    e90 = histogram_distrib_stderr(hist, m, 0.9)
    print('median distrib[0.9] stderr', e90)
    print('median distrib[0.9] stdev', math.sqrt(e90))
    e95 = histogram_distrib_stderr(hist, m, 0.95)
    print('median distrib[0.95] stderr', e95)
    print('median distrib[0.95] stdev', math.sqrt(e95))
    e = histogram_distrib_stderr_correct_zero_bias(hist, m, 0.95)
    print('median zero corrected[0.95] stderr', e)
    print('median zero corrected[0.95] stdev', math.sqrt(e))

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
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e100)), color='g', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e95)), color='g', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e90)), color='g', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e80)), color='g', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e)), color='y', linewidth=2)
        plt.savefig('insert_sizes_histogram_chr22.svg')
        print('wrote figure: insert_sizes_histogram_chr22.svg')
    except ImportError:
        print('cannot import matplotlib, will not generate figure')
