#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import argparse
import os
from structural_variant import __version__
import multiprocessing
import pysam
import math
import numpy

__prog__ = os.path.basename(os.path.realpath(__file__))
MAX_POOL_SIZE = 5
MIN_MAPPING_QUALITY = 20
BAM = None


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
        '-o', '--output',
        help='path to the output directory', required=True
    )

    args = parser.parse_args()
    return args


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
        return (values[m] + values[n]) / 2
    else:
        m = len(values) // 2 + 1
        return values[m]


def histogram_median_err(hist, median):
    err = []
    for v, f in hist.items():
        for i in range(0, f):
            err.append(abs(v - median))
    err.sort()
    if len(err) % 2 == 0:
        m = len(err) // 2
        n = len(err) // 2 + 1
        return (err[m] + err[n]) / 2
    else:
        m = len(err) // 2 + 1
        return err[m]


def histogram_distrib_median_stderr(hist, median, fraction):
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

def histogram_stderr(hist, avg=None):
    if avg is None:
        avg = histogram_average(hist)
    running_total = 0
    N = 0
    for v, f in hist.items():
        running_total += math.pow(v - avg, 2) * f
        N += f
    return running_total / N


def histogram_abserr(hist, avg=None):
    if avg is None:
        avg = histogram_average(hist)
    running_total = 0
    N = 0
    for v, f in hist.items():
        running_total += abs(v - avg) * f
        N += f
    return running_total / N


def sum_histograms(*hist):
    sum_hist = {}
    for h in hist:
        for v, f in h.items():
            sum_hist[v] = sum_hist.get(v, 0) + f
    return sum_hist


def profile_chr(chr):
    print('profiling chr', chr)
    hist = {}
    for read in BAM.fetch(chr, multiple_iterators=True):
        if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                or read.next_reference_id != read.reference_id or read.is_secondary:
            continue
        i = abs(read.template_length)
        hist[i] = hist.get(i, 0) + 1
    a = histogram_average(hist)
    print(chr, 'average', a)
    s = histogram_stderr(hist, a)
    print(chr, 'stderr', s)
    ae = histogram_abserr(hist, a)
    print(chr, 'abserr', ae)
    print(chr, 'stdev', math.sqrt(s))
    m = histogram_median(hist)
    print(chr, 'median', m)
    print(chr, 'median error', histogram_median_err(hist, m))
    e = histogram_distrib_median_stderr(hist, m, 0.75)
    print(chr, 'median distrib[0.75] error', e)
    print(chr, 'median distrib[0.75] stdev', math.sqrt(e))
    e = histogram_distrib_median_stderr(hist, m, 0.9)
    print(chr, 'median distrib[0.9] error', e)
    print(chr, 'median distrib[0.9] stdev', math.sqrt(e))
    e = histogram_distrib_median_stderr(hist, m, 0.99)
    print(chr, 'median distrib[0.99] error', e)
    print(chr, 'median distrib[0.99] stdev', math.sqrt(e))
    return chr, hist


def main():
    global BAM
    insert_size_by_chr = {}
    args = parse_arguments()
    BAM = pysam.AlignmentFile(args.bamfile)
    insert_sizes = None
    with multiprocessing.Pool(MAX_POOL_SIZE) as p:
        insert_sizes = p.map(profile_chr, [r for r in BAM.references if not r.startswith('GL')])
    hist = sum_histograms([v for k, v in insert_sizes])
    a = histogram_average(hist)
    print('\nFINAL')
    print('total', 'average', a)
    s = histogram_stderr(hist, a)
    print('total', 'stderr', s)
    ae = histogram_abserr(hist, a)
    print('total', 'abserr', ae)
    print('total', 'stdev', math.sqrt(s))
    m = histogram_median(hist)
    print('total', 'median', histogram_median(hist))
    print('total', 'median err', histogram_median_err(hist, m))
    print('total', 'most freq', sorted(hist.items(), key = lambda x: x[1])[:2])

main()
