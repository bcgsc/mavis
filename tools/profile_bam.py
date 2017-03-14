#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis import __version__
from mavis.constants import PROTOCOL
from mavis.annotate.file_io import load_reference_genes
from mavis.bam.stats import compute_transcriptome_bam_stats, compute_genome_bam_stats
import pysam
from configparser import ConfigParser
import mavis.pipeline.config as pconf
from mavis.pipeline.util import log
import re


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output', help='path to the output file')
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        '-l', '--library', nargs=3, metavar=('<library>', '<protocol>', '</path/to/bam>'),
        help='libraries to be profiled', required=True, action='append'
    )
    parser.add_argument('--bins', help='number of bins to create for sampling', type=int, default=100)
    parser.add_argument('--cap', help='max amount of valid reads to parse for any given sample', default=3000, type=int)
    parser.add_argument(
        '-d', '--distribution_fraction', default=0.99, type=float, 
        help='fraction of the fragment sizes distribution to use in calculating the stdev')
    g = parser.add_argument_group('genome only arguments')
    g.add_argument('--bin_size', help='size of the bins to create', default=1000, type=int)
    
    g = parser.add_argument_group('transcriptome only arguments')
    g.add_argument(
        '--annotations', help='path to the annotations file', default=os.environ.get('MAVIS_ANNOTATIONS', None))
    g.add_argument(
        '--best_transcripts_only', default=False, action='store_true',
        help='when set to true only best transcripts are used in computing the fragment size')
    
    args = parser.parse_args()
    return args


def main(args):
    # assert the inputs looks ok
    for lib, protocol, bam_file in args.library:
        PROTOCOL.enforce(protocol)
        if not os.path.exists(bam_file):
            raise UserWarning('input bam file does not exist', bam_file)
    
    # now build the configuration file
    config = ConfigParser()
    for sec in ['DEFAULTS', 'reference', 'qsub', 'illustrate']:
        config[sec] = {}
    
    for tag in pconf.REFERENCE_REQUIRED_FILES:
        config['reference'][tag] = pconf.REFERENCE_DEFAULTS[tag] if pconf.REFERENCE_DEFAULTS[tag] else ''
    for tag, val in pconf.QSUB_TAGS.items():
        config['qsub'][tag] = str(val)
    for tag, val in pconf.LIBRARY_DEFAULT_TAGS.items():
        config['DEFAULTS'][tag] = str(val)
    for tag, val in pconf.ILLUSTRATION_DEFAULTS.__dict__.items():
        config['illustrate'][tag] = str(val)
    
    for sec in config:
        for tag in config[sec]:
            if '_regex_' in tag:
                config[sec][tag] = re.sub('\$', '$$', config[sec][tag])
    
    # load the annotations file if given
    annotations = None
    if args.annotations:
        log('loading reference annotations:', args.annotations)
        annotations = load_reference_genes(args.annotations, args.best_transcripts_only)
    
    # compute the required stats for the input libraries
    for lib, protocol, bam_file in args.library:
        print('profiling:', lib, protocol, 'from', bam_file)
        bam = None
        config[lib] = {}
        config[lib]['protocol'] = protocol
        config[lib]['bam_file'] = bam_file
        try:
            bam = pysam.AlignmentFile(bam_file, 'rb')
            bamstats = None
            if protocol == PROTOCOL.TRANS:
                bamstats = compute_transcriptome_bam_stats(
                    bam, 
                    annotations=annotations, 
                    sample_size=args.bins,
                    sample_cap=args.cap, 
                    distribution_fraction=args.distribution_fraction,
                    log=log
                )
            elif protocol == PROTOCOL.GENOME:
                bamstats = compute_genome_bam_stats(
                    bam, 
                    sample_size=args.bins,
                    sample_bin_size=args.bin_size, 
                    sample_cap=args.cap, 
                    distribution_fraction=args.distribution_fraction,
                    log=log
                )
            else:
                raise ValueError('unrecognized value for protocol', protocol)
            print(
                'median', bamstats.median_fragment_size, 
                'stdev', bamstats.stdev_fragment_size, 
                'read length', bamstats.read_length)
            config[lib]['median_fragment_size'] = str(int(bamstats.median_fragment_size))
            config[lib]['stdev_fragment_size'] = str(int(bamstats.stdev_fragment_size))
            config[lib]['read_length'] = str(int(bamstats.read_length))
        finally:
            try:
                bam.close()
            except AttributeError:
                pass
    with open(args.output, 'w') as configfile:
        log('writing:', args.output)
        config.write(configfile)
    """
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
        ax.hist(x, normed=True, bins=range(min(hist), max(hist) + binwidth, binwidth), color='0.5')
        ax.set_xlabel('abs fragment size')
        ax.set_ylabel('frequency')
        ax.set_xticks([0, m, max(x)])
        ax.set_title('histogram of absolute fragment sizes')
        plt.grid(False)
        x = numpy.linspace(min(hist), max(hist), 100)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e100)), color='k', linewidth=2)
        ax.plot(x, mlab.normpdf(x, m, math.sqrt(e95)), color='b', linewidth=2)
        plt.axvline(x=m, color='k')
        stdev = math.sqrt(e95)
        for i in list(range(1, 4)) + list(range(-3, 0)):
            plt.axvline(x=m + stdev * i, color='k', linewidth=0.5)

        plt.savefig('fragment_sizes_histogram.svg')
        print('wrote figure: fragment_sizes_histogram.svg')
    except ImportError:
        print('cannot import matplotlib, will not generate figure')
    """

if __name__ == '__main__':
    main(parse_arguments())
