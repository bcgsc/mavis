#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
from ..constants import STRAND
from ..util import devnull
import math
import numpy as np
import statistics as stats
import warnings


class BamStats:
    def __init__(self, median_fragment_size, stdev_fragment_size, read_length):
        self.median_fragment_size = median_fragment_size
        self.stdev_fragment_size = stdev_fragment_size
        self.read_length = read_length


class Histogram(dict):
    def add(self, item, freq=1):
        """
        add a key to the histogram with a default frequency of 1
        """
        self[item] = self.get(item, 0) + freq

    def median(self):
        """
        flattens the histogram to compute the median value
        """
        values = []
        for v, f in self.items():
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

    def distribution_stderr(self, median, fraction, error_function=lambda x, y: math.pow(x - y, 2)):
        values = []
        for v, f in self.items():
            err = error_function(v, median)
            values.extend([err for i in range(0, f)])
        values.sort()

        end = int(len(values) * fraction)
        return sum(values[0:end]) / end

    def __add__(self, other):
        """
        sum two histograms and return the result as a new histogram

        Example:
            >>> x, y = Histogram(), Histogram()
            >>> x.add('item')
            >>> y.add('item')
            >>> x + y
            {'item': 2}
        """
        result = Histogram()
        result.update(self)
        for k, v in other.items():
            result.add(k, v)
        return result


def compute_transcriptome_bam_stats(
    bam_file_handle,
    annotations,
    sample_size,
    log=devnull,
    best_transcripts_only=False,
    min_mapping_quality=1,
    stranded=True,
    sample_cap=10000,
    distribution_fraction=0.99
):
    """
    computes various statistical measures relating the input bam file

    Args:
        bam_file_handle (pysam.AlignmentFile): the input bam file handle
        annotations (object): see :func:`~mavis.annotate.load_reference_genes`
        sample_size (int): the number of genes to compute stats over
        log (callable): outputs logging information
        best_transcripts_only (bool): limit transcripts to those flagged best only
        min_mapping_quality (int): the minimum mapping quality for a read to be used
        stranded (bool): if True then reads must match the gene strand
        sample_cap (int): maximum number of reads to collect for any given sample region
        distribution_fraction (float): the proportion of the distribution to use in computing stdev
    
    Returns:
        BamStats: the fragment size median, stdev and the read length in a object
    """
    total_annotations = []
    for chr, anns_list in annotations.items():
        if chr in bam_file_handle.references:
            total_annotations.extend(anns_list)
    
    genes = total_annotations
    if len(total_annotations) > sample_size:
        randoms = [int(n * len(total_annotations)) for n in np.random.rand(sample_size)]
        genes = [total_annotations[r] for r in randoms]
    else:
        warnings.warn(
            'insufficient annotations to match requested sample size. requested {}, but only {} annotations'.format(
                sample_size, len(total_annotations)))

    fragment_hist = Histogram()
    read_lengths = []
    for gene in genes:
        count = 0
        log('current bin (gene)', gene.name, gene.chr, gene.start, gene.end)
        for read in bam_file_handle.fetch(gene.chr, gene.start, gene.end):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < min_mapping_quality \
                    or read.next_reference_id != read.reference_id or read.is_secondary \
                    or not read.is_proper_pair:
                continue
            if stranded:
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
            read_lengths.append(len(read.query_sequence))

            if read.reference_end > read.next_reference_start:
                continue

            for t in gene.spliced_transcripts:
                try:
                    c1 = t.convert_genomic_to_cdna(read.reference_start)
                    c2 = t.convert_genomic_to_cdna(read.next_reference_start)
                    f = abs(c1 - c2) - 2
                    fragment_hist.add(abs(c1 - c2) - 2)
                except IndexError:
                    pass
            count += 1
            if count >= sample_cap:
                break
    read_length = stats.median(read_lengths)
    result = Histogram()
    for k, v in fragment_hist.items():
        result.add(k + read_length, v)
    median = result.median()
    err = result.distribution_stderr(median, distribution_fraction)
    return BamStats(median, math.sqrt(err), read_length)


def compute_genome_bam_stats(
    bam_file_handle,
    sample_bin_size,
    sample_size,
    log=devnull,
    min_mapping_quality=1,
    sample_cap=10000,
    distribution_fraction=0.99
):
    """
    computes various statistical measures relating the input bam file

    Args:
        bam_file_handle (pysam.AlignmentFile): the input bam file handle
        sample_bin_size (int): how large to make the sample bin (in bp)
        sample_size (int): the number of genes to compute stats over
        log (callable): outputs logging information
        min_mapping_quality (int): the minimum mapping quality for a read to be used
        sample_cap (int): maximum number of reads to collect for any given sample region
        distribution_fraction (float): the proportion of the distribution to use in computing stdev

    Returns:
        BamStats: the fragment size median, stdev and the read length in a object
    """
    lengths = {r: l for r, l in zip(bam_file_handle.references, bam_file_handle.lengths)}
    total = sum([l - sample_bin_size for l in bam_file_handle.lengths])
    bins = []
    randoms = [int(n * (total - 1) + 1) for n in np.random.rand(sample_size)]
    for pos in randoms:
        template_index = 0
        for c, l in zip(bam_file_handle.references, bam_file_handle.lengths):
            if pos > l - sample_bin_size:
                pos -= (l - sample_bin_size)
                template_index += 1
            else:
                break
        bins.append((bam_file_handle.references[template_index], pos, pos + sample_bin_size))

    hist = Histogram()
    read_lengths = []
    for bin_chr, bin_start, bin_end in bins:
        log('current bin', bin_chr, bin_start, bin_end)
        count = 0
        for read in bam_file_handle.fetch(bin_chr, bin_start, bin_end):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < min_mapping_quality \
                    or read.next_reference_id != read.reference_id or read.is_secondary \
                    or not read.is_proper_pair:
                continue
            elif count > sample_cap:
                break
            count += 1
            hist[abs(read.template_length)] = hist.get(abs(read.template_length), 0) + 1
            read_lengths.append(len(read.query_sequence))
    median = hist.median()
    err = hist.distribution_stderr(median, distribution_fraction)

    return BamStats(median, math.sqrt(err), np.median(read_lengths))