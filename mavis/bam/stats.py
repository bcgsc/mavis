#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import math
import statistics as stats
import warnings

import numpy as np

from .read import sequenced_strand
from ..constants import STRAND


class BamStats:

    def __init__(self, median_fragment_size, stdev_fragment_size, read_length):
        self.median_fragment_size = median_fragment_size
        self.stdev_fragment_size = stdev_fragment_size
        self.read_length = read_length
        self.stranded = False
        self.strand_determining_read = 2
        self.sdr_percent_support = None

    def add_stranded_information(self, strand_hist):
        self.stranded = True
        if strand_hist[1] > strand_hist[2]:
            self.strand_determining_read = 1
            self.sdr_percent_support = strand_hist[1] * 100 / (strand_hist[1] + strand_hist[2])
        elif strand_hist[1] < strand_hist[2]:
            self.strand_determining_read = 2
            self.sdr_percent_support = strand_hist[2] * 100 / (strand_hist[1] + strand_hist[2])
        else:
            raise AssertionError(
                'stranded values are equal. either this information was not collected or is not deterministic',
                strand_hist.items())

    def __str__(self):
        result = 'BamStats(fragment_size={0.median_fragment_size}+/-' \
                 '{0.stdev_fragment_size:.4}, read_length={0.read_length}'.format(self)
        if self.stranded:
            result += ', strand_determining_read={0.strand_determining_read}[{0.sdr_percent_support:.4}]'.format(self)
        result += ')'
        return result


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
        for val, freq in self.items():
            for i in range(0, freq):
                values.append(val)
        values.sort()
        if len(values) % 2 == 0:
            low_center = len(values) // 2
            high_center = len(values) // 2 + 1
            return (values[low_center - 1] + values[high_center - 1]) / 2
        else:
            center = len(values) // 2 + 1
            return values[center - 1]

    def distribution_stderr(self, median, fraction, error_function=lambda x, y: math.pow(x - y, 2)):
        values = []
        for val, freq in self.items():
            err = error_function(val, median)
            values.extend([err for i in range(0, freq)])
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
        for val, freq in other.items():
            result.add(val, freq)
        return result


def compute_transcriptome_bam_stats(
    bam_cache,
    annotations,
    sample_size,
    min_mapping_quality=1,
    stranded=True,
    sample_cap=10000,
    distribution_fraction=0.97
):
    """
    computes various statistical measures relating the input bam file

    Args:
        bam_file_handle (BamCache): the input bam file handle
        annotations (object): see :func:`~mavis.annotate.load_reference_genes`
        sample_size (int): the number of genes to compute stats over
        log (callable): outputs logging information
        min_mapping_quality (int): the minimum mapping quality for a read to be used
        stranded (bool): if True then reads must match the gene strand
        sample_cap (int): maximum number of reads to collect for any given sample region
        distribution_fraction (float): the proportion of the distribution to use in computing stdev

    Returns:
        BamStats: the fragment size median, stdev and the read length in a object
    """
    total_annotations = []
    for chr, anns_list in annotations.items():
        if bam_cache.valid_chr(chr):
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
    read_strand_verification = Histogram()
    read_strand_verification[1] = 0
    read_strand_verification[2] = 0

    read_lengths = []
    for gene in genes:
        for read in bam_cache.fetch(gene.chr, gene.start, gene.end, cache_if=lambda x: False, limit=sample_cap):
            if any([
                read.is_unmapped,
                read.mate_is_unmapped,
                read.mapping_quality < min_mapping_quality,
                read.next_reference_id != read.reference_id,
                read.is_secondary,
                not read.is_proper_pair
            ]):
                continue
            if stranded:
                try:
                    strand = sequenced_strand(read, 1)
                    if strand == gene.get_strand():
                        read_strand_verification.add(1)
                    else:
                        read_strand_verification.add(2)
                except ValueError:
                    pass
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
            current_frags = set()
            for spl_tx in gene.spliced_transcripts:
                try:
                    cdna_pos1 = spl_tx.convert_genomic_to_cdna(read.reference_start)
                    cdna_pos2 = spl_tx.convert_genomic_to_cdna(read.next_reference_start)
                    current_frags.add(abs(cdna_pos1 - cdna_pos2) - 2)
                except IndexError:
                    pass
            if current_frags:
                fragment_hist.add(sum(current_frags) / len(current_frags))
    read_length = stats.median(read_lengths)
    result = Histogram()
    for val, freq in fragment_hist.items():
        result.add(val + read_length, freq)
    median = result.median()
    err = result.distribution_stderr(median, distribution_fraction)
    bamstats = BamStats(median, math.sqrt(err), read_length)
    if stranded:
        bamstats.add_stranded_information(read_strand_verification)
    return bamstats


def compute_genome_bam_stats(
    bam_file_handle,
    sample_bin_size,
    sample_size,
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
    total = sum([l - sample_bin_size for l in bam_file_handle.fh.lengths])
    bins = []
    randoms = [int(n * (total - 1) + 1) for n in np.random.rand(sample_size)]
    for pos in randoms:
        template_index = 0
        for template_length in bam_file_handle.fh.lengths:
            if pos > template_length - sample_bin_size:
                pos -= (template_length - sample_bin_size)
                template_index += 1
            else:
                break
        bins.append((bam_file_handle.fh.references[template_index], pos, pos + sample_bin_size))

    hist = Histogram()
    read_lengths = []
    for bin_chr, bin_start, bin_end in bins:
        for read in bam_file_handle.fetch(bin_chr, bin_start, bin_end, limit=sample_cap, cache_if=lambda x: False):
            if any([
                read.is_unmapped,
                read.mate_is_unmapped,
                read.mapping_quality < min_mapping_quality,
                read.next_reference_id != read.reference_id,
                read.is_secondary,
                not read.is_proper_pair
            ]):
                continue
            hist[abs(read.template_length)] = hist.get(abs(read.template_length), 0) + 1
            read_lengths.append(len(read.query_sequence))
    median = hist.median()
    err = hist.distribution_stderr(median, distribution_fraction)

    return BamStats(median, math.sqrt(err), np.median(read_lengths))
