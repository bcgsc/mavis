#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
from ..constants import STRAND
import math
import numpy as np
import itertools
import statistics as stats


class BamStats:
    def __init__(self, median_fragment_size, stdev_fragment_size, read_length):
        self.median_fragment_size = median_fragment_size
        self.stdev_fragment_size = stdev_fragment_size
        self.read_length = read_length


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


def histogram_distrib_stderr(hist, median, fraction, error_function=lambda x, y: math.pow(x - y, 2)):
    values = []
    for v, f in hist.items():
        for i in range(0, f):
            values.append(error_function(v, median))
    values.sort()

    end = int(len(values) * fraction)
    return sum(values[0:end]) / end


def sum_histograms(hist):
    sum_hist = {}
    for h in hist:
        for v, f in h.items():
            sum_hist[v] = sum_hist.get(v, 0) + f
    return sum_hist


def compute_transcriptome_bam_stats(
    bam_file_handle, 
    annotations, 
    sample_size, 
    log=lambda *pos, **kwargs: None,
    best_transcripts_only=False, 
    MIN_MAPPING_QUALITY=1, 
    stranded=True,
    sample_cap=10000, 
    distribution_fraction=0.99, 
    **kwargs
):
    total_annotations = []
    for chr, anns_list in annotations.items():
        if chr in bam_file_handle.references:
            total_annotations.extend(anns_list)
    
    randoms = [int(n * len(total_annotations)) for n in np.random.rand(sample_size)]
    genes = [total_annotations[r] for r in randoms]
    
    fragment_hist = dict()
    read_lengths = []
    for gene in genes:
        count = 0
        log('current bin (gene)', gene.name, gene.chr, gene.start, gene.end)
        for read in bam_file_handle.fetch(gene.chr, gene.start, gene.end):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
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
                    fragment_hist[f] = fragment_hist.get(f, 0) + 1
                except IndexError:
                    pass
            count += 1
            if count >= sample_cap:
                break
    read_length = stats.median(read_lengths)
    result = {}
    for k, v in fragment_hist.items():
        result[k + read_length] = v
    median = histogram_median(result)
    err = histogram_distrib_stderr(result, median, distribution_fraction)
    return BamStats(median, math.sqrt(err), read_length)


def compute_genome_bam_stats(
    bam_file_handle, 
    sample_bin_size, 
    sample_size,
    log=lambda *pos, **kwargs: None, 
    MIN_MAPPING_QUALITY=1, 
    sample_cap=10000, 
    distribution_fraction=0.99, 
    **kwargs
):
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
    
    hist = {}
    read_lengths = []
    for bin_chr, bin_start, bin_end in bins:
        log('current bin', bin_chr, bin_start, bin_end)
        count = 0
        for read in bam_file_handle.fetch(bin_chr, bin_start, bin_end):
            if read.is_unmapped or read.mate_is_unmapped or read.mapping_quality < MIN_MAPPING_QUALITY \
                    or read.next_reference_id != read.reference_id or read.is_secondary \
                    or not read.is_proper_pair:
                continue
            elif count > sample_cap:
                break
            count += 1
            hist[abs(read.template_length)] = hist.get(abs(read.template_length), 0) + 1
            read_lengths.append(len(read.query_sequence))
    median = histogram_median(hist)
    err = histogram_distrib_stderr(hist, median, distribution_fraction)
    
    return BamStats(median, math.sqrt(err), np.median(read_lengths))
