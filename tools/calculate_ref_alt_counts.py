"""
Module to calculate the ref and alt counts of small indels from a mavis output file

"""
import argparse
import logging
import statistics as stats

import pysam

from mavis.annotate.file_io import load_reference_genome
from mavis.constants import SVTYPE
from mavis.util import LOG as log
from mavis.util import output_tabbed_file, read_inputs
from mavis.validate.call import EventCall


def get_read_length(bam_cache, sample_cap=20):
    """
    Determines read length from the first few reads in the bam file. 
    Assumes consistent read length in bam file.
    """
    read_lengths = []
    counter = 0
    for read in bam_cache.fetch():
        if counter == sample_cap:
            break
        read_lengths.append(len(read.query_sequence))
        counter += 1

    read_length = stats.median(read_lengths)
    return read_length


# Requires BAM files and events list
def calculate_ref_count(bpp, read_length, reference_genome, bam_cache, buffer=1):
    """
    Attempt to calculate the amount of reference reads that support the breakpoint pair
    by creating the ref/alt sequence for the breakpoint pair and substring matching against reads
    found in the breakpoint region. This sequence is extended by the amount of local repeat sequences
    to make the match more unique.

    Args:
        bpp (BreakpointPair): Breakpoint pair to calculate the ref/alt count for
        read_length (int): The read length for the reads in the bam file
        bam_cache (BamCache): the bam cache to collect the reads from
        buffer (int): The amount of sequence to add to both sides of the substring that is generated.

    Returns:
        ref_reads (set): The reads that only contained the ref sequence
        alt_reads (set): The reads that only contained the alt sequence
        ignored_reads (set): The reads that did not contain either the ref or alt sequence
        multimap_reads (set): The reads that contained both ref and alt sequences
        ref_sequence (string): The ref sequence that was generated
        alt_sequence (string): The alt sequence that was generated
    """
    alt_reads = set()
    ref_reads = set()
    ignored_reads = set()
    multimap_reads = set()

    ref = reference_genome[bpp.break1.chr]
    try:
        count, repeat = EventCall.characterize_repeat_region(bpp, reference_genome)
    except ValueError:
        # indels cannot be characterized
        count, repeat = (0, '')
    # Build the Ref and Alt sequences - This will fail if there are any SNPs in the sample that also occur near the breakpoints.
    # TODO: Support sequences from complex events
    ins_seq = bpp.untemplated_seq if bpp.untemplated_seq else ''
    if bpp.event_type == SVTYPE.DUP:
        repeat_start = bpp.break1.start - len(repeat) * count
        ref_sequence = str(ref[repeat_start - buffer: bpp.break2.end + buffer].seq)
        dup_seq = str(ref[bpp.break1.start - 1: bpp.break2.end].seq) + ins_seq
        alt_sequence = str(ref[repeat_start - buffer: bpp.break2.end].seq) + dup_seq + str(ref[bpp.break2.end:bpp.break2.end + buffer].seq)
    elif bpp.event_type == SVTYPE.INS or bpp.event_type == SVTYPE.DEL:
        repeat_start = min(bpp.break1.start - 1, bpp.break2.end - 1 - len(repeat) * (count + 1))
        ref_sequence = str(ref[repeat_start - buffer: bpp.break2.end + buffer].seq)
        alt_sequence = str(ref[repeat_start - buffer: bpp.break1.start].seq) + ins_seq + str(ref[bpp.break2.end - 1: bpp.break2.end + buffer].seq)
    else:
        raise ValueError("Cannot determine ref and alt counts for non dup/ins/del event types")
    all_reads = bam_cache.fetch(bpp.break1.chr, bpp.break1.start - read_length, bpp.break2.end + read_length)

    for read in all_reads:
        # compare the actual sequence
        ref_align = ref_sequence in read.seq
        alt_align = alt_sequence in read.seq

        if ref_align and alt_align:
            multimap_reads.add(read)

        elif ref_align:
            ref_reads.add(read)

        elif alt_align:
            alt_reads.add(read)

        else:
            ignored_reads.add(read)
    return ref_reads, alt_reads, ignored_reads, multimap_reads, ref_sequence, alt_sequence


class RefAltCalculator():
    """
    Class to calculate the ref and alt counts of a bpp for a given list of bams

    Args:
        input_bams (list of tuples): List of (library name, path to bam) tuples
        reference_genome (str or dict): Path to the reference fasta file or a dictionary created by load_reference_genome
        max_event_size (int): The maximum size of a indel event to calculate the ref/alt counts for
        buffer (int): The amount of overhang (accounting for repeats) a read must have in order to be considered
    """

    def __init__(self, input_bams, reference_genome, max_event_size=6, buffer=1):
        if isinstance(reference_genome, str):
            log('loading:', reference_genome, time_stamp=True)
            self.reference_genome = load_reference_genome(reference_genome)
        else:
            self.reference_genome = reference_genome
        self._load_bams(input_bams)
        self.bpp_cache = dict()
        self.max_event_size = max_event_size
        self.buffer = buffer

    def _load_bams(self, input_bams):
        """
        Loads the input bams with pysam and get determines their read lengths
        """
        input_bams = [(i, pysam.AlignmentFile(j, 'rb')) for i, j in input_bams]
        input_bams = [(i, get_read_length(j), j) for i, j in input_bams]
        self.input_bams = input_bams

    def calculate_ref_counts(self, bpp):
        """
        Calculates the ref and alt count for a single BreakPointPair object
        """
        if len(bpp.break1) + len(bpp.break2) > 2 or \
            max(bpp.break2.end - bpp.break1.start,
                len(bpp.untemplated_seq if bpp.untemplated_seq else '')) > self.max_event_size:
            raise ValueError("Cannot determine ref and alt count for non precise breakpoint pairs")

        if bpp not in self.bpp_cache:
            log("processing {}".format(bpp))
            data = dict()
            for name, read_length, bam in self.input_bams:
                ref, alt, ign, mul, ref_sequence, alt_sequence = calculate_ref_count(
                    bpp, read_length, self.reference_genome, bam, self.buffer)
                log(bpp, name)
                log('Calculated counts: Ref: {}, Alt: {}, Mul: {}, Ignored: {} '.format(
                    len(ref), len(alt), len(mul), len(ign)))
                log('Ref_probe: {}, Alt_probe: {}'.format(ref_sequence, alt_sequence))
                info = {'{}_ref_count'.format(name): len(ref),
                        '{}_alt_count'.format(name): len(alt),
                        '{}_ignored_count'.format(name): len(ign)}
                for key, value in info.items():
                    data[key] = value
                self.bpp_cache[bpp] = data
        for key, value in self.bpp_cache[bpp].items():
            bpp.data[key] = value
        return bpp

    def calculate_all_counts(self, input_files, output_file):
        """
        Helper method to calculate the ref and alt counts for all bpps in a file

        Args:
            input_files (list): List of mavis formatted files to use as input
            output_file (str): Path to the desired output file
        """
        processed_bpps = {}
        filtered_events = []

        bpps = read_inputs(input_files, add_default={'stranded': False})

        for bpp in bpps:
            # only use precise bpps that are within a certain event size
            try:
                processed_bpps[bpp.product_id] = self.calculate_ref_counts(bpp)
            except ValueError:
                # wrong event type to calculate a ref/alt count
                filtered_events.append(bpp)
                continue

        log('filtered {} events'.format(len(filtered_events)))

        output_tabbed_file(processed_bpps.values(), output_file)
        return processed_bpps, filtered_events


def parse_arguments():
    """
    parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Calculates the ref and alt count for small indels found in a mavis output file',
                                     add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        '-o', '--output',
        help='Path to the output file', required=True, metavar='FILEPATH')
    required.add_argument(
        '-n', '--input', required=True, metavar='FILEPATH', nargs='+',
        help='Path to the Input mavis summary file')
    required.add_argument(
        '-b', '--bam', action='append', nargs=2, default=[],
        help='Name for the library and the path to its bam file', required=True, metavar=('<name>', 'FILEPATH'))
    required.add_argument(
        '-r', '--reference', required=True, metavar='FILEPATH',
        help='Path to the Input reference genome fasta file')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('--event_size', type=int, default=6,
                          help='The maximum size of a indel event to calculate the ref/alt counts')
    optional.add_argument('--buffer', type=int, default=6,
                          help="The amount of overhang (accounting for repeats) a read must have in order to be considered")

    args = parser.parse_args()
    return args


def main():
    """
    main entry point
    """
    log_conf = {'format': '{message}', 'style': '{', 'level': 1}
    logging.basicConfig(**log_conf)
    args = parse_arguments()
    ref_alt_calculator = RefAltCalculator(args.bam, args.reference,
                                          max_event_size=args.event_size, buffer=args.buffer)
    ref_alt_calculator.calculate_all_counts(args.input, args.output)


if __name__ == "__main__":
    main()
