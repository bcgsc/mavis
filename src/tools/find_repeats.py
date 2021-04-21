"""
Script used in finding potential masking regions within a genome
"""
import argparse
import os

from mavis.annotate.base import BioInterval
from mavis.annotate.file_io import load_reference_genome
from mavis.util import log


def parse_arguments():
    """
    parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output',
        help='path to the output file', required=True, metavar='FILEPATH'
    )
    parser.add_argument(
        '-n', '--input', required=True, metavar='FILEPATH',
        help='Path to the Input reference genome fasta file'
    )
    parser.add_argument(
        '--min_length', default=20, type=int, help='Minimum total length of the repeat region to find', metavar='INT')
    parser.add_argument(
        '--repeat_seq', default='N', type=str, help='Repeat sequence to look for. Case insensitive', nargs='+')
    args = parser.parse_args()
    if args.min_length < 2:
        parser.error('argument --min_length: cannot specify a shorter repeat than 2 bases')
    if not os.path.exists(args.input):
        parser.error('argument --input: File does not exist')
    return args


def main():
    args = parse_arguments()
    repeat_sequences = sorted(list(set([s.lower() for s in args.repeat_seq])))
    log('loading:', args.input)
    reference_genome = load_reference_genome(args.input)
    comments = [
        os.path.basename(__file__),
        'input: {}'.format(args.input),
        'min_length: {}'.format(args.min_length),
        'repeat_seq: {}'.format(', '.join(args.repeat_seq))
    ]
    log('writing:', args.output)
    with open(args.output, 'w') as fh:
        for comment in comments:
            fh.write('## {}\n'.format(comment))
        fh.write('chr\tstart\tend\tname\n')
        visited = set()
        for chrom, seq in sorted(reference_genome.items()):
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            seq = str(seq.seq).lower()
            if seq in visited:
                continue
            else:
                visited.add(seq)
            spans = []
            for repseq in repeat_sequences:
                log('finding {}_repeat (min_length: {}), for chr{} (length: {})'.format(repseq, args.min_length, chrom, len(seq)))
                index = 0
                while index < len(seq):
                    next_n = seq.find(repseq, index)
                    if next_n < 0:
                        break
                    index = next_n
                    while index + len(repseq) <= len(seq) and seq[index:index + len(repseq)] == repseq:
                        index += len(repseq)
                    span = BioInterval(chrom, next_n + 1, index, name='repeat_{}'.format(repseq))
                    if len(span) >= args.min_length and len(span) >= 2 * len(repseq):
                        spans.append(span)
            log('found', len(spans), 'spans', time_stamp=False)
            for span in spans:
                fh.write('{}\t{}\t{}\t{}\n'.format(span.reference_object, span.start, span.end, span.name))

if __name__ == '__main__':
    main()
