from structural_variant.annotate import load_reference_genome
from structural_variant.interval import Interval
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output',
        help='path to the output file', required=True
    )
    parser.add_argument(
        '-n', '--input', default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
        help='path to the input file'
    )
    g = parser.add_mutually_exclusive_group()
    g.add_argument('--repeats', nargs=2, metavar=('<min size>', '<max size>'), type=int)
    g.add_argument('--nspans', nargs=1, metavar='<min length>', type=int)
    args = parser.parse_args()
    if not (args.repeats or args.nspans):
        print('error: --repeats exclusive or --nspans is required')
        parser.print_help()
        exit(1)
    return args

class Repeat(Interval):
    def __init__(self, chr, start, end, seq):
        Interval.__init__(self, start, end)
        self.seq = seq
        self.chr = chr

    def __hash__(self):
        return hash((self.chr, self.start, self.end, self.seq))


def main():
    args = parse_arguments()
    genome = load_reference_genome(args.input)
    print('loaded genome', args.input)
    with open(args.output, 'w') as fh:
        fh.write('chr\tstart\tend\tname\n')
        print('opened', args.output, 'for writing')
        if args.nspans:
            for chr in ['1']:
                print('finding nspans for chr{}'.format(chr))
                intervals = []
                last_was_N = False
                for i, c in enumerate(genome[chr].seq):
                    if not last_was_N:
                        if c == 'N':
                            last_was_N = True
                            intervals.append(Interval(i, i))
                    else: # last thing we saw was an N
                        if c == 'N':
                            intervals[-1].end = i
                        else:
                            last_was_N = False
                for i in intervals:
                    if len(i) >= args.nspans[0]:
                        fh.write('{}\t{}\t{}\tnspan\t{}\n'.format(chr, i.start, i.end, genome[chr].seq[i.start:i.end + 1]))
                        print(chr, i)
        else: # repeats
            repeats = []
            for chr in genome:
                l = len(genome[chr])
                for repeat_size in range(args.repeats[0], args.repeats[1] + 1):
                    for start in range(0, repeat_size):
                        pos = start
                        last = None
                        while pos + repeat_size <= l:
                            mer = Repeat(chr, pos, pos + repeat_size - 1, genome[chr].seq[pos:pos + repeat_size])
                            if last is not None:
                                if last.seq == mer.seq:
                                    last.end = mer.end
                                else:
                                    if len(last) > repeat_size:
                                        repeats.append(last)
                                        print(last.chr, last, last.seq)
                                    last = mer
                            else:
                                last = mer
                            pos += repeat_size
            for r in repeats:
                fh.write('{}\t{}\t{}\t{}\n'.format(r.chr, r.start, r.end, r.seq))


main()
