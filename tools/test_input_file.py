import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from mavis.breakpoint import read_bpp_from_input_file
from mavis.constants import COLUMNS, PROTOCOL, STRAND


def check_input_file(filename):
    print('loading:', args.input)
    bpps = read_bpp_from_input_file(
        filename,
        validate={
            COLUMNS.tools: '^(\S+_[^\s;]+)(;\S+_[^\s;]+)*$',
            COLUMNS.library: '^[\w-]+$'
        },
        in_={
            COLUMNS.protocol: PROTOCOL
        },
        add={COLUMNS.library: 'None'}
    )
    for bpp in bpps:
        if any([
            not bpp.stranded and bpp.break1.strand != STRAND.NS,
            not bpp.stranded and bpp.break2.strand != STRAND.NS
        ]):
            raise UserWarning('Error in input file. Cannot specify the strand if the pair is not stranded', bpp)
    print('loaded:', len(bpps), 'breakpoint pairs')
    print('OK! no errors were detected')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='input file to check')
    args = parser.parse_args()
    check_input_file(args.input) 
