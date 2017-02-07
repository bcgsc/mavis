from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.constants import COLUMNS, PROTOCOL
import argparse
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='input file to check')
    args = parser.parse_args()
    print('loading:', args.input)
    bpps = read_bpp_from_input_file(
        args.input,
        validate={
            COLUMNS.tools: '^(\S+_v?\d+\.\d+\.\d+)(;\S+_v?\d+\.\d+\.\d+)*$',
            COLUMNS.library: '^[\w-]+$'
        },
        _in={
            COLUMNS.protocol: PROTOCOL
        }
    )
    
    print('OK! no errors were detected')
