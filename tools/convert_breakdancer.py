import itertools
import os
import argparse
import TSV
from mavis.constants import COLUMNS, ORIENT, SVTYPE, PROTOCOL
from mavis.error import InvalidRearrangement
from mavis.util import output_tabbed_file
from mavis.breakpoint import BreakpointPair, Breakpoint
from mavis import __version__

__prog__ = os.path.basename(os.path.realpath(__file__))
default_version = '1.3.7'

SUPPORTED_EVENT_TYPES = {
    'DEL': [SVTYPE.DEL],
    'INS': [SVTYPE.INS],
    'ITX': [SVTYPE.DUP],
    'CTX': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'INV': [SVTYPE.INV]
}


def main():
    parser = argparse.ArgumentParser(
        description='Converts breakdancer file(s) to mavis-compatible tab file.',
        add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-n', '--input', required=True, help='Breakdancer *.max output', nargs='+')
    required.add_argument('-o', '--output', required=True, help='Full name of the converted output file')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='outputs version number')
    optional.add_argument(
        '--tool-version', default=default_version,
        help='the version of breakdander that was used in the analysis')

    args = parser.parse_args()

    for f in args.input:
        if not os.path.isfile(f):
            print("ERROR: Cannot find file: " + f)
            print("Exiting.")
            exit(1)

    events = []

    for bd_file in args.input:
        print("Now processing: " + bd_file)
        header, rows = TSV.read_file(bd_file, cast={'Type': lambda x: SUPPORTED_EVENT_TYPES[x]})

        for row in rows:
            data = {}
            data[COLUMNS.protocol] = PROTOCOL.GENOME
            data[COLUMNS.tools] = 'breakdancer_v' + args.tool_version
            data[COLUMNS.event_type] = row['Type']
            count = 0
            for orient1, orient2, opposing in itertools.product(
                    [ORIENT.LEFT, ORIENT.RIGHT], [ORIENT.LEFT, ORIENT.RIGHT], [True, False]):
                try:
                    bpp = BreakpointPair(
                        Breakpoint(row['Chr1'], row['Pos1'], orient=orient1),
                        Breakpoint(row['Chr2'], row['Pos2'], orient=orient2),
                        opposing_strands=opposing, stranded=False, data=data
                    )
                    if bpp.event_type in BreakpointPair.classify(bpp):
                        events.append(bpp)
                        count += 1
                except InvalidRearrangement:
                    pass
            if count < 1:
                raise UserWarning('did not find a single valid event combination for', row)

    output_tabbed_file(events, args.output)


if __name__ == '__main__':
    main()
