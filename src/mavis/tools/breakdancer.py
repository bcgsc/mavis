import re

import tab
from argparse import Namespace


def convert_file(input_file):
    bam_to_lib = {}
    with open(input_file, 'r') as fh:
        # comments in breakdancer are marked with a single # so they need to be discarded before reading
        lines = fh.readlines()
        header = 0
        while header < len(lines) and lines[header].startswith('#'):
            metadata_match = re.match(r'^#(\S+)\t.*\tlibrary:(\S+)\t.*', lines[header])
            if metadata_match:
                bam_to_lib[metadata_match.group(1)] = metadata_match.group(2)
            header += 1
        lines = lines[header - 1 :]
        input_file = Namespace(readlines=lambda: lines)
    header, rows = tab.read_file(input_file, allow_short=True, require=['num_Reads_lib'])
    for row in rows:
        for bam, lib in bam_to_lib.items():
            row['num_Reads_lib'] = row['num_Reads_lib'].replace(bam, lib)
    return rows
