import re

import pandas as pd


def convert_file(input_file):
    bam_to_lib = {}

    # read comments
    with open(input_file, 'r') as fh:
        # comments in breakdancer are marked with a single # so they need to be discarded before reading
        lines = fh.readlines()
        line_index = 0
        while line_index < len(lines) and lines[line_index].startswith('#'):
            metadata_match = re.match(r'^#(\S+)\t.*\tlibrary:(\S+)\t.*', lines[line_index])
            if metadata_match:
                bam_to_lib[metadata_match.group(1)] = metadata_match.group(2)
            line_index += 1
        header = [c.strip() for c in re.sub(r'^#', '', lines[line_index - 1]).split('\t')]
    # read the main file
    df = pd.read_csv(
        input_file,
        names=header,
        sep='\t',
        comment='#',
        dtype={
            'num_Reads_lib': str,
            'Pos1': int,
            'Pos2': int,
            'Chr1': str,
            'Chr2': str,
            'Type': str,
        },
    )
    if 'num_Reads_lib' not in df:
        raise KeyError('missing required column: num_Reads_lib')

    for bam, lib in bam_to_lib.items():
        df['num_Reads_lib'] = df['num_Reads_lib'].str.replace(bam, lib)
    return df.to_dict('records')
