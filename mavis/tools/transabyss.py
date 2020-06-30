import re

from ..constants import COLUMNS

from .constants import SUPPORTED_TOOL, TRACKING_COLUMN


def convert_row(row):
    """
    transforms the transabyss output into the common format for expansion.
    Maps the input column names to column names which MAVIS can read
    """
    std_row = {}
    if TRACKING_COLUMN not in row:
        std_row[TRACKING_COLUMN] = '{}-{}'.format(SUPPORTED_TOOL.TA, row['id'])

    std_row[COLUMNS.event_type] = row.get('rearrangement', row['type'])
    for retained_column in ['genes', 'gene']:
        if retained_column in row:
            std_row['{}_{}'.format(SUPPORTED_TOOL.TA, retained_column)] = row[retained_column]
    if std_row[COLUMNS.event_type] in ['LSR', 'translocation']:
        del std_row[COLUMNS.event_type]
    if 'breakpoint' in row:
        std_row[COLUMNS.break1_orientation], std_row[COLUMNS.break2_orientation] = row[
            'orientations'
        ].split(',')
        match = re.match(
            r'^(?P<chr1>[^:]+):(?P<pos1_start>\d+)\|(?P<chr2>[^:]+):(?P<pos2_start>\d+)$',
            row['breakpoint'],
        )
        if not match:
            raise OSError(
                'file format error: the breakpoint column did not satisfy the expected pattern', row
            )
        for group, col in zip(
            ['chr1', 'pos1_start', 'chr2', 'pos2_start'],
            [
                COLUMNS.break1_chromosome,
                COLUMNS.break1_position_start,
                COLUMNS.break2_chromosome,
                COLUMNS.break2_position_start,
            ],
        ):
            std_row[col] = match[group]
    else:
        std_row.update(
            {
                COLUMNS.break1_chromosome: row['chr'],
                COLUMNS.break1_position_start: int(row['chr_start']),
                COLUMNS.break2_position_start: int(row['chr_end']),
            }
        )
        if std_row[COLUMNS.event_type] == 'del':
            std_row[COLUMNS.break1_position_start] -= 1
            std_row[COLUMNS.break2_position_start] += 1
        elif std_row[COLUMNS.event_type] == 'ins':
            std_row[COLUMNS.break2_position_start] += 1

        # add the untemplated sequence where appropriate
        if std_row[COLUMNS.event_type] == 'del':
            assert row['alt'] == 'na'
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row[COLUMNS.event_type] in ['dup', 'ITD']:
            length = (
                std_row[COLUMNS.break2_position_start] - std_row[COLUMNS.break1_position_start] + 1
            )
            if len(row['alt']) != length:
                raise AssertionError(
                    'expected alternate sequence to be equal to the length of the event',
                    len(row['alt']),
                    length,
                    row,
                    std_row,
                )
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row[COLUMNS.event_type] == 'ins':
            std_row[COLUMNS.untemplated_seq] = row['alt'].upper()
        else:
            raise NotImplementedError('unexpected indel type', std_row[COLUMNS.event_type])
    return std_row
