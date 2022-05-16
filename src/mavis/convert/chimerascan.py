from typing import Dict

from ..constants import COLUMNS, ORIENT
from .constants import SUPPORTED_TOOL, TRACKING_COLUMN


def convert_row(row: Dict) -> Dict:
    """
    transforms the chimerscan output into the common format for expansion. Maps the input column
    names to column names which MAVIS can read
    """
    std_row = {}
    for retained_column in ['genes5p', 'genes3p']:
        if retained_column in row:
            std_row['{}_{}'.format(SUPPORTED_TOOL.CHIMERASCAN, retained_column)] = row[
                retained_column
            ]
    if TRACKING_COLUMN not in row:
        std_row[TRACKING_COLUMN] = '{}-{}'.format(
            SUPPORTED_TOOL.CHIMERASCAN, row['chimera_cluster_id']
        )

    std_row.update(
        {COLUMNS.break1_chromosome: row['chrom5p'], COLUMNS.break2_chromosome: row['chrom3p']}
    )
    if row['strand5p'] == '+':
        std_row[COLUMNS.break1_position_start] = row['end5p']
        std_row[COLUMNS.break1_orientation] = ORIENT.LEFT
    else:
        std_row[COLUMNS.break1_position_start] = row['start5p']
        std_row[COLUMNS.break1_orientation] = ORIENT.RIGHT
    if row['strand3p'] == '+':
        std_row[COLUMNS.break2_position_start] = row['start3p']
        std_row[COLUMNS.break2_orientation] = ORIENT.RIGHT
    else:
        std_row[COLUMNS.break2_position_start] = row['end3p']
        std_row[COLUMNS.break2_orientation] = ORIENT.LEFT
    std_row[COLUMNS.opposing_strands] = row['strand5p'] != row['strand3p']
    return std_row
