from typing import Dict

from ..constants import COLUMNS, SVTYPE


def convert_row(row: Dict) -> Dict:
    """
    Converts the fields from DGV output into MAVIS definitions of an SV

    Files can be obtained from: http://dgv.tcag.ca/dgv/app/downloads?ref=GRCh37/hg19

    Extracted BED Columns
    - variantaccession: id for tracking information
    - chrom: chromosome name
    - start: start coordinate of locus
    - end: end coordinate of locus
    """
    return {
        COLUMNS.break1_chromosome: row['chr'],
        COLUMNS.break2_chromosome: row['chr'],
        COLUMNS.break1_position_start: row['start'],
        COLUMNS.break1_position_end: row['start'],
        COLUMNS.break2_position_start: row['end'],
        COLUMNS.break2_position_end: row['end'],
        COLUMNS.untemplated_seq: None,
        COLUMNS.tracking_id: row['variantaccession'],
    }
