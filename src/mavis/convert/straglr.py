from typing import Dict

from ..constants import COLUMNS, SVTYPE


def convert_row(row: Dict) -> Dict:
    """
    Converts the fields from the original STRAGLR BED output into MAVIS definitions of an SV
    Since STRAGLR defines regions where short tandem repeats exist we make the definitions here fairly
    non-specific

    See their github page for more details: https://github.com/bcgsc/straglr

    BED Columns
    - chrom: chromosome name
    - start: start coordinate of locus
    - end: end coordinate of locus
    - repeat_unit: repeat motif
    - allele<N>.size: where N={1,2,3...} depending on --max_num_clusters e.g. N={1,2} if --max_num_clusters==2 (default)
    - allele<N>.copy_number
    - allele<N>.support
    """
    return {
        COLUMNS.break1_chromosome: row['chrom'],
        COLUMNS.break2_chromosome: row['chrom'],
        COLUMNS.break1_position_start: row['start'],
        COLUMNS.break1_position_end: row['end'],
        COLUMNS.break2_position_start: row['start'],
        COLUMNS.break2_position_end: row['end'],
        COLUMNS.untemplated_seq: None,
        COLUMNS.event_type: SVTYPE.INS,
    }
