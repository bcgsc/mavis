from ..constants import COLUMNS, ORIENT, STRAND

from .constants import TRACKING_COLUMN, SUPPORTED_TOOL


def get_orient(string):
    if string == "downstream":
        return ORIENT.LEFT
    elif string == "upstream":
        return ORIENT.RIGHT
    return ORIENT.NS


def convert_row(row):
    """
    transforms the aribba output into the common format for expansion. Maps the input column
    names to column names which MAVIS can read
    """
    std_row = {}

    try:
        std_row[COLUMNS.break1_chromosome], b1_start = row["breakpoint1"].split(":")
        std_row[COLUMNS.break2_chromosome], b2_start = row["breakpoint2"].split(":")

        std_row[COLUMNS.break1_strand] = row["strand1(gene/fusion)"].split("/")[1]
        std_row[COLUMNS.break2_strand] = row["strand2(gene/fusion)"].split("/")[1]
        std_row[COLUMNS.event_type] = row["type"].split("/")[0]
        std_row[COLUMNS.break1_orientation] = get_orient(row["direction1"])
        std_row[COLUMNS.break2_orientation] = get_orient(row["direction2"])

        std_row[COLUMNS.break1_position_start] = std_row[COLUMNS.break1_position_end] = b1_start
        std_row[COLUMNS.break2_position_start] = std_row[COLUMNS.break2_position_end] = b2_start
    except (ValueError, TypeError):
        raise AssertionError(
            "Could not parse the breakpoint from the Arriba row: {}, {}".format(
                row["breakpoint1"], row["breakpoint2"]
            )
        )
    return std_row
