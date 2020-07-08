from ..constants import ORIENT


def convert_row(row):
    """
    transforms the starfusion output into the common format for expansion. Maps the input column
    names to column names that MAVIS can read
    """
    std_row = {}
    try:
        std_row['break1_chromosome'], b1_start, std_row['break1_strand'] = row[
            'LeftBreakpoint'
        ].split(':')
        std_row['break2_chromosome'], b2_start, std_row['break2_strand'] = row[
            'RightBreakpoint'
        ].split(':')
    except (ValueError, TypeError):
        raise AssertionError(
            'Could not parse the breakpoint from the starfusion row: {}, {}'.format(
                row['LeftBreakpoint'], row['RightBreakpoint']
            )
        )
    std_row['break1_position_start'] = std_row['break1_position_end'] = b1_start
    std_row['break2_position_start'] = std_row['break2_position_end'] = b2_start

    std_row['break1_orientation'] = std_row['break2_orientation'] = ORIENT.NS

    return std_row
