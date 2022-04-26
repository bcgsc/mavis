#!/usr/bin/env python
import argparse
import logging
import itertools
from typing import Dict, Tuple, List
from mavis.constants import COLUMNS, ORIENT, STRAND, SVTYPE
import pandas as pd
from mavis.error import InvalidRearrangement

from mavis.annotate.file_io import parse_annotations_json
from mavis.annotate.base import ReferenceName
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.convert import TRACKING_COLUMN
from mavis.util import output_tabbed_file

"""
Converts existing dgv tab file to a MAVIS format. 
DGV files can be obtained from: http://dgv.tcag.ca/dgv/app/downloads?ref=GRCh37/hg19
"""


def _convert_tool_row(
    row: Dict, file_type: str, stranded: bool, assume_no_untemplated: bool = True
) -> List[BreakpointPair]:
    """
    converts a row parsed from an input file to the appropriate column names for it to be converted to MAVIS style row
    """
    std_row = {}
    try:
        std_row[TRACKING_COLUMN] = row.get(TRACKING_COLUMN, '')
    except AttributeError:
        try:
            std_row[TRACKING_COLUMN] = getattr(row, TRACKING_COLUMN)
        except AttributeError:
            pass
    std_row[COLUMNS.break1_orientation] = std_row[COLUMNS.break2_orientation] = ORIENT.NS
    std_row[COLUMNS.break1_strand] = std_row[COLUMNS.break2_strand] = STRAND.NS
    result = []

    # convert the specified file type to a standard format
    std_row.update(
        {
            COLUMNS.break1_chromosome: row['chr'],
            COLUMNS.break2_chromosome: row['chr'],
            COLUMNS.break1_position_start: row['start'],
            COLUMNS.break1_position_end: row['start'],
            COLUMNS.break2_position_start: row['end'],
            COLUMNS.break2_position_end: row['end'],
            COLUMNS.untemplated_seq: None,
            COLUMNS.tracking_id: row['variantaccession'],
        }
    )

    if stranded:
        std_row[COLUMNS.break1_strand] = STRAND.expand(std_row[COLUMNS.break1_strand])
        std_row[COLUMNS.break2_strand] = STRAND.expand(std_row[COLUMNS.break2_strand])
    else:
        std_row[COLUMNS.break1_strand] = [STRAND.NS]
        std_row[COLUMNS.break2_strand] = [STRAND.NS]

    if not std_row.get(TRACKING_COLUMN, None):
        std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, std_row.get('id', uuid()))

    combinations = list(
        itertools.product(
            ORIENT.expand(std_row[COLUMNS.break1_orientation]),
            ORIENT.expand(std_row[COLUMNS.break2_orientation]),
            std_row[COLUMNS.break1_strand],
            std_row[COLUMNS.break2_strand],
            TOOL_SVTYPE_MAPPING[std_row[COLUMNS.event_type]]
            if COLUMNS.event_type in std_row
            else [None],
            [True, False]
            if std_row.get(COLUMNS.opposing_strands, None) is None
            else [std_row[COLUMNS.opposing_strands]],
        )
    )
    # # add the product of all uncertainties as breakpoint pairs
    for orient1, orient2, strand1, strand2, event_type, oppose in combinations:
        try:
            untemplated_seq = std_row.get(COLUMNS.untemplated_seq, None)
            if assume_no_untemplated and event_type != SVTYPE.INS and not untemplated_seq:
                untemplated_seq = ''
            break1 = Breakpoint(
                std_row[COLUMNS.break1_chromosome],
                std_row[COLUMNS.break1_position_start],
                std_row.get(COLUMNS.break1_position_end, std_row[COLUMNS.break1_position_start]),
                orient=orient1,
                strand=strand1,
            )
            break2 = Breakpoint(
                std_row.get(COLUMNS.break2_chromosome, std_row[COLUMNS.break1_chromosome]),
                std_row[COLUMNS.break2_position_start],
                std_row.get(COLUMNS.break2_position_end, std_row[COLUMNS.break2_position_start]),
                orient=orient2,
                strand=strand2,
            )
            if (
                len(break1) == 1
                and len(break2) == 1
                and event_type == SVTYPE.DEL
                and abs(break1.start - break2.start) < 2
            ):
                break1 = Breakpoint(
                    break1.chr,
                    break1.start - 1,
                    break1.end - 1,
                    orient=break1.orient,
                    strand=break1.strand,
                )
                break2 = Breakpoint(
                    break2.chr,
                    break2.start + 1,
                    break2.end + 1,
                    orient=break2.orient,
                    strand=break2.strand,
                )
            bpp = BreakpointPair(
                break1,
                break2,
                opposing_strands=oppose,
                untemplated_seq=untemplated_seq,
                event_type=event_type,
                stranded=stranded,
                **{COLUMNS.tools: file_type, COLUMNS.tracking_id: std_row[COLUMNS.tracking_id]},
            )

            for col, value in std_row.items():
                if col not in COLUMNS and col not in bpp.data:
                    bpp.data[col] = value
            if not event_type or event_type in BreakpointPair.classify(bpp):
                result.append(bpp)
        except (InvalidRearrangement, AssertionError):
            pass

    return result


def convert_to_dictionary(filepaths: str) -> Dict:
    """
    reads a file of regions. The expect input format for the file is tab-delimited and
    the header should contain the following columns

    - chr: the chromosome
    - start: start of the region, 1-based inclusive
    - end: end of the region, 1-based inclusive
    - name: the name/label of the region

    For example:
        #chr    start       end         name
        chr20   25600000    27500000    centromere

    Args:
        filepath: path to the input tab-delimited file
    Returns:
        Dictionary by chromosome name with values of lists of regions on the chromosome
    """
    col_list = ['chr', 'start', 'end', 'variantaccession']
    result = []
    rows = None
    df = pd.read_csv(
        filepaths,
        sep='\t',
        dtype={'chr': str, 'start': int, 'end': int, 'variantaccession': str},
        usecols=col_list,
    )
    rows = df.where(df.notnull(), None).to_dict('records')
    if rows:
        for row in rows:
            try:
                std_rows = _convert_tool_row(row, "bed", False)
            except Exception as err:
                raise err
            else:
                result.extend(std_rows)
    return result


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


def main():

    parser = argparse.ArgumentParser(
        description="Convert the dgv file format into a MAVIS file.",
        add_help=False,
    )
    req_parser = parser.add_argument_group("required arguments")
    req_parser.add_argument("-n", "--input", required=True, help="input DGV file")
    req_parser.add_argument("-o", "--output", required=True, help="output file name")
    args = parser.parse_args()
    output_tabbed_file(convert_to_dictionary(args.input), args.output)


if __name__ == "__main__":
    main()
