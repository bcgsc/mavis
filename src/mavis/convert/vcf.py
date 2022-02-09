import gzip
import logging
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd

try:
    # TypedDict added to typing package directly in later versions
    from typing import TypedDict
except ImportError:
    from typing_extensions import TypedDict

from ..constants import COLUMNS, ORIENT, SVTYPE
from ..util import logger

PANDAS_DEFAULT_NA_VALUES = [
    '-1.#IND',
    '1.#QNAN',
    '1.#IND',
    '-1.#QNAN',
    '#N/A',
    'N/A',
    'NA',
    '#NA',
    'NULL',
    'NaN',
    '-NaN',
    'nan',
    '-nan',
]


class VcfInfoType(TypedDict, total=False):
    SVTYPE: str
    CHR2: str
    CIPOS: Tuple[int, int]
    CIEND: Tuple[int, int]
    CT: str
    END: Optional[int]
    PRECISE: bool


@dataclass
class VcfRecordType:
    id: str
    pos: int
    chrom: str
    alts: List[Optional[str]]
    info: VcfInfoType
    ref: str

    @property
    def stop(self) -> Optional[int]:
        return self.info.get('END', self.pos)


def parse_bnd_alt(alt: str) -> Tuple[str, int, str, str, str, str]:
    """
    parses the alt statement from vcf files using the specification in vcf 4.2/4.2.

    Assumes that the reference base is always the outermost base (this is based on the spec and also manta results as
    the spec was missing some cases)

    r = reference base/seq
    u = untemplated sequence/alternate sequence
    p = chromosome:position

    | alt format   | orients |
    | ------------ | ------- |
    | ru[p[        | LR      |
    | [p[ur        | RR      |
    | ]p]ur        | RL      |
    | ru]p]        | LL      |
    """
    # ru[p[
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)\[(?P<chr>[^:]+):(?P<pos>\d+)\[$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.LEFT,
            ORIENT.RIGHT,
            match.group('ref'),
            match.group('useq'),
        )
    # [p[ur
    match = re.match(r'^\[(?P<chr>[^:]+):(?P<pos>\d+)\[(?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.RIGHT,
            ORIENT.RIGHT,
            match.group('ref'),
            match.group('useq'),
        )
    # ]p]ur
    match = re.match(r'^\](?P<chr>[^:]+):(?P<pos>\d+)\](?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.RIGHT,
            ORIENT.LEFT,
            match.group('ref'),
            match.group('useq'),
        )
    # ru]p]
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)\](?P<chr>[^:]+):(?P<pos>\d+)\]$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.LEFT,
            ORIENT.LEFT,
            match.group('ref'),
            match.group('useq'),
        )
    else:
        raise NotImplementedError('alt specification in unexpected format', alt)


def convert_record(record: VcfRecordType) -> List[Dict]:
    """
    converts a vcf record

    Note:
        CT = connection type, If given this field will be used in determining the orientation at the breakpoints.
        From https://groups.google.com/forum/#!topic/delly-users/6Mq2juBraRY, we can expect certain CT types for
        certain event types
            - translocation/inverted translocation: 3to3, 3to5, 5to3, 5to5
            - inversion: 3to3, 5to5
            - deletion: 3to5
            - duplication: 5to3
    """
    records = []

    for alt in record.alts if record.alts else [None]:
        info = {}
        for key in record.info.keys():
            try:
                value = record.info[key]
            except UnicodeDecodeError as err:
                logger.warning(f'Ignoring invalid INFO field {key} with error: {err}')
            else:
                try:
                    value = value[0] if len(value) == 1 else value
                except TypeError:
                    pass  # anything non-tuple
            info[key] = value

        std_row = {}
        if record.id and record.id != 'N':  # to account for NovoBreak N in the ID field
            std_row['id'] = record.id

        if info.get('SVTYPE') == 'BND':
            chr2, end, orient1, orient2, ref, alt = parse_bnd_alt(alt)
            std_row[COLUMNS.break1_orientation] = orient1
            std_row[COLUMNS.break2_orientation] = orient2
            std_row[COLUMNS.untemplated_seq] = alt
            if record.ref != ref:
                raise AssertionError(
                    'Expected the ref specification in the vcf record to match the sequence '
                    'in the alt string: {} vs {}'.format(record.ref, ref)
                )
        else:
            chr2 = info.get('CHR2', record.chrom)
            end = record.stop
            if (
                alt
                and record.ref
                and re.match(r'^[A-Z]+$', alt)
                and re.match(r'^[A-Z]+', record.ref)
            ):
                std_row[COLUMNS.untemplated_seq] = alt[1:]
                size = len(alt) - len(record.ref)
                if size > 0:
                    std_row[COLUMNS.event_type] = SVTYPE.INS
                elif size < 0:
                    std_row[COLUMNS.event_type] = SVTYPE.DEL
        std_row.update({COLUMNS.break1_chromosome: record.chrom, COLUMNS.break2_chromosome: chr2})
        if info.get(
            'PRECISE', False
        ):  # DELLY CI only apply when split reads were not used to refine the breakpoint which is then flagged
            std_row.update(
                {
                    COLUMNS.break1_position_start: record.pos,
                    COLUMNS.break1_position_end: record.pos,
                    COLUMNS.break2_position_start: end,
                    COLUMNS.break2_position_end: end,
                }
            )
        else:
            std_row.update(
                {
                    COLUMNS.break1_position_start: max(
                        1, record.pos + info.get('CIPOS', (0, 0))[0]
                    ),
                    COLUMNS.break1_position_end: record.pos + info.get('CIPOS', (0, 0))[1],
                    COLUMNS.break2_position_start: max(1, end + info.get('CIEND', (0, 0))[0]),
                    COLUMNS.break2_position_end: end + info.get('CIEND', (0, 0))[1],
                }
            )
        if std_row['break1_position_end'] == 0 and std_row['break1_position_start'] == 1:
            # addresses cases where pos = 0 and telomeric BND alt syntax https://github.com/bcgsc/mavis/issues/294
            std_row.update({'break1_position_end': 1})
        if std_row['break2_position_end'] == 0 and std_row['break2_position_start'] == 1:
            std_row.update({'break2_position_end': 1})

        if 'SVTYPE' in info:
            std_row[COLUMNS.event_type] = info['SVTYPE']

        try:
            orient1, orient2 = info['CT'].split('to')
            connection_type = {'3': ORIENT.LEFT, '5': ORIENT.RIGHT, 'N': ORIENT.NS}
            std_row[COLUMNS.break1_orientation] = connection_type[orient1]
            std_row[COLUMNS.break2_orientation] = connection_type[orient2]
        except KeyError:
            pass
        std_row.update(
            {k: v for k, v in info.items() if k not in {'CHR2', 'SVTYPE', 'CIPOS', 'CIEND', 'CT'}}
        )
        records.append(std_row)
    return records


def convert_pandas_rows_to_variants(df: pd.DataFrame) -> List[VcfRecordType]:
    def parse_info(info_field):
        info = {}
        for pair in info_field.split(';'):
            if '=' in pair:
                key, value = pair.split('=', 1)
                info[key] = value
            else:
                info[pair] = True

        # convert info types
        for key in info:
            if key in {'CIPOS', 'CIEND'}:
                ci_start, ci_end = info[key].split(',')
                info[key] = (int(ci_start), int(ci_end))
            elif key == 'END':
                info[key] = int(info[key])

        return info

    df['info'] = df['INFO'].apply(parse_info)
    df['alts'] = df['ALT'].apply(lambda a: a.split(','))

    rows = []
    for _, row in df.iterrows():

        rows.append(
            VcfRecordType(
                id=row['ID'],
                pos=row['POS'],
                info=VcfInfoType(row['info']),
                chrom=row['CHROM'],
                ref=row['REF'],
                alts=row['alts'],
            )
        )
    return rows


def pandas_vcf(input_file: str) -> Tuple[List[str], pd.DataFrame]:
    """
    Read a standard vcf file into a pandas dataframe
    """
    # read the comment/header information
    try:
        header_lines = []
        with open(input_file, 'r') as fh:
            line = '##'
            while line.startswith('##'):
                header_lines.append(line)
                line = fh.readline().strip()
            header_lines = header_lines[1:]
    except UnicodeDecodeError:
        header_lines = []
        with gzip.open(input_file, 'rt') as fh:
            line = '##'
            while line.startswith('##'):
                header_lines.append(line)
                line = fh.readline().strip()
            header_lines = header_lines[1:]
    # read the data
    df = pd.read_csv(
        input_file,
        sep='\t',
        skiprows=len(header_lines),
        dtype={
            'CHROM': str,
            'POS': int,
            'ID': str,
            'INFO': str,
            'FORMAT': str,
            'REF': str,
            'ALT': str,
        },
        na_values=PANDAS_DEFAULT_NA_VALUES + ['.'],
    )
    df = df.rename(columns={df.columns[0]: df.columns[0].replace('#', '')})
    required_columns = ['CHROM', 'INFO', 'POS', 'REF', 'ALT', 'ID']
    for col in required_columns:
        if col not in df.columns:
            raise KeyError(f'Missing required column: {col}')
    # convert the format fields using the header
    return header_lines, df


def convert_file(input_file: str) -> List[Dict]:
    """process a VCF file

    Args:
        input_file: the input file name

    Raises:
        err: [description]
    """
    rows = []

    _, data = pandas_vcf(input_file)

    for variant_record in convert_pandas_rows_to_variants(data):
        try:
            rows.extend(convert_record(variant_record))
        except NotImplementedError as err:
            logging.warning(str(err))
    return rows
