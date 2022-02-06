import errno
import itertools
import logging
import os
import re
import time
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Set

import pandas as pd
from mavis_config import bash_expands
from shortuuid import uuid

from .breakpoint import Breakpoint, BreakpointPair
from .constants import (
    COLUMNS,
    FLOAT_COLUMNS,
    INTEGER_COLUMNS,
    ORIENT,
    PROTOCOL,
    STRAND,
    SUMMARY_LIST_COLUMNS,
    SVTYPE,
    sort_columns,
)
from .error import InvalidRearrangement
from .interval import Interval

if TYPE_CHECKING:
    from mavis.annotate.base import BioInterval

ENV_VAR_PREFIX = 'MAVIS_'

logger = logging.getLogger('mavis')


def filepath(path):
    try:
        file_list = bash_expands(path)
    except FileNotFoundError:
        raise TypeError('File does not exist', path)
    else:
        if not file_list:
            raise TypeError('File not found', path)
        elif len(file_list) > 1:
            raise TypeError('File pattern match multiple files and expected only one', path)
    return file_list[0]


class NullableType:
    def __init__(self, callback_func):
        self.callback_func = callback_func

    def __call__(self, item):
        if str(item).lower() == 'none':
            return None
        else:
            return self.callback_func(item)


def cast_null(input_value):
    value = str(input_value).lower()
    if value in ['none', 'null']:
        return None
    raise TypeError('casting to null/None failed', input_value)


def cast_boolean(input_value):
    value = str(input_value).lower()
    if value in ['t', 'true', '1', 'y', 'yes', '+']:
        return True
    elif value in ['f', 'false', '0', 'n', 'no', '-']:
        return False
    raise TypeError('casting to boolean failed', input_value)


def cast(value, cast_func):
    """
    cast a value to a given type

    Example:
        >>> cast('1', int)
        1
    """
    if cast_func == bool:
        value = cast_boolean(value)
    else:
        value = cast_func(value)
    return value


def soft_cast(value, cast_type):
    """
    cast a value to a given type, if the cast fails, cast to null

    Example:
        >>> cast(None, int)
        None
        >>> cast('', int)
        None
    """
    try:
        return cast(value, cast_type)
    except (TypeError, ValueError):
        pass
    return cast_null(value)


def log_arguments(args):
    """
    output the arguments to the console

    Args:
        args (Namespace): the namespace to print arguments for
    """
    logger.info('arguments')

    indent = ' '

    for arg, val in sorted(args.__dict__.items()):
        if isinstance(val, list):
            if len(val) <= 1:
                logger.info(f'{indent}{arg} = {val}')
                continue
            logger.info(f'{indent}{arg} = [')
            for v in val:
                logger.info(f'{indent * 2}{repr(v)}')
            logger.info(f'{indent}]')
        elif any([isinstance(val, typ) for typ in [str, int, float, bool, tuple]]) or val is None:
            logger.info(f'{indent}{arg}= {repr(val)}')
        else:
            logger.info(f'{arg} = {object.__repr__(val)}')


def mkdirp(dirname):
    """
    Make a directory or path of directories. Suppresses the error that is normally raised when the directory already exists
    """
    logger.info(f"creating output directory: '{dirname}'")
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise exc
    return dirname


def filter_on_overlap(
    bpps: List[BreakpointPair], regions_by_reference_name: Dict[str, List['BioInterval']]
):
    """
    filter a set of breakpoint pairs based on overlap with a set of genomic regions

    Args:
        bpps: list of breakpoint pairs to be filtered
        regions_by_reference_name: regions to filter against
    """
    logger.info(f'filtering from {len(bpps)} using overlaps with regions filter')
    failed = []
    passed = []
    for bpp in bpps:
        overlaps = False
        for r in regions_by_reference_name.get(bpp.break1.chr, []):
            if Interval.overlaps(r, bpp.break1):
                overlaps = True
                bpp.data[COLUMNS.filter_comment] = 'overlapped masked region: ' + str(r)
                break
        for r in regions_by_reference_name.get(bpp.break2.chr, []):
            if overlaps:
                break
            if Interval.overlaps(r, bpp.break2):
                overlaps = True
                bpp.data[COLUMNS.filter_comment] = 'overlapped masked region: ' + str(r)
        if overlaps:
            failed.append(bpp)
        else:
            passed.append(bpp)
    logger.info(f'filtered from {len(bpps)} down to {len(passed)} (removed {len(failed)})')
    return passed, failed


def read_inputs(
    inputs: List[str], required_columns: List[str] = [], **kwargs
) -> List[BreakpointPair]:
    bpps = []

    for finput in bash_expands(*inputs):
        logger.info(f'loading: {finput}')
        bpps.extend(
            read_bpp_from_input_file(
                finput, required_columns=[COLUMNS.protocol, *required_columns], **kwargs
            )
        )
    logger.info(f'loaded {len(bpps)} breakpoint pairs')
    return bpps


def output_tabbed_file(bpps: List[BreakpointPair], filename: str, header=None):
    if header is None:
        custom_header = False
        header = set()
    else:
        custom_header = True
    rows = []
    for row in bpps:
        if not isinstance(row, dict):
            row = row.flatten()
        rows.append(row)
        if not custom_header:
            header.update(row.keys())  # type: ignore
    header = sort_columns(header)
    logger.info(f'writing: {filename}')
    df = pd.DataFrame.from_records(rows, columns=header)
    df = df.fillna('None')
    df.to_csv(filename, columns=header, index=False, sep='\t')


def write_bed_file(filename, bed_rows):
    logger.info(f'writing: {filename}')
    with open(filename, 'w') as fh:
        for bed in bed_rows:
            fh.write('\t'.join([str(c) for c in bed]) + '\n')


def get_connected_components(adj_matrix):
    """
    for a dictionary representing an adjacency matrix of undirected edges returns the connected components
    """
    nodes_visited = set()
    components = []
    for node in adj_matrix:
        if node in nodes_visited:
            continue
        component = {node}
        unvisited = adj_matrix[node]
        while unvisited:
            curr = unvisited.pop()
            component.add(curr)
            nodes_visited.add(curr)
            unvisited.update(adj_matrix.get(curr, set()))
            unvisited.difference_update(nodes_visited)
        components.append(component)
    return components


def generate_complete_stamp(
    output_dir: str, prefix: str = 'MAVIS.', start_time: Optional[int] = None
) -> str:
    """
    writes a complete stamp, optionally including the run time if start_time is given

    Args:
        output_dir: path to the output dir the stamp should be written in
        prefix: prefix for the stamp name
        start_time: the start time

    Return:
        path to the complete stamp

    Example:
        >>> generate_complete_stamp('some_output_dir')
        'some_output_dir/MAVIS.COMPLETE'
    """
    stamp = os.path.join(output_dir, str(prefix) + 'COMPLETE')
    logger.info(f'complete: {stamp}')
    with open(stamp, 'w') as fh:
        if start_time is not None:
            duration = int(time.time()) - start_time
            hours = duration - duration % 3600
            minutes = duration - hours - (duration - hours) % 60
            seconds = duration - hours - minutes
            fh.write(
                'run time (hh/mm/ss): {}:{:02d}:{:02d}\n'.format(
                    hours // 3600, minutes // 60, seconds
                )
            )
            fh.write('run time (s): {}\n'.format(duration))
    return stamp


def filter_uninformative(annotations_by_chr, breakpoint_pairs, max_proximity=5000):
    result = []
    filtered = []
    for bpp in breakpoint_pairs:
        # loop over the annotations
        overlaps_gene = False
        window1 = Interval(bpp.break1.start - max_proximity, bpp.break1.end + max_proximity)
        window2 = Interval(bpp.break2.start - max_proximity, bpp.break2.end + max_proximity)
        for gene in annotations_by_chr.get(bpp.break1.chr, []):
            if Interval.overlaps(gene, window1):
                overlaps_gene = True
                break
        for gene in annotations_by_chr.get(bpp.break2.chr, []):
            if Interval.overlaps(gene, window2):
                overlaps_gene = True
                break
        if overlaps_gene:
            result.append(bpp)
        else:
            filtered.append(bpp)
    return result, filtered


def read_bpp_from_input_file(
    filename: str,
    expand_orient: bool = False,
    expand_strand: bool = False,
    expand_svtype: bool = False,
    integer_columns: Set[str] = INTEGER_COLUMNS,
    float_columns: Set[str] = FLOAT_COLUMNS,
    required_columns: Set[str] = set(),
    add_default: Dict[str, Any] = {},
    summary: bool = False,
    apply: Dict[str, Callable] = {},
    overwrite: Dict[str, Any] = {},
) -> List[BreakpointPair]:
    """
    reads a file using the tab module. Each row is converted to a breakpoint pair and
    other column data is stored in the data attribute

    Args:
        filename: path to the input file
        expand_ns: expand not specified orient/strand settings to all specific version (for strand this is only applied if the bam itself is stranded)
        explicit_strand: used to stop unstranded breakpoint pairs from losing input strand information
        summary: the input is post-summary so some float/int columns have been merged and delimited with semi-colons
        overwrite: set column values for all breakpoints, if the column exists overwrite its current value

    Returns:
        a list of pairs
    """

    def soft_null_cast(value):
        try:
            cast_null(value)
        except TypeError:
            return value

    if summary:
        integer_columns = integer_columns - SUMMARY_LIST_COLUMNS
        float_columns = float_columns - SUMMARY_LIST_COLUMNS

    try:
        df = pd.read_csv(
            filename,
            dtype={
                **{col: pd.Int64Dtype() for col in integer_columns},
                **{col: float for col in float_columns},
                **{
                    col: str
                    for col in COLUMNS.keys()
                    if col not in (float_columns | integer_columns)
                },
            },
            sep='\t',
            comment='#',
            na_values=['None', 'none', 'N/A', 'n/a', 'null', 'NULL', 'Null', 'nan', '<NA>', 'NaN'],
        )
        df = df.where(pd.notnull(df), None)
    except pd.errors.EmptyDataError:
        return []

    for col in required_columns:
        if col not in df and col not in add_default:
            raise KeyError(f'missing required column: {col}')

    if COLUMNS.opposing_strands in df:
        df[COLUMNS.opposing_strands] = df[COLUMNS.opposing_strands].apply(
            lambda x: None if x == '?' else soft_cast(x, cast_type=bool)
        )
    else:
        df[COLUMNS.opposing_strands] = None

    if COLUMNS.stranded in df:
        df[COLUMNS.stranded] = df[COLUMNS.stranded].apply(cast_boolean)
    else:
        df[COLUMNS.stranded] = None

    if COLUMNS.untemplated_seq in df:
        df[COLUMNS.untemplated_seq] = df[COLUMNS.untemplated_seq].apply(soft_null_cast)
    else:
        df[COLUMNS.untemplated_seq] = None

    for col in [COLUMNS.break1_chromosome, COLUMNS.break2_chromosome]:
        df[col] = df[col].apply(lambda v: re.sub(r'^chr', '', v))

    if COLUMNS.tracking_id not in df:
        df[COLUMNS.tracking_id] = ''
    else:
        df[COLUMNS.tracking_id] = df[COLUMNS.tracking_id].fillna(str(uuid()))

    # add default values
    for col, default_value in add_default.items():
        if col in df:
            df[col] = df[col].fillna(default_value)
        else:
            df[col] = default_value

    # run the custom functions
    for col, func in apply.items():
        df[col] = df[col].apply(func)

    # set overwriting defaults
    for col, value in overwrite.items():
        df[col] = value

    # enforce controlled vocabulary
    for vocab, cols in [
        (ORIENT, [COLUMNS.break1_orientation, COLUMNS.break2_orientation]),
        (STRAND, [COLUMNS.break1_strand, COLUMNS.break2_strand]),
        (PROTOCOL, [COLUMNS.protocol]),
    ]:
        for col in cols:
            if col in df:
                df[col].apply(lambda c: vocab.enforce(c))  # type: ignore
            elif hasattr(vocab, 'NS'):
                df[col] = vocab.NS  # type: ignore

    def validate_pipeline_id(value):
        if not re.match(r'^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$', value):
            raise AssertionError(
                'All mavis pipeline step ids must satisfy the regex:',
                '^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$',
                value,
            )

    for col in [COLUMNS.cluster_id, COLUMNS.annotation_id, COLUMNS.validation_id]:
        if col in df:
            try:
                df[col].apply(validate_pipeline_id)
            except AssertionError as err:
                raise AssertionError(f'error in column ({col}): {err}')

    rows = df.where(df.notnull(), None).to_dict('records')
    non_data_columns = {
        COLUMNS.break1_chromosome,
        COLUMNS.break1_position_start,
        COLUMNS.break1_position_end,
        COLUMNS.break1_strand,
        COLUMNS.break1_orientation,
        COLUMNS.break2_chromosome,
        COLUMNS.break2_position_start,
        COLUMNS.break2_position_end,
        COLUMNS.break2_strand,
        COLUMNS.break2_orientation,
        COLUMNS.stranded,
        COLUMNS.opposing_strands,
        COLUMNS.untemplated_seq,
    }
    pairs: List[BreakpointPair] = []

    for line_index, row in enumerate(rows):
        row['line_no'] = line_index + 1

        if '_index' in row:
            del row['_index']
        for attr, val in row.items():
            row[attr] = soft_null_cast(val)

        stranded = row[COLUMNS.stranded]

        strand1 = row[COLUMNS.break1_strand] if stranded else STRAND.NS
        strand2 = row[COLUMNS.break2_strand] if stranded else STRAND.NS

        temp = []
        expand_strand = stranded and expand_strand
        event_type = [None]
        if not pd.isnull(row.get(COLUMNS.event_type)):
            try:
                event_type = row[COLUMNS.event_type].split(';')
                for putative_event_type in event_type:
                    SVTYPE.enforce(putative_event_type)
            except KeyError:
                pass

        for orient1, orient2, strand1, strand2, putative_event_type in itertools.product(
            ORIENT.expand(row[COLUMNS.break1_orientation])
            if expand_orient
            else [row[COLUMNS.break1_orientation]],
            ORIENT.expand(row[COLUMNS.break2_orientation])
            if expand_orient
            else [row[COLUMNS.break2_orientation]],
            STRAND.expand(strand1) if expand_strand and stranded else [strand1],
            STRAND.expand(strand2) if expand_strand and stranded else [strand2],
            event_type,
        ):
            try:
                break1 = Breakpoint(
                    row[COLUMNS.break1_chromosome],
                    row[COLUMNS.break1_position_start],
                    row[COLUMNS.break1_position_end],
                    strand=strand1,
                    orient=orient1,
                )
                break2 = Breakpoint(
                    row[COLUMNS.break2_chromosome],
                    row[COLUMNS.break2_position_start],
                    row[COLUMNS.break2_position_end],
                    strand=strand2,
                    orient=orient2,
                )

                data = {k: v for k, v in row.items() if k not in non_data_columns}
                bpp = BreakpointPair(
                    break1,
                    break2,
                    opposing_strands=row[COLUMNS.opposing_strands],
                    untemplated_seq=row[COLUMNS.untemplated_seq],
                    stranded=row[COLUMNS.stranded],
                )
                bpp.data.update(data)
                if putative_event_type:
                    bpp.data[COLUMNS.event_type] = putative_event_type
                    if putative_event_type not in BreakpointPair.classify(bpp):
                        raise InvalidRearrangement(
                            'error: expected one of',
                            BreakpointPair.classify(bpp),
                            'but found',
                            putative_event_type,
                            str(bpp),
                            row,
                        )
                if expand_svtype and not putative_event_type:
                    for svtype in BreakpointPair.classify(
                        bpp, distance=lambda x, y: Interval(y - x)
                    ):
                        new_bpp = bpp.copy()
                        new_bpp.data[COLUMNS.event_type] = svtype
                        temp.append(new_bpp)
                else:
                    temp.append(bpp)
            except InvalidRearrangement as err:
                if not any([expand_strand, expand_svtype, expand_orient]):
                    raise err
            except AssertionError as err:
                if not expand_strand:
                    raise err
        if not temp:
            raise InvalidRearrangement('could not produce a valid rearrangement', row)
        else:
            pairs.extend(temp)
    return pairs
