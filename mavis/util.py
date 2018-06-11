from argparse import Namespace
from datetime import datetime
import errno
from functools import partial
from glob import glob
import itertools
import os
import re
import time
import logging
import sys

from braceexpand import braceexpand
from tab import tab
from shortuuid import uuid

from .breakpoint import Breakpoint, BreakpointPair
from .constants import COLUMNS, ORIENT, PROTOCOL, sort_columns, STRAND, SVTYPE, MavisNamespace
from .error import InvalidRearrangement
from .interval import Interval

ENV_VAR_PREFIX = 'MAVIS_'


class Log:
    """
    wrapper aroung the builtin logging to make it more readable
    """
    def __init__(self, indent_str='  ', indent_level=0, level=logging.INFO):
        self.indent_str = indent_str
        self.indent_level = indent_level
        self.level = level

    def __call__(self, *pos, time_stamp=False, level=None, indent_level=0, **kwargs):
        if level is None and self.level is None:
            return
        elif self.level is not None:
            level = self.level

        stamp = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]') if time_stamp else ' ' * 21
        indent_prefix = self.indent_str * (self.indent_level + indent_level)
        message = '{} {}{}'.format(stamp, indent_prefix, ' '.join([str(p) for p in pos]))
        logging.log(level, message, **kwargs)

    def indent(self):
        return Log(self.indent_str, self.indent_level + 1, self.level)

    def dedent(self):
        return Log(self.indent_str, max(0, self.indent_level - 1), self.level)

    def __enter__(self):
        return self

    def __exit__(self, *pos):
        pass


LOG = Log()
DEVNULL = Log(level=None)


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


def cast(value, cast_func):
    """
    cast a value to a given type

    Example:
        >>> cast('1', int)
        1
    """
    if cast_func == bool:
        value = tab.cast_boolean(value)
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
    return tab.cast_null(value)


def get_env_variable(arg, default, cast_type=None):
    """
    Args:
        arg (str): the argument/variable name
    Returns:
        the setting from the environment variable if given, otherwise the default value
    """
    if cast_type is None:
        cast_type = type(default)
    name = ENV_VAR_PREFIX + str(arg).upper()
    result = os.environ.get(name, None)
    if result is not None:
        return cast(result, cast_type)
    return default


class WeakMavisNamespace(MavisNamespace):

    def is_env_overwritable(self, attr):
        return True


def bash_expands(*expressions):
    """
    expand a file glob expression, allowing bash-style brackets.

    Returns:
        list: a list of files

    Example:
        >>> bash_expands('./{test,doc}/*py')
        [...]
    """
    result = []
    for expression in expressions:
        eresult = []
        for name in braceexpand(expression):
            for fname in glob(name):
                eresult.append(fname)
        if not eresult:
            raise FileNotFoundError('The expression does not match any files', expression)
        result.extend(eresult)
    return [os.path.abspath(f) for f in result]


def log_arguments(args):
    """
    output the arguments to the console

    Args:
        args (Namespace): the namespace to print arguments for
    """
    LOG('arguments', time_stamp=True)
    with LOG.indent() as log:
        for arg, val in sorted(args.items()):
            if isinstance(val, list):
                if len(val) <= 1:
                    log(arg, '= {}'.format(val))
                    continue
                log(arg, '= [')
                for v in val:
                    log(repr(v), indent_level=1)
                log(']')
            elif any([isinstance(val, typ) for typ in [str, int, float, bool, tuple]]) or val is None:
                log(arg, '=', repr(val))
            else:
                log(arg, '=', object.__repr__(val))


def mkdirp(dirname):
    """
    Make a directory or path of directories. Suppresses the error that is normally raised when the directory already exists
    """
    LOG("creating output directory: '{}'".format(dirname))
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise exc
    return dirname


def filter_on_overlap(bpps, regions_by_reference_name):
    """
    filter a set of breakpoint pairs based on overlap with a set of genomic regions

    Args:
        bpps (:class:`list` of :class:`~mavis.breakpoint.BreakpointPair`): list of breakpoint pairs to be filtered
        regions_by_reference_name (:class:`dict` of :class:`list` of :class:`~mavis.annotate.base.BioInterval` by :class:`str`): regions to filter against
    """
    LOG('filtering from', len(bpps), 'using overlaps with regions filter')
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
    LOG('filtered from', len(bpps), 'down to', len(passed), '(removed {})'.format(len(failed)))
    return passed, failed


def read_inputs(inputs, **kwargs):
    bpps = []
    kwargs.setdefault('require', [])
    kwargs['require'] = list(set(kwargs['require'] + [COLUMNS.protocol]))
    kwargs.setdefault('in_', {})
    kwargs['in_'][COLUMNS.protocol] = PROTOCOL.values()
    for finput in bash_expands(*inputs):
        try:
            LOG('loading:', finput)
            bpps.extend(read_bpp_from_input_file(
                finput,
                **kwargs
            ))
        except tab.EmptyFileError:
            LOG('ignoring empty file:', finput)
    LOG('loaded', len(bpps), 'breakpoint pairs')
    return bpps


def output_tabbed_file(bpps, filename, header=None):
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
            header.update(row.keys())
    header = sort_columns(header)

    with open(filename, 'w') as fh:
        LOG('writing:', filename)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')


def write_bed_file(filename, bed_rows):
    LOG('writing:', filename)
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


def generate_complete_stamp(output_dir, log=DEVNULL, prefix='MAVIS.', start_time=None):
    """
    writes a complete stamp, optionally including the run time if start_time is given

    Args:
        output_dir (str): path to the output dir the stamp should be written in
        log (function): function to print logging messages to
        prefix (str): prefix for the stamp name
        start_time (int): the start time

    Return:
        str: path to the complete stamp

    Example:
        >>> generate_complete_stamp('some_output_dir')
        'some_output_dir/MAVIS.COMPLETE'
    """
    stamp = os.path.join(output_dir, str(prefix) + 'COMPLETE')
    log('complete:', stamp)
    with open(stamp, 'w') as fh:
        if start_time is not None:
            duration = int(time.time()) - start_time
            hours = duration - duration % 3600
            minutes = duration - hours - (duration - hours) % 60
            seconds = duration - hours - minutes
            fh.write('run time (hh/mm/ss): {}:{:02d}:{:02d}\n'.format(hours // 3600, minutes // 60, seconds))
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


def unique_exists(pattern, allow_none=False, get_newest=False):
    result = bash_expands(pattern)
    if len(result) == 1:
        return result[0]
    elif result:
        if get_newest:
            return max(result, key=lambda x: os.stat(x).st_mtime)
        raise OSError('duplicate results:', result)
    elif allow_none:
        return None
    else:
        raise OSError('no result found', pattern)


def read_bpp_from_input_file(filename, expand_orient=False, expand_strand=False, expand_svtype=False, **kwargs):
    """
    reads a file using the tab module. Each row is converted to a breakpoint pair and
    other column data is stored in the data attribute

    Args:
        filename (str): path to the input file
        expand_ns (bool): expand not specified orient/strand settings to all specific version
            (for strand this is only applied if the bam itself is stranded)
        explicit_strand (bool): used to stop unstranded breakpoint pairs from losing input strand information
    Returns:
        :class:`list` of :any:`BreakpointPair`: a list of pairs

    Example:
        >>> read_bpp_from_input_file('filename')
        [BreakpointPair(), BreakpointPair(), ...]

    One can also validate other expected columns that will go in the data attribute using the usual arguments
    to the tab.read_file function

    .. code-block:: python

        >>> read_bpp_from_input_file('filename', cast={'index': int})
        [BreakpointPair(), BreakpointPair(), ...]
    """
    def soft_null_cast(value):
        try:
            tab.cast_null(value)
        except TypeError:
            return value
    kwargs['require'] = set() if 'require' not in kwargs else set(kwargs['require'])
    kwargs['require'].update({COLUMNS.break1_chromosome, COLUMNS.break2_chromosome})
    kwargs.setdefault('cast', {}).update(
        {
            COLUMNS.break1_position_start: int,
            COLUMNS.break1_position_end: int,
            COLUMNS.break2_position_start: int,
            COLUMNS.break2_position_end: int,
            COLUMNS.opposing_strands: lambda x: None if x == '?' else soft_cast(x, cast_type=bool),
            COLUMNS.stranded: tab.cast_boolean,
            COLUMNS.untemplated_seq: soft_null_cast,
            COLUMNS.break1_chromosome: lambda x: re.sub('^chr', '', x),
            COLUMNS.break2_chromosome: lambda x: re.sub('^chr', '', x),
            COLUMNS.tracking_id: lambda x: x if x else str(uuid())
        })
    kwargs.setdefault('add_default', {}).update({
        COLUMNS.untemplated_seq: None,
        COLUMNS.break1_orientation: ORIENT.NS,
        COLUMNS.break1_strand: STRAND.NS,
        COLUMNS.break2_orientation: ORIENT.NS,
        COLUMNS.break2_strand: STRAND.NS,
        COLUMNS.opposing_strands: None,
        COLUMNS.tracking_id: ''
    })
    kwargs.setdefault('in_', {}).update(
        {
            COLUMNS.break1_orientation: ORIENT.values(),
            COLUMNS.break1_strand: STRAND.values(),
            COLUMNS.break2_orientation: ORIENT.values(),
            COLUMNS.break2_strand: STRAND.values()
        })
    _, rows = tab.read_file(
        filename, suppress_index=True,
        **kwargs
    )
    restricted = [
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
        COLUMNS.untemplated_seq
    ]
    pairs = []
    for line_index, row in enumerate(rows):
        row['line_no'] = line_index + 1
        if '_index' in row:
            del row['_index']
        for attr, val in row.items():
            row[attr] = soft_null_cast(val)
        for attr in row:
            if attr in [COLUMNS.cluster_id, COLUMNS.annotation_id, COLUMNS.validation_id]:
                if not re.match('^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$', row[attr]):
                    raise AssertionError(
                        'error in column', attr, 'All mavis pipeline step ids must satisfy the regex:',
                        '^([A-Za-z0-9-]+|)(;[A-Za-z0-9-]+)*$', row[attr])
        stranded = row[COLUMNS.stranded]

        strand1 = row[COLUMNS.break1_strand] if stranded else STRAND.NS
        strand2 = row[COLUMNS.break2_strand] if stranded else STRAND.NS

        temp = []
        expand_strand = stranded and expand_strand
        event_type = [None]
        if row.get(COLUMNS.event_type, None) not in [None, 'None']:
            try:
                event_type = row[COLUMNS.event_type].split(';')
                for putative_event_type in event_type:
                    SVTYPE.enforce(putative_event_type)
            except KeyError:
                pass

        for orient1, orient2, strand1, strand2, putative_event_type in itertools.product(
            ORIENT.expand(row[COLUMNS.break1_orientation]) if expand_orient else [row[COLUMNS.break1_orientation]],
            ORIENT.expand(row[COLUMNS.break2_orientation]) if expand_orient else [row[COLUMNS.break2_orientation]],
            STRAND.expand(strand1) if expand_strand and stranded else [strand1],
            STRAND.expand(strand2) if expand_strand and stranded else [strand2],
            event_type
        ):
            try:
                break1 = Breakpoint(
                    row[COLUMNS.break1_chromosome],
                    row[COLUMNS.break1_position_start],
                    row[COLUMNS.break1_position_end],
                    strand=strand1,
                    orient=orient1
                )
                break2 = Breakpoint(
                    row[COLUMNS.break2_chromosome],
                    row[COLUMNS.break2_position_start],
                    row[COLUMNS.break2_position_end],
                    strand=strand2,
                    orient=orient2
                )

                data = {k: v for k, v in row.items() if k not in restricted}
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
                            'error: expected one of', BreakpointPair.classify(bpp),
                            'but found', putative_event_type, str(bpp), row)
                if expand_svtype and not putative_event_type:
                    for svtype in BreakpointPair.classify(bpp, distance=lambda x, y: Interval(y - x)):
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
