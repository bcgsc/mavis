from argparse import Namespace
from datetime import datetime
import errno
from glob import glob
import os
import re

from braceexpand import braceexpand
from TSV.TSV import EmptyFileError, tsv_boolean

from .breakpoint import read_bpp_from_input_file
from .constants import COLUMNS, PROTOCOL, sort_columns
from .interval import Interval

ENV_VAR_PREFIX = 'MAVIS_'


def cast(value, cast_func):
    if cast_func == bool:
        value = tsv_boolean(value)
    else:
        value = cast_func(value)
    return value


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
    else:
        return default


class MavisNamespace(Namespace):

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __add__(self, other):
        d = {}
        d.update(self.__dict__)
        d.update(other.__dict__)
        return MavisNamespace(**d)

    def update(self, other):
        self.__dict__.update(other.__dict__)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, val):
        self.__dict__[key] = val

    def flatten(self):
        d = {}
        d.update(self.items())
        return d

    def get(self, key, default):
        try:
            return self[key]
        except AttributeError:
            return default

    def keys(self):
        return self.__dict__.keys()


class WeakMavisNamespace(MavisNamespace):

    def __getattribute__(self, attr):
        return get_env_variable(attr, object.__getattribute__(self, attr))


class ChrListString(list):

    def __init__(self, string):
        if not isinstance(string, str):
            for item in string:
                self.append(item)
        else:
            delim = '\s+' if ';' not in string else ';'
            items = [i for i in re.split(delim, string) if i]
            for item in items:
                self.append(item)

    def __contains__(self, item):
        if list.__len__(self) == 0:
            return True
        else:
            return list.__contains__(self, item)


def bash_expands(expression):
    result = []
    for name in braceexpand(expression):
        for fname in glob(name):
            result.append(fname)
    return result


def log_arguments(args):
    log('arguments')
    for arg, val in sorted(args.items()):
        if isinstance(val, list):
            log(arg, '= [', time_stamp=False)
            for v in val:
                log('\t', repr(v), time_stamp=False)
            log(']', time_stamp=False)
        elif any([isinstance(val, typ) for typ in [str, int, float, bool, tuple]]) or val is None:
            log(arg, '=', repr(val), time_stamp=False)
        else:
            log(arg, '=', object.__repr__(val), time_stamp=False)


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def devnull(*pos, **kwargs):
    pass


def mkdirp(dirname):
    log("creating output directory: '{}'".format(dirname))
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise exc
    return dirname


def filter_on_overlap(bpps, regions_by_reference_name):
    log('filtering from', len(bpps), 'using overlaps with regions filter')
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
    log('filtered from', len(bpps), 'down to', len(passed), '(removed {})'.format(len(failed)))
    return passed, failed


def read_inputs(inputs, **kwargs):
    bpps = []
    kwargs.setdefault('require', [])
    kwargs['require'] = list(set(kwargs['require'] + [COLUMNS.protocol]))
    kwargs.setdefault('in_', {})
    kwargs['in_'][COLUMNS.protocol] = PROTOCOL
    for expr in inputs:
        for finput in bash_expands(expr):
            try:
                log('loading:', finput)
                bpps.extend(read_bpp_from_input_file(
                    finput,
                    **kwargs
                ))
            except EmptyFileError:
                log('ignoring empty file:', finput)
    log('loaded', len(bpps), 'breakpoint pairs')
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
        log('writing:', filename)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row.get(c, None)) for c in header]) + '\n')


def write_bed_file(filename, bed_rows):
    log('writing:', filename)
    with open(filename, 'w') as fh:
        for bed in bed_rows:
            fh.write('\t'.join([str(c) for c in bed]) + '\n')


def generate_complete_stamp(output_dir, log=devnull, prefix='MAVIS.'):
    stamp = os.path.join(output_dir, str(prefix) + 'COMPLETE')
    log('complete:', stamp)
    with open(stamp, 'w'):
        pass
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
