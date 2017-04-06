from datetime import datetime
import errno
import math
import os
import random
from .breakpoint import read_bpp_from_input_file
from .constants import PROTOCOL, COLUMNS, sort_columns
from .interval import Interval
from argparse import Namespace
import subprocess
from TSV.TSV import EmptyHeaderError


class MavisNamespace(Namespace):
    def items(self):
        return self.__dict__.items()

    def __add__(self, other):
        d = {}
        d.update(self.__dict__)
        d.update(other.__dict__)
        return MavisNamespace(**d)

    def update(self, other):
        self.__dict__.update(other.__dict__)

    def __getitem__(self, key):
        return getattr(self, key)


def get_version():
    v = subprocess.check_output('cd {}; git describe'.format(os.path.dirname(__file__)), shell=True)
    v = v.decode('UTF8')
    v = v.strip()
    return v


def build_batch_id(prefix='', suffix='', size=6):
    date = datetime.now()
    m = int(math.pow(10, size) - 1)
    return '{prefix}batch{date.year}{date.month:02d}{date.day:02d}r{r:06d}{suffix}'.format(
        prefix=prefix, suffix=suffix, date=date, r=random.randint(1, m))


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
    log('filtering', len(bpps), 'on overlaps with regions')
    failed = []
    passed = []
    for bpp in bpps:
        overlaps = False
        for r in regions_by_reference_name.get(bpp.break1.chr, []):
            if Interval.overlaps(r, bpp.break1):
                overlaps = True
                bpp.data['failure_comment'] = 'overlapped masked region: ' + str(r)
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
    log('filtered', len(bpps), 'to', len(passed))
    return passed, failed


def read_inputs(inputs, force_stranded=False, **kwargs):
    bpps = []
    kwargs.setdefault('require', [])
    kwargs['require'] = list(set(kwargs['require'] + [COLUMNS.library, COLUMNS.protocol]))
    kwargs.setdefault('in_', {})
    kwargs['in_'][COLUMNS.protocol] = PROTOCOL
    for finput in inputs:
        try:
            log('loading:', finput)
            bpps.extend(read_bpp_from_input_file(
                finput, force_stranded=force_stranded,
                **kwargs
            ))
        except EmptyHeaderError:
            log('ignoring empty file:', finput)
    log('loaded', len(bpps), 'breakpoint pairs')
    return bpps


def output_tabbed_file(bpps, filename):
    header = set()
    rows = []
    for row in bpps:
        try:
            row = row.flatten()
        except AttributeError:
            pass
        rows.append(row)
        header.update(row)

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
