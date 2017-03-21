from datetime import datetime
import errno
import math
import os
import random
from mavis.breakpoint import read_bpp_from_input_file
from mavis.constants import PROTOCOL, COLUMNS, sort_columns
from mavis.interval import Interval
from vocab import Vocab


PIPELINE_STEP = Vocab(
    ANNOTATE='annotate',
    VALIDATE='validate',
    PIPELINE='pipeline',
    CLUSTER='cluster',
    PAIR='pairing',
    SUMMARY='summary'
)


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
        log('loading:', finput)
        bpps.extend(read_bpp_from_input_file(
            finput, force_stranded=force_stranded,
            **kwargs
        ))
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
