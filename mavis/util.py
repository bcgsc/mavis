from datetime import datetime
import errno
import math
import os
import random
from .breakpoint import read_bpp_from_input_file
from .constants import PROTOCOL, COLUMNS, sort_columns
from .interval import Interval
from argparse import Namespace
from TSV.TSV import EmptyFileError


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
    for finput in inputs:
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


def generate_complete_stamp(output_dir, log=devnull, prefix='MAVIS.'):
    stamp = os.path.join(output_dir, str(prefix) + 'COMPLETE')
    log('complete:', stamp)
    with open(stamp, 'w') as fh:
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
