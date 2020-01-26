"""
Convert between VCF and MAVIS outputs
"""
import re

from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import COLUMNS, ORIENT, SVTYPE
from ..util import DEVNULL


def parse_bnd_alt(alt):
    """
    parses the alt statement from vcf files using the specification in vcf 4.2/4.2.

    Assumes that the reference base is always the outermost base (this is based on the spec and also manta results as
    the spec was missing some cases)
    """
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)\[(?P<chr>[^:]+):(?P<pos>\d+)\[$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.RIGHT,
            match.group('ref'),
            match.group('useq'),
        )
    match = re.match(r'^\[(?P<chr>[^:]+):(?P<pos>\d+)\[(?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.RIGHT,
            match.group('ref'),
            match.group('useq'),
        )
    match = re.match(r'^\](?P<chr>[^:]+):(?P<pos>\d+)\](?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.LEFT,
            match.group('ref'),
            match.group('useq'),
        )
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)](?P<chr>[^:]+):(?P<pos>\d+)\]$', alt)
    if match:
        return (
            match.group('chr'),
            int(match.group('pos')),
            ORIENT.LEFT,
            match.group('ref'),
            match.group('useq'),
        )
    else:
        raise NotImplementedError('alt specification in unexpected format', alt)


def generate_bnd_alt(bpp):
    """convert a mavis breakpoint pair to the VCF BND Alt string specification equivalents

    Args:
        bpp (Breakpoint): the breakpoint to be converted

    Raises:
        NotImplementedError: [description]

    Cases (VCF 4.2)
    REF ALT Meaning
    1. s t[p[ piece extending to the right of p is joined after t
    2. s t]p] reverse comp piece extending left of p is joined after t
    3. s ]p]t piece extending to the left of p is joined before t
    4. s [p[t reverse comp piece extending right of p is joined before t

    These translate to MAVIS current and mate breakpoint orientations of
    1. LR
    2. LL
    3. RL
    4. RR
    """
    if bpp.LR:
        return (
            f'N{bpp.untemplated_seq}[{bpp.break2.chr}:{bpp.break2.pos}[',
            f']{bpp.break1.chr}:{bpp.break1.pos}]{bpp.untemplated_seq}N',
        )
    elif bpp.LL:
        pass
    elif bpp.RL:
        pass
    elif bpp.RR:
        pass
    else:
        raise NotImplementedError(
            'Cannot convert breakpoint to BND alt notation with unknown orientation'
        )


def parse_record(record, log=DEVNULL):
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
                log('Ignoring invalid INFO field {} with error: {}'.format(key, err))
            else:
                try:
                    value = value[0] if len(value) == 1 else value
                except TypeError:
                    pass  # anything non-tuple
            info[key] = value

        std_row = {}
        if record.id and record.id != 'N':  # to account for NovoBreak N in the ID field
            std_row['id'] = record.id

        if info.get('SVTYPE', None) == 'BND':
            chr2, end, orient2, ref, alt = parse_bnd_alt(alt)
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
