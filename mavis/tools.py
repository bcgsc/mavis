import glob
import itertools
import re
import warnings

from braceexpand import braceexpand
from shortuuid import uuid
import tab
from pysam import VariantFile
from argparse import Namespace

from .breakpoint import Breakpoint, BreakpointPair
from .constants import COLUMNS, MavisNamespace, ORIENT, STRAND, SVTYPE
from .error import InvalidRearrangement
from .util import devnull, read_bpp_from_input_file

SUPPORTED_TOOL = MavisNamespace(
    MANTA='manta',
    DELLY='delly',
    TA='transabyss',
    PINDEL='pindel',
    CHIMERASCAN='chimerascan',
    MAVIS='mavis',
    DEFUSE='defuse',
    BREAKDANCER='breakdancer',
    VCF='vcf'
)
"""
Supported Tools used to call SVs and then used as input into MAVIS

- chimerascan [Iyer-2011]_
- defuse [McPherson-2011]_
- delly [Rausch-2012]_
- manta [Chen-2016]_
- pindel [Ye-2009]_
- transabyss [Robertson-2010]_
"""

TOOL_SVTYPE_MAPPING = {v: [v] for v in SVTYPE.values()}
TOOL_SVTYPE_MAPPING.update({
    'DEL': [SVTYPE.DEL],
    'INS': [SVTYPE.INS],
    'ITX': [SVTYPE.DUP],
    'CTX': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'INV': [SVTYPE.INV],
    'BND': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'TRA': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'CNV': [SVTYPE.DUP],
    'RPL': [SVTYPE.INS],
    'DUP:TANDEM': [SVTYPE.DUP],
    'DUP': [SVTYPE.DUP],
    'interchromosomal': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'eversion': [SVTYPE.DUP],
    'translocation': [SVTYPE.TRANS, SVTYPE.ITRANS],
    'ins': [SVTYPE.INS],
    'del': [SVTYPE.DEL],
    'dup': [SVTYPE.DUP],
    'ITD': [SVTYPE.DUP]
})

TRACKING_COLUMN = 'tracking_id'


def convert_tool_output(fnames, file_type=SUPPORTED_TOOL.MAVIS, stranded=False, log=devnull, collapse=True, assume_no_untemplated=True):
    """
    Reads output from a given SV caller and converts to a set of MAVIS breakpoint pairs. Also collapses duplicates
    """
    result = []
    for fname in fnames:
        result.extend(_convert_tool_output(fname, file_type, stranded, log, assume_no_untemplated=assume_no_untemplated))
    if collapse:
        collapse_mapping = {}
        for bpp in result:
            collapse_mapping.setdefault(bpp, []).append(bpp)
        log('collapsed', len(result), 'to', len(collapse_mapping), 'calls')
        result = list(collapse_mapping.keys())
    return result


def _parse_transabyss(row):
    """
    transforms the transabyss output into the common format for expansion. Maps the input column
    names to column names which MAVIS can read
    """
    std_row = {}
    if TRACKING_COLUMN not in row:
        std_row[TRACKING_COLUMN] = '{}-{}'.format(SUPPORTED_TOOL.TA, row['id'])

    std_row[COLUMNS.event_type] = row.get('rearrangement', row['type'])
    for retained_column in ['genes', 'gene']:
        if retained_column in row:
            std_row['{}_{}'.format(SUPPORTED_TOOL.TA, retained_column)] = row[retained_column]
    if std_row[COLUMNS.event_type] in ['LSR', 'translocation']:
        del std_row[COLUMNS.event_type]
    if 'breakpoint' in row:
        std_row[COLUMNS.break1_orientation], std_row[COLUMNS.break2_orientation] = row['orientations'].split(',')
        match = re.match(
            r'^(?P<chr1>[^:]+):(?P<pos1_start>\d+)\|(?P<chr2>[^:]+):(?P<pos2_start>\d+)$', row['breakpoint'])
        if not match:
            raise OSError(
                'file format error: the breakpoint column did not satisfy the expected pattern', row)
        for group, col in zip(
            ['chr1', 'pos1_start', 'chr2', 'pos2_start'],
            [COLUMNS.break1_chromosome, COLUMNS.break1_position_start, COLUMNS.break2_chromosome, COLUMNS.break2_position_start]
        ):
            std_row[col] = match[group]
    else:
        std_row.update({
            COLUMNS.break1_chromosome: row['chr'],
            COLUMNS.break1_position_start: int(row['chr_start']),
            COLUMNS.break2_position_start: int(row['chr_end'])
        })
        if std_row[COLUMNS.event_type] == 'del':
            std_row[COLUMNS.break1_position_start] -= 1
            std_row[COLUMNS.break2_position_start] += 1
        elif std_row[COLUMNS.event_type] == 'ins':
            std_row[COLUMNS.break2_position_start] += 1

        # add the untemplated sequence where appropriate
        if std_row[COLUMNS.event_type] == 'del':
            assert row['alt'] == 'na'
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row[COLUMNS.event_type] in ['dup', 'ITD']:
            length = std_row[COLUMNS.break2_position_start] - std_row[COLUMNS.break1_position_start] + 1
            if len(row['alt']) != length:
                raise AssertionError(
                    'expected alternate sequence to be equal to the length of the event',
                    len(row['alt']), length, row, std_row)
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row[COLUMNS.event_type] == 'ins':
            std_row[COLUMNS.untemplated_seq] = row['alt'].upper()
        else:
            raise NotImplementedError('unexpected indel type', std_row[COLUMNS.event_type])
    return std_row


def _parse_chimerascan(row):
    """
    transforms the chimerscan output into the common format for expansion. Maps the input column
    names to column names which MAVIS can read
    """
    std_row = {}
    for retained_column in ['genes5p', 'genes3p']:
        if retained_column in row:
            std_row['{}_{}'.format(SUPPORTED_TOOL.CHIMERASCAN, retained_column)] = row[retained_column]
    if TRACKING_COLUMN not in row:
        std_row[TRACKING_COLUMN] = '{}-{}'.format(SUPPORTED_TOOL.CHIMERASCAN, row['chimera_cluster_id'])

    std_row.update({COLUMNS.break1_chromosome: row['chrom5p'], COLUMNS.break2_chromosome: row['chrom3p']})
    if row['strand5p'] == '+':
        std_row[COLUMNS.break1_position_start] = row['end5p']
        std_row[COLUMNS.break1_orientation] = ORIENT.LEFT
    else:
        std_row[COLUMNS.break1_position_start] = row['start5p']
        std_row[COLUMNS.break1_orientation] = ORIENT.RIGHT
    if row['strand3p'] == '+':
        std_row[COLUMNS.break2_position_start] = row['start3p']
        std_row[COLUMNS.break2_orientation] = ORIENT.RIGHT
    else:
        std_row[COLUMNS.break2_position_start] = row['end3p']
        std_row[COLUMNS.break2_orientation] = ORIENT.LEFT
    std_row[COLUMNS.opposing_strands] = row['strand5p'] != row['strand3p']
    return std_row


def _parse_bnd_alt(alt):
    """
    parses the alt statement from vcf files using the specification in vcf 4.2/4.2.

    Assumes that the reference base is always the outermost base (this is based on the spec and also manta results as
    the spec was missing some cases)
    """
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)\[(?P<chr>[^:]+):(?P<pos>\d+)\[$', alt)
    if match:
        return (match.group('chr'), int(match.group('pos')), ORIENT.RIGHT, match.group('ref'), match.group('useq'))
    match = re.match(r'^\[(?P<chr>[^:]+):(?P<pos>\d+)\[(?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (match.group('chr'), int(match.group('pos')), ORIENT.RIGHT, match.group('ref'), match.group('useq'))
    match = re.match(r'^\](?P<chr>[^:]+):(?P<pos>\d+)\](?P<useq>\w*)(?P<ref>\w)$', alt)
    if match:
        return (match.group('chr'), int(match.group('pos')), ORIENT.LEFT, match.group('ref'), match.group('useq'))
    match = re.match(r'^(?P<ref>\w)(?P<useq>\w*)](?P<chr>[^:]+):(?P<pos>\d+)\]$', alt)
    if match:
        return (match.group('chr'), int(match.group('pos')), ORIENT.LEFT, match.group('ref'), match.group('useq'))
    else:
        raise NotImplementedError('alt specification in unexpected format', alt)


def _parse_vcf_record(record):
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
        for entry in record.info.items():
            info[entry[0]] = entry[1:] if len(entry[1:]) > 1 else entry[1]
        std_row = {}
        if record.id and record.id != 'N':  # to account for NovoBreak N in the ID field
            std_row['id'] = record.id

        if info.get('SVTYPE', None) == 'BND':
            chr2, end, orient2, ref, alt = _parse_bnd_alt(alt)
            std_row[COLUMNS.break2_orientation] = orient2
            std_row[COLUMNS.untemplated_seq] = alt
            if record.ref != ref:
                raise AssertionError(
                    'Expected the ref specification in the vcf record to match the sequence '
                    'in the alt string: {} vs {}'.format(record.ref, ref))
        else:
            chr2 = info.get('CHR2', record.chrom)
            end = record.stop
            if alt and record.ref and re.match(r'^[A-Z]+$', alt) and re.match(r'^[A-Z]+', record.ref):
                std_row[COLUMNS.untemplated_seq] = alt[1:]
                size = len(alt) - len(record.ref)
                if size > 0:
                    std_row[COLUMNS.event_type] = SVTYPE.INS
                elif size < 0:
                    std_row[COLUMNS.event_type] = SVTYPE.DEL
        std_row.update({COLUMNS.break1_chromosome: record.chrom, COLUMNS.break2_chromosome: chr2})
        if info.get('PRECISE', False):  # DELLY CI only apply when split reads were not used to refine the breakpoint which is then flagged
            std_row.update({
                COLUMNS.break1_position_start: record.pos,
                COLUMNS.break1_position_end: record.pos,
                COLUMNS.break2_position_start: end,
                COLUMNS.break2_position_end: end
            })
        else:
            std_row.update({
                COLUMNS.break1_position_start: max(1, record.pos + info.get('CIPOS', (0, 0))[0]),
                COLUMNS.break1_position_end: record.pos + info.get('CIPOS', (0, 0))[1],
                COLUMNS.break2_position_start: max(1, end + info.get('CIEND', (0, 0))[0]),
                COLUMNS.break2_position_end: end + info.get('CIEND', (0, 0))[1]
            })

        if 'SVTYPE' in info:
            std_row[COLUMNS.event_type] = info['SVTYPE']

        try:
            orient1, orient2 = info['CT'].split('to')
            connection_type = {'3': ORIENT.LEFT, '5': ORIENT.RIGHT, 'N': ORIENT.NS}
            std_row[COLUMNS.break1_orientation] = connection_type[orient1]
            std_row[COLUMNS.break2_orientation] = connection_type[orient2]
        except KeyError:
            pass
        records.append(std_row)
    return records


def _convert_tool_row(row, file_type, stranded, assume_no_untemplated=True):
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
    if file_type in [SUPPORTED_TOOL.DELLY, SUPPORTED_TOOL.MANTA, SUPPORTED_TOOL.PINDEL, SUPPORTED_TOOL.VCF]:

        std_row.update(row)

    elif file_type == SUPPORTED_TOOL.CHIMERASCAN:

        std_row.update(_parse_chimerascan(row))

    elif file_type == SUPPORTED_TOOL.DEFUSE:

        std_row[COLUMNS.break1_orientation] = ORIENT.LEFT if row['genomic_strand1'] == STRAND.POS else ORIENT.RIGHT
        std_row[COLUMNS.break2_orientation] = ORIENT.LEFT if row['genomic_strand2'] == STRAND.POS else ORIENT.RIGHT
        std_row.update({
            COLUMNS.break1_chromosome: row['gene_chromosome1'],
            COLUMNS.break2_chromosome: row['gene_chromosome2'],
            COLUMNS.break1_position_start: row['genomic_break_pos1'],
            COLUMNS.break2_position_start: row['genomic_break_pos2']
        })
        if TRACKING_COLUMN in row:
            std_row[TRACKING_COLUMN] = row[TRACKING_COLUMN]
        else:
            std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, row['cluster_id'])

    elif file_type == SUPPORTED_TOOL.TA:

        std_row.update(_parse_transabyss(row))

    elif file_type == SUPPORTED_TOOL.BREAKDANCER:

        std_row.update({
            COLUMNS.event_type: row['Type'],
            COLUMNS.break1_chromosome: row['Chr1'],
            COLUMNS.break2_chromosome: row['Chr2'],
            COLUMNS.break1_position_start: row['Pos1'],
            COLUMNS.break2_position_start: row['Pos2'],
        })

    else:
        raise NotImplementedError('unsupported file type', file_type)

    if stranded:
        std_row[COLUMNS.break1_strand] = STRAND.expand(std_row[COLUMNS.break1_strand])
        std_row[COLUMNS.break2_strand] = STRAND.expand(std_row[COLUMNS.break2_strand])
    else:
        std_row[COLUMNS.break1_strand] = [STRAND.NS]
        std_row[COLUMNS.break2_strand] = [STRAND.NS]

    if not std_row.get(TRACKING_COLUMN, None):
        std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, std_row.get('id', uuid()))

    combinations = list(itertools.product(
        ORIENT.expand(std_row[COLUMNS.break1_orientation]), ORIENT.expand(std_row[COLUMNS.break2_orientation]),
        std_row[COLUMNS.break1_strand], std_row[COLUMNS.break2_strand],
        TOOL_SVTYPE_MAPPING[std_row[COLUMNS.event_type]] if COLUMNS.event_type in std_row else [None],
        [True, False] if std_row.get(COLUMNS.opposing_strands, None) is None else [std_row[COLUMNS.opposing_strands]]
    ))
    # add the product of all uncertainties as breakpoint pairs
    for orient1, orient2, strand1, strand2, event_type, oppose in combinations:
        try:
            untemplated_seq = std_row.get(COLUMNS.untemplated_seq, None)
            if assume_no_untemplated and event_type != SVTYPE.INS and not untemplated_seq:
                untemplated_seq = ''
            break1 = Breakpoint(
                std_row[COLUMNS.break1_chromosome],
                std_row[COLUMNS.break1_position_start],
                std_row.get(COLUMNS.break1_position_end, std_row[COLUMNS.break1_position_start]),
                orient=orient1, strand=strand1
            )
            break2 = Breakpoint(
                std_row.get(COLUMNS.break2_chromosome, std_row[COLUMNS.break1_chromosome]),
                std_row[COLUMNS.break2_position_start],
                std_row.get(COLUMNS.break2_position_end, std_row[COLUMNS.break2_position_start]),
                orient=orient2, strand=strand2
            )
            if len(break1) == 1 and len(break2) == 1 and event_type == SVTYPE.DEL and abs(break1.start - break2.start) < 2:
                break1 = Breakpoint(break1.chr, break1.start - 1, break1.end - 1, orient=break1.orient, strand=break1.strand)
                break2 = Breakpoint(break2.chr, break2.start + 1, break2.end + 1, orient=break2.orient, strand=break2.strand)
            bpp = BreakpointPair(
                break1, break2,
                opposing_strands=oppose,
                untemplated_seq=untemplated_seq,
                event_type=event_type,
                data={
                    COLUMNS.tools: file_type,
                    COLUMNS.tracking_id: std_row[COLUMNS.tracking_id]
                },
                stranded=stranded
            )

            bpp.data.update({k: std_row[k] for k in std_row if k.startswith(file_type)})
            if not event_type or event_type in BreakpointPair.classify(bpp):
                result.append(bpp)

        except (InvalidRearrangement, AssertionError):
            pass
    if not result:
        raise UserWarning(
            'row failed to create any breakpoint pairs. This generally indicates an input formatting error',
            row, std_row, combinations)
    return result


def _convert_tool_output(input_file, file_type=SUPPORTED_TOOL.MAVIS, stranded=False, log=devnull, assume_no_untemplated=True):
    log('reading:', input_file)
    result = []
    rows = None
    if file_type == SUPPORTED_TOOL.MAVIS:
        result = read_bpp_from_input_file(input_file, expand_orient=True, expand_svtype=True, add_default={'stranded': stranded})
    elif file_type in [SUPPORTED_TOOL.DELLY, SUPPORTED_TOOL.MANTA, SUPPORTED_TOOL.PINDEL]:
        rows = []
        for vcf_record in VariantFile(input_file).fetch():
            rows.extend(_parse_vcf_record(vcf_record))
    elif file_type == SUPPORTED_TOOL.BREAKDANCER:
        with open(input_file, 'r') as fh:
            # comments in breakdancer are marked with a single # so they need to be discarded before reading
            lines = fh.readlines()
            header = 0
            while header < len(lines) and lines[header].startswith('#'):
                header += 1
            lines = lines[header - 1:]
            input_file = Namespace(readlines=lambda: lines)
        _, rows = tab.read_file(input_file, allow_short=True)
    else:
        _, rows = tab.read_file(input_file)
    if rows:
        log('found', len(rows), 'rows')
        for row in rows:
            std_rows = _convert_tool_row(row, file_type, stranded, assume_no_untemplated=assume_no_untemplated)
            result.extend(std_rows)
    log('generated', len(result), 'breakpoint pairs')
    return result
