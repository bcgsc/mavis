import glob
import itertools
import re
import warnings

from braceexpand import braceexpand
from shortuuid import uuid
import tab
from pysam import VariantFile

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
    DEFUSE='defuse'
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


def convert_tool_output(input_file, file_type=SUPPORTED_TOOL.MAVIS, stranded=False, log=devnull, collapse=True, assume_no_untemplated=True):
    """
    Reads output from a given SV caller and converts to a set of MAVIS breakpoint pairs. Also collapses duplicates
    """
    result = []
    fnames = []
    for name in braceexpand(input_file):
        for subname in glob.glob(name):
            fnames.append(subname)
    if not fnames:
        raise OSError('no such file', input_file)
    for fname in fnames:
        result.extend(_convert_tool_output(fname, file_type, stranded, log, assume_no_untemplated=assume_no_untemplated))
    if collapse:
        collapse_mapping = {}
        for bpp in result:
            collapse_mapping.setdefault(bpp, []).append(bpp)
        log('collapsed', len(result), 'to', len(collapse_mapping), 'calls')
        result = list(collapse_mapping.keys())
    return result


def _parse_breakdancer(row):
    pass  # TODO: add breakdancer support


def _parse_transabyss(row):
    """
    transforms the transabyss output into the common format for expansion. Maps the input column
    names to column names which MAVIS can read
    """
    std_row = {}
    if TRACKING_COLUMN not in row:
        std_row[TRACKING_COLUMN] = '{}-{}'.format(SUPPORTED_TOOL.TA, row['id'])

    std_row['event_type'] = row.get('rearrangement', row['type'])
    for retained_column in ['genes', 'gene']:
        if retained_column in row:
            std_row['{}_{}'.format(SUPPORTED_TOOL.TA, retained_column)] = row[retained_column]
    if std_row['event_type'] in ['LSR', 'translocation']:
        del std_row['event_type']
    if 'breakpoint' in row:
        std_row['orient1'], std_row['orient2'] = row['orientations'].split(',')
        match = re.match(
            r'^(?P<chr1>[^:]+):(?P<pos1_start>\d+)\|(?P<chr2>[^:]+):(?P<pos2_start>\d+)$', row['breakpoint'])
        if not match:
            raise OSError(
                'file format error: the breakpoint column did not satisfy the expected pattern', row)
        std_row.update({k: match[k] for k in ['chr1', 'pos1_start', 'chr2', 'pos2_start']})
    else:
        std_row.update({
            'chr1': row['chr'], 'pos1_start': int(row['chr_start']), 'pos2_start': int(row['chr_end'])
        })
        if std_row['event_type'] == 'del':
            std_row['pos1_start'] -= 1
            std_row['pos2_start'] += 1
        elif std_row['event_type'] == 'ins':
            std_row['pos2_start'] += 1

        # add the untemplated sequence where appropriate
        if std_row['event_type'] == 'del':
            assert row['alt'] == 'na'
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row['event_type'] in ['dup', 'ITD']:
            length = std_row['pos2_start'] - std_row['pos1_start'] + 1
            if len(row['alt']) != length:
                raise AssertionError(
                    'expected alternate sequence to be equal to the length of the event',
                    len(row['alt']), length, row, std_row)
            std_row[COLUMNS.untemplated_seq] = ''
        elif std_row['event_type'] == 'ins':
            std_row[COLUMNS.untemplated_seq] = row['alt'].upper()
        else:
            raise NotImplementedError('unexpected indel type', std_row['event_type'])
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

    std_row.update({'chr1': row['chrom5p'], 'chr2': row['chrom3p']})
    if row['strand5p'] == '+':
        std_row['pos1_start'] = row['end5p']
        std_row['orient1'] = ORIENT.LEFT
    else:
        std_row['pos1_start'] = row['start5p']
        std_row['orient1'] = ORIENT.RIGHT
    if row['strand3p'] == '+':
        std_row['pos2_start'] = row['start3p']
        std_row['orient2'] = ORIENT.RIGHT
    else:
        std_row['pos2_start'] = row['end3p']
        std_row['orient2'] = ORIENT.LEFT
    std_row['opposing_strands'] = row['strand5p'] != row['strand3p']
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


def _parse_vcf_record(row):
    """
    converts a vcf record
    """
    records = []
    for alt in row.alts if row.alts else [None]:
        info = {}
        for entry in row.info.items():
            info[entry[0]] = entry[1:] if len(entry[1:]) > 1 else entry[1]
        std_row = {}
        if row.id:
            std_row['id'] = row.id

        if info.get('SVTYPE', None) == 'BND':
            chr2, end, orient2, ref, alt = _parse_bnd_alt(alt)
            std_row['orient2'] = orient2
            std_row[COLUMNS.untemplated_seq] = alt
            if row.ref != ref:
                raise AssertionError(
                    'Expected the ref specification in the vcf row to match the sequence '
                    'in the alt string: {} vs {}'.format(row.ref, ref))
        else:
            chr2 = info.get('CHR2', row.chrom)
            end = row.stop
            if alt and row.ref and re.match(r'^[A-Z]+$', alt) and re.match(r'^[A-Z]+', row.ref):
                std_row[COLUMNS.untemplated_seq] = alt[1:]
                size = len(alt) - len(row.ref)
                if size > 0:
                    std_row['event_type'] = SVTYPE.INS
                elif size < 0:
                    std_row['event_type'] = SVTYPE.DEL

        std_row.update({
            'chr1': row.chrom, 'chr2': chr2,
            'pos1_start': max(1, row.pos + info.get('CIPOS', (0, 0))[0]),
            'pos1_end': row.pos + info.get('CIPOS', (0, 0))[1],
            'pos2_start': max(1, end + info.get('CIEND', (0, 0))[0]),
            'pos2_end': end + info.get('CIEND', (0, 0))[1]
        })
        if 'SVTYPE' in info:
            std_row['event_type'] = info['SVTYPE']

        try:
            orient1, orient2 = info['CT'].split('to')
            connection_type = {'3': ORIENT.LEFT, '5': ORIENT.RIGHT, 'N': ORIENT.NS}
            std_row['orient1'] = connection_type[orient1]
            std_row['orient2'] = connection_type[orient2]
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
    std_row['orient1'] = std_row['orient2'] = ORIENT.NS
    std_row['strand1'] = std_row['strand2'] = STRAND.NS
    result = []
    # convert the specified file type to a standard format
    if file_type in [SUPPORTED_TOOL.DELLY, SUPPORTED_TOOL.MANTA, SUPPORTED_TOOL.PINDEL]:

        std_row.update(row)

    elif file_type == SUPPORTED_TOOL.CHIMERASCAN:

        std_row.update(_parse_chimerascan(row))

    elif file_type == SUPPORTED_TOOL.DEFUSE:

        std_row['orient1'] = ORIENT.LEFT if row['genomic_strand1'] == STRAND.POS else ORIENT.RIGHT
        std_row['orient2'] = ORIENT.LEFT if row['genomic_strand2'] == STRAND.POS else ORIENT.RIGHT
        std_row.update({
            'chr1': row['gene_chromosome1'], 'chr2': row['gene_chromosome2'],
            'pos1_start': row['genomic_break_pos1'], 'pos2_start': row['genomic_break_pos2']
        })
        if TRACKING_COLUMN in row:
            std_row[TRACKING_COLUMN] = row[TRACKING_COLUMN]
        else:
            std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, row['cluster_id'])

    elif file_type == SUPPORTED_TOOL.TA:

        std_row.update(_parse_transabyss(row))

    else:
        raise NotImplementedError('unsupported file type', file_type)

    if stranded:
        std_row['strand1'] = STRAND.expand(std_row['strand1'])
        std_row['strand2'] = STRAND.expand(std_row['strand2'])
    else:
        std_row['strand1'] = [STRAND.NS]
        std_row['strand2'] = [STRAND.NS]

    if not std_row.get(TRACKING_COLUMN, None):
        std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, std_row.get('id', uuid()))
    if assume_no_untemplated and not std_row.get(COLUMNS.untemplated_seq, None):
        std_row[COLUMNS.untemplated_seq] = ''

    combinations = list(itertools.product(
        ORIENT.expand(std_row['orient1']), ORIENT.expand(std_row['orient2']),
        std_row['strand1'], std_row['strand2'], TOOL_SVTYPE_MAPPING[std_row['event_type']] if 'event_type' in std_row else [None],
        [True, False] if 'opposing_strands' not in std_row else [std_row['opposing_strands']]
    ))
    # add the product of all uncertainties as breakpoint pairs
    for orient1, orient2, strand1, strand2, event_type, oppose in combinations:
        try:

            bpp = BreakpointPair(
                Breakpoint(
                    std_row['chr1'],
                    std_row['pos1_start'],
                    std_row.get('pos1_end', std_row['pos1_start']),
                    orient=orient1, strand=strand1
                ),
                Breakpoint(
                    std_row.get('chr2', std_row['chr1']),
                    std_row['pos2_start'],
                    std_row.get('pos2_end', std_row['pos2_start']),
                    orient=orient2, strand=strand2
                ),
                opposing_strands=oppose,
                untemplated_seq=std_row.get('untemplated_seq', None),
                event_type=event_type,
                data={
                    COLUMNS.tools: file_type,
                    COLUMNS.tracking_id: std_row['tracking_id']
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
    if file_type == SUPPORTED_TOOL.MAVIS:
        result = read_bpp_from_input_file(input_file)
    else:
        if file_type in [SUPPORTED_TOOL.DELLY, SUPPORTED_TOOL.MANTA, SUPPORTED_TOOL.PINDEL]:
            rows = []
            for vcf_record in VariantFile(input_file).fetch():
                rows.extend(_parse_vcf_record(vcf_record))
        else:
            dummy, rows = tab.read_file(input_file)
        log('found', len(rows), 'rows')
        for row in rows:
            std_rows = _convert_tool_row(row, file_type, stranded, assume_no_untemplated=assume_no_untemplated)
            result.extend(std_rows)
    log('generated', len(result), 'breakpoint pairs')
    return result
