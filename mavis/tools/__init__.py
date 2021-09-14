import itertools

from shortuuid import uuid
import tab

from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import COLUMNS, ORIENT, STRAND, SVTYPE
from ..error import InvalidRearrangement
from ..util import DEVNULL, read_bpp_from_input_file

from .constants import SUPPORTED_TOOL, TRACKING_COLUMN, TOOL_SVTYPE_MAPPING
from .transabyss import convert_row as _parse_transabyss
from .cnvnator import convert_row as _parse_cnvnator
from .vcf import convert_file as read_vcf
from .breakdancer import convert_file as _convert_breakdancer_file
from .starfusion import convert_row as _parse_starfusion
from .chimerascan import convert_row as _parse_chimerascan
from .sniffles import convert_row as _parse_sniffles


def convert_tool_output(
    fnames,
    file_type=SUPPORTED_TOOL.MAVIS,
    stranded=False,
    log=DEVNULL,
    collapse=True,
    assume_no_untemplated=True,
):
    """
    Reads output from a given SV caller and converts to a set of MAVIS breakpoint pairs. Also collapses duplicates
    """
    result = []
    for fname in fnames:
        result.extend(
            _convert_tool_output(
                fname, file_type, stranded, log, assume_no_untemplated=assume_no_untemplated
            )
        )
    if collapse:
        collapse_mapping = {}
        for bpp in result:
            collapse_mapping.setdefault(bpp, []).append(bpp)
        log('collapsed', len(result), 'to', len(collapse_mapping), 'calls')
        result = []
        temp_sets = set()
        for bpp, bpp_list in collapse_mapping.items():
            for otherbpp in bpp_list:
                for col, val in otherbpp.data.items():
                    if val is None:
                        continue
                    if col not in bpp.data or not bpp.data[col]:
                        bpp.data[col] = val
                    elif isinstance(bpp.data[col], set):
                        bpp.data[col].add(val)
                    elif bpp.data[col] != val and bpp.data[col]:
                        bpp.data[col] = {bpp.data[col], val}
                        temp_sets.add(col)
            result.append(bpp)
        for bpp in result:
            for col, val in bpp.data.items():
                if isinstance(val, set) and col in temp_sets:
                    bpp.data[col] = ';'.join(sorted([str(v) for v in val]))
    return result

def _parse_sniffles_vcf(record, log=DEVNULL):
    """
    converts a vcf record from sniffles
    Note:
    Derived from _parse_vcf_record
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
            chr2, end, orient2, ref, alt = _parse_bnd_alt(alt)
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
        # Only include this line if the SVTYPE is recognised
        # This is to remove sniffles weird hybrid SVs eg DEL/INV
        if(info['SVTYPE'] in TOOL_SVTYPE_MAPPING.keys()):
            records.append(std_row)
        else:
            log('Ignoring record with unrecognised SVTYPE {}'.format(info['SVTYPE']))

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
    if file_type in [
        SUPPORTED_TOOL.DELLY,
        SUPPORTED_TOOL.MANTA,
        SUPPORTED_TOOL.PINDEL,
        SUPPORTED_TOOL.VCF,
        SUPPORTED_TOOL.BREAKSEQ,
        SUPPORTED_TOOL.STRELKA,
        SUPPORTED_TOOL.MUTECT,
    ]:

        std_row.update(row)

    elif file_type == SUPPORTED_TOOL.CHIMERASCAN:

        std_row.update(_parse_chimerascan(row))

    elif file_type == SUPPORTED_TOOL.CNVNATOR:

        std_row.update(_parse_cnvnator(row))

    elif file_type == SUPPORTED_TOOL.STARFUSION:

        std_row.update(_parse_starfusion(row))

    elif file_type == SUPPORTED_TOOL.DEFUSE:

        std_row[COLUMNS.break1_orientation] = (
            ORIENT.LEFT if row['genomic_strand1'] == STRAND.POS else ORIENT.RIGHT
        )
        std_row[COLUMNS.break2_orientation] = (
            ORIENT.LEFT if row['genomic_strand2'] == STRAND.POS else ORIENT.RIGHT
        )
        std_row.update(
            {
                COLUMNS.break1_chromosome: row['gene_chromosome1'],
                COLUMNS.break2_chromosome: row['gene_chromosome2'],
                COLUMNS.break1_position_start: row['genomic_break_pos1'],
                COLUMNS.break2_position_start: row['genomic_break_pos2'],
            }
        )
        if TRACKING_COLUMN in row:
            std_row[TRACKING_COLUMN] = row[TRACKING_COLUMN]
        else:
            std_row[TRACKING_COLUMN] = '{}-{}'.format(file_type, row['cluster_id'])

    elif file_type == SUPPORTED_TOOL.TA:

        std_row.update(_parse_transabyss(row))

    elif file_type == SUPPORTED_TOOL.BREAKDANCER:

        std_row.update(
            {
                COLUMNS.event_type: row['Type'],
                COLUMNS.break1_chromosome: row['Chr1'],
                COLUMNS.break2_chromosome: row['Chr2'],
                COLUMNS.break1_position_start: row['Pos1'],
                COLUMNS.break2_position_start: row['Pos2'],
            }
        )
        std_row.update(
            {k: v for k, v in row.items() if k not in {'Type', 'Chr1', 'Chr2', 'Pos1', 'Pos2'}}
        )
    elif file_type == SUPPORTED_TOOL.SNIFFLES:

        print("sniffles work")
        std_row.update(_parse_sniffles(row))

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

    combinations = list(
        itertools.product(
            ORIENT.expand(std_row[COLUMNS.break1_orientation]),
            ORIENT.expand(std_row[COLUMNS.break2_orientation]),
            std_row[COLUMNS.break1_strand],
            std_row[COLUMNS.break2_strand],
            TOOL_SVTYPE_MAPPING[std_row[COLUMNS.event_type]]
            if COLUMNS.event_type in std_row
            else [None],
            [True, False]
            if std_row.get(COLUMNS.opposing_strands, None) is None
            else [std_row[COLUMNS.opposing_strands]],
        )
    )
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
                orient=orient1,
                strand=strand1,
            )
            break2 = Breakpoint(
                std_row.get(COLUMNS.break2_chromosome, std_row[COLUMNS.break1_chromosome]),
                std_row[COLUMNS.break2_position_start],
                std_row.get(COLUMNS.break2_position_end, std_row[COLUMNS.break2_position_start]),
                orient=orient2,
                strand=strand2,
            )
            if (
                len(break1) == 1
                and len(break2) == 1
                and event_type == SVTYPE.DEL
                and abs(break1.start - break2.start) < 2
            ):
                break1 = Breakpoint(
                    break1.chr,
                    break1.start - 1,
                    break1.end - 1,
                    orient=break1.orient,
                    strand=break1.strand,
                )
                break2 = Breakpoint(
                    break2.chr,
                    break2.start + 1,
                    break2.end + 1,
                    orient=break2.orient,
                    strand=break2.strand,
                )
            bpp = BreakpointPair(
                break1,
                break2,
                opposing_strands=oppose,
                untemplated_seq=untemplated_seq,
                event_type=event_type,
                data={COLUMNS.tools: file_type, COLUMNS.tracking_id: std_row[COLUMNS.tracking_id]},
                stranded=stranded,
            )

            for col, value in std_row.items():
                if col not in COLUMNS and col not in bpp.data:
                    bpp.data[col] = value
            if not event_type or event_type in BreakpointPair.classify(bpp):
                result.append(bpp)

        except (InvalidRearrangement, AssertionError):
            pass
    if not result:
        raise UserWarning(
            'row failed to create any breakpoint pairs. This generally indicates an input formatting error',
            row,
            std_row,
            combinations,
        )
    return result


def _convert_tool_output(
    input_file,
    file_type=SUPPORTED_TOOL.MAVIS,
    stranded=False,
    log=DEVNULL,
    assume_no_untemplated=True,
):
    log('reading:', input_file)
    result = []
    rows = None
    if file_type == SUPPORTED_TOOL.MAVIS:
        result = read_bpp_from_input_file(
            input_file, expand_orient=True, expand_svtype=True, add_default={'stranded': stranded}
        )
    elif file_type == SUPPORTED_TOOL.CNVNATOR:
        _, rows = tab.read_file(
            input_file,
            header=[
                'event_type',
                'coordinates',
                'size',
                'normalized_RD',
                'e-val1',
                'e-val2',
                'e-val3',
                'e-val4',
                'q0',
            ],
        )
    elif file_type in [
        SUPPORTED_TOOL.DELLY,
        SUPPORTED_TOOL.MANTA,
        SUPPORTED_TOOL.PINDEL,
        SUPPORTED_TOOL.VCF,
        SUPPORTED_TOOL.BREAKSEQ,
        SUPPORTED_TOOL.STRELKA,
        SUPPORTED_TOOL.MUTECT,
    ]:
        rows = read_vcf(input_file, file_type, log)
    elif file_type == SUPPORTED_TOOL.BREAKDANCER:
        rows = _convert_breakdancer_file(input_file)
    else:
        _, rows = tab.read_file(input_file)
    if rows:
        log('found', len(rows), 'rows')
        for row in rows:
            try:
                std_rows = _convert_tool_row(
                    row, file_type, stranded, assume_no_untemplated=assume_no_untemplated
                )
            except Exception as err:
                log('Error in converting row', row)
                raise err
            else:
                result.extend(std_rows)
    log('generated', len(result), 'breakpoint pairs')
    return result
