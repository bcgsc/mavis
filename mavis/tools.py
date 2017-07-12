import TSV
import re
import vcf
import itertools
import warnings
from .constants import STRAND, SVTYPE, ORIENT
from vocab import Vocab
from .breakpoint import read_bpp_from_input_file, Breakpoint, BreakpointPair
from .error import InvalidRearrangement
from .util import devnull

SUPPORTED_TOOL = Vocab(
    MANTA='manta',
    DELLY='delly',
    TA='transabyss',
    BREAKDANCER='breakdancer',
    PINDEL='pindel',
    CHIMERASCAN='chimerascan',
    MAVIS='mavis',
    DEFUSE='defuse'
)


def default_chr_filter(chr):
    ok = [str(x) for x in range(1, 23)]
    ok.extend(['X', 'Y', 'MT'])
    ok.extend(['chr' + x for x in ok])
    return True if chr in ok else False


def convert_tool_output(input_file, file_type=SUPPORTED_TOOL.MAVIS, stranded=False, log=devnull):
    log('reading:', input_file)
    result = []
    if file_type == SUPPORTED_TOOL.TA:
        warnings.warn('currently assuming that trans-abyss is calling the strand exactly opposite and swapping them')
    if file_type == SUPPORTED_TOOL.MAVIS:
        result = read_bpp_from_input_file(input_file)
    else:
        if file_type in [SUPPORTED_TOOL.DELLY, SUPPORTED_TOOL.MANTA]:
            rows = vcf.Reader(filename=input_file)
        else:
            header, rows = TSV.read_file(input_file)
        log('found', len(rows), 'rows')
        for row in rows:
            types = {k: [v] for k, v in SVTYPE.items()}
            types.update({
                'DEL': [SVTYPE.DEL],
                'INS': [SVTYPE.INS],
                'ITX': [SVTYPE.DUP],
                'CTX': [SVTYPE.TRANS, SVTYPE.ITRANS],
                'INV': [SVTYPE.INV],
                'BND': [SVTYPE.TRANS],
                'TRA': [SVTYPE.TRANS, SVTYPE.ITRANS],
                'CNV': [SVTYPE.DUP],
                'RPL': [SVTYPE.INS],
                'DUP:TANDEM': [SVTYPE.DUP],
                'interchromosomal': [SVTYPE.TRANS, SVTYPE.ITRANS],
                'eversion': [SVTYPE.DUP],
                'translocation': [SVTYPE.TRANS, SVTYPE.ITRANS]
            })
            std_row = {}
            orient1 = orient2 = ORIENT.NS
            strand1 = strand2 = STRAND.NS
            event_type = None

            # convert the specified file type to a standard format

            if file_type == SUPPORTED_TOOL.DELLY or file_type == SUPPORTED_TOOL.MANTA:

                if row.INFO['SVTYPE'] == 'BND':
                    chr2 = row.ALT[0].chr
                    end = row.ALT[0].pos
                else:
                    chr2 = row.INFO.get('CHR2', row.CHROM)
                    end = row.INFO.get('END', row.POS)

                std_row.update({
                    'chr1': row.CHROM, 'chr2': row.INFO['CHR2'],
                    'pos1_start': max(1, row.POS + row.INFO.get('CIPOS', (0, 0))[0]),
                    'pos1_end': row.POS + row.INFO.get('CIPOS', (0, 0))[1],
                    'pos2_start': max(1, row.INFO['END'] + row.get('CIEND', (0, 0))[0]),
                    'pos2_end': row.INFO['END'] + row.INFO.get('CIEND', (0, 0))[1]
                })
                event_type = row.INFO['SVTYPE']
                try:
                    o1, o2 = row.INFO['CT'].split('to')
                    CT = {'3': ORIENT.LEFT, '5': ORIENT.RIGHT, 'N': ORIENT.NS}
                    orient1 = CT[o1]
                    orient2 = CT[o2]
                except KeyError:
                    pass

            elif file_type == SUPPORTED_TOOL.BREAKDANCER:
                std_row.update({
                    'chr1': row['Chr1'], 'chr2': row['Chr2'],
                    'pos1_start': row['Pos1'], 'pos2_start': row['Pos2'],
                    'event_type': row['Type']
                })

            elif file_type == SUPPORTED_TOOL.PINDEL:
                info = {k: v for k, v in ['='.split(x) for x in ';'.split(row['info'])]}
                std_row.update({'chr1': row['CHROM'], 'pos1_start': row['POS'], 'pos2_start': info['END']})

            elif file_type == SUPPORTED_TOOL.DEFUSE:
                orient1 = ORIENT.LEFT if row['genomic_strand1'] == STRAND.POS else ORIENT.RIGHT
                orient2 = ORIENT.LEFT if row['genomic_strand2'] == STRAND.POS else ORIENT.RIGHT
                std_row.update({
                    'chr1': row['gene_chromosome1'], 'chr2': row['gene_chromosome2'],
                    'pos1_start': row['genomic_break_pos1'], 'pos2_start': row['genomic_break_pos2']
                })

            elif file_type == SUPPORTED_TOOL.TA:
                std_row['event_type'] = row['rearrangement']
                if 'breakpoint' in row:
                    orient1, orient2 = row['orientations'].split(',')
                    if stranded:
                        strand1, strand2 = row['strands'].split(',')
                        strand1 = STRAND.POS if strand1 == STRAND.NEG else STRAND.POS
                        strand2 = STRAND.POS if strand2 == STRAND.NEG else STRAND.POS
                    m = re.match(
                        '^(?P<chr1>[^:]+):(?P<pos1_start>\d+)\|(?P<chr2>[^:]+):(?P<pos2_start>\d+)$', row['breakpoint'])
                    if not m:
                        raise OSError(
                            'file format error: the breakpoint column did not satisfy the expected pattern', row)
                    std_row.update({k: m[k] for k in ['chr1', 'pos1_start', 'chr2', 'pos1_start']})
                else:
                    std_row.update({
                        'chr1': row['chr'], 'pos1_start': row['chr_start'], 'pos2_start': row['chr_end'] + 1})
                    if stranded:
                        strand1 = strand2 = row['ctg_strand']
            else:
                raise NotImplementedError('unsupported file type', file_type)

            if stranded:
                strand1 = STRAND.expand(strand1)
                strand2 = STRAND.expand(strand2)
            else:
                strand1 = [strand1]
                strand2 = [strand2]

            one_or_more = False
            combinations = itertools.product(
                ORIENT.expand(orient1), ORIENT.expand(orient2),
                strand1, strand2, types[row['event_type']] if 'event_type' in row else SVTYPE.values(), [True, False]
            )
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
                            std_row['chr1'],
                            std_row['pos1_start'],
                            std_row.get('pos1_end', std_row['pos1_start']),
                            orient=orient1, strand=strand1
                        ),
                        opposing_strands=oppose,
                        event_type=event_type
                    )
                    if event_type in BreakpointPair.classify(bpp):
                        result.append(bpp)
                        one_or_more = True
                except (InvalidRearrangement, AssertionError):
                    pass
            if not one_or_more:
                raise UserWarning(
                    'row failed to create any breakpoint pairs. This generally indicates an input formatting error',
                    row, combinations)
    log('generated', len(result), 'breakpoint pairs')
    return result
