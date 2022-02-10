from ..constants import SVTYPE, MavisNamespace


class SUPPORTED_TOOL(MavisNamespace):
    """
    Supported Tools used to call SVs and then used as input into MAVIS

    Attributes:
        CHIMERASCAN: chimerascan [Iyer-2011]_
        DEFUSE: defuse [McPherson-2011]_
        DELLY: delly [Rausch-2012]_
        MANTA: manta [Chen-2016]_
        PINDEL: pindel [Ye-2009]_
        TA: transabyss [Robertson-2010]_
    """

    MANTA = 'manta'
    DELLY = 'delly'
    TA = 'transabyss'
    PINDEL = 'pindel'
    CHIMERASCAN = 'chimerascan'
    MAVIS = 'mavis'
    DEFUSE = 'defuse'
    BREAKDANCER = 'breakdancer'
    VCF = 'vcf'
    BREAKSEQ = 'breakseq'
    CNVNATOR = 'cnvnator'
    STRELKA = 'strelka'
    STARFUSION = 'starfusion'


TOOL_SVTYPE_MAPPING = {v: [v] for v in SVTYPE.values()}  # type: ignore
TOOL_SVTYPE_MAPPING.update(
    {
        'DEL': [SVTYPE.DEL],
        'INS': [SVTYPE.INS],
        'ITX': [SVTYPE.DUP],
        'CTX': [SVTYPE.TRANS, SVTYPE.ITRANS],
        'INV': [SVTYPE.INV],
        'BND': [SVTYPE.TRANS, SVTYPE.ITRANS, SVTYPE.DUP, SVTYPE.INS, SVTYPE.DEL, SVTYPE.INV],
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
        'ITD': [SVTYPE.DUP],
        'IDP': [SVTYPE.INS],
        'DEL/INV': [SVTYPE.DEL, SVTYPE.INV],
        'DUP/INS': [SVTYPE.DUP, SVTYPE.INS],
        'INVDUP': [SVTYPE.INV, SVTYPE.DUP, SVTYPE.INS],
        'INV/INVDUP': [SVTYPE.INV, SVTYPE.DUP, SVTYPE.INS],
    }
)

TRACKING_COLUMN = 'tracking_id'
