from ..constants import MavisNamespace, SVTYPE

SUPPORTED_TOOL = MavisNamespace(
    MANTA='manta',
    DELLY='delly',
    TA='transabyss',
    PINDEL='pindel',
    CHIMERASCAN='chimerascan',
    MAVIS='mavis',
    DEFUSE='defuse',
    BREAKDANCER='breakdancer',
    VCF='vcf',
    BREAKSEQ='breakseq',
    CNVNATOR='cnvnator',
    STRELKA='strelka',
    STARFUSION='starfusion',
    MUTECT='mutect',
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
    }
)

TRACKING_COLUMN = 'tracking_id'
