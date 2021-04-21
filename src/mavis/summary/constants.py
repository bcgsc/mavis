from ..constants import MavisNamespace

HOMOPOLYMER_MIN_LENGTH = 3


class PAIRING_STATE(MavisNamespace):
    EXP = 'expressed'
    NO_EXP = 'not expressed'
    SOMATIC = 'somatic'
    GERMLINE = 'germline'
    CO_EXP = 'co-expressed'
    GERMLINE_EXP = 'germline expression'
    SOMATIC_EXP = 'somatic expression'
    MATCH = 'matched'
    NO_MATCH = 'not matched'
    GENOMIC = 'genomic support'
    NO_GENOMIC = 'no genomic support'
