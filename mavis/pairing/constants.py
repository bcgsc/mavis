from ..util import MavisNamespace
from vocab import Vocab


DEFAULTS = MavisNamespace(
    flanking_call_distance=0,
    split_call_distance=10,
    contig_call_distance=0,
    spanning_call_distance=5
)

PAIRING_STATE = Vocab(
    EXP='expressed',
    NO_EXP='not expressed',
    SOMATIC='somatic',
    GERMLINE='germline',
    CO_EXP='co-expressed',
    GERMLINE_EXP='germline expression',
    SOMATIC_EXP='somatic expression',
    MATCH='matched',
    NO_MATCH='not matched',
    GENOMIC='genomic support',
    NO_GENOMIC='no genomic support'
)
