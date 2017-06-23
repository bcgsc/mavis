from ..util import MavisNamespace
from vocab import Vocab


DEFAULTS = MavisNamespace(
    flanking_call_distance=0,
    split_call_distance=10,
    contig_call_distance=0,
    spanning_call_distance=5
)
""":class:`MavisNamespace`: default settings for the pairing module

.. glossary::
    :sorted:

    flanking_call_distance:
        the maximum distance allowed between breakpoint pairs (called by flanking pairs) in order for them to pair 

    split_call_distance:
        the maximum distance allowed between breakpoint pairs (called by split reads) in order for them to pair

    contig_call_distance:
        the maximum distance allowed between breakpoint pairs (called by contig) in order for them to pair 

    spanning_call_distance:
        the maximum distance allowed between breakpoint pairs (called by spanning reads) in order for them to pair
"""

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
