import re

import tab

from ..constants import MavisNamespace, float_fraction

PASS_FILENAME = 'annotations.tab'


class SPLICE_TYPE(MavisNamespace):
    """
    holds controlled vocabulary for allowed splice type classification values

    Attributes:
        RETAIN: an intron was retained
        SKIP: an exon was skipped
        NORMAL: no exons were skipped and no introns were retained. the normal/expected splicing pattern was followed
        MULTI_RETAIN: multiple introns were retained
        MULTI_SKIP: multiple exons were skipped
        COMPLEX: some combination of exon skipping and intron retention
    """

    RETAIN: str = 'retained intron'
    SKIP: str = 'skipped exon'
    NORMAL: str = 'normal'
    MULTI_RETAIN: str = 'retained multiple introns'
    MULTI_SKIP: str = 'skipped multiple exons'
    COMPLEX: str = 'complex'


class SPLICE_SITE_TYPE(MavisNamespace):
    DONOR: int = 3
    ACCEPTOR: int = 5


SPLICE_SITE_RADIUS = 2
"""int: number of bases away from an exon boundary considered to be part of the splice site such that if it were altered
        the splice site would be considered to be abrogated.
"""

# splice site sequences based on: http://www.nature.com/nrg/journal/v17/n7/fig_tab/nrg.2016.46_F5.html?foxtrotcallback=true

DONOR_SEQ = [
    re.compile('(AG)(GT[AG]AG)'),
    re.compile('([CA]AG)(GTA)'),
]

ACCEPTOR_SEQ = [
    re.compile('([TC]{8}[ATCG]CAG)([GA][ATCG])'),
    re.compile('([TC]{9}TAG)([GA][ATCG])'),
    re.compile('([TC]{8}[ATCG]AAG)([GA][ATCG])'),
]