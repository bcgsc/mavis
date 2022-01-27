import re

from ..constants import MavisNamespace

PASS_FILENAME = 'annotations.tab'


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
