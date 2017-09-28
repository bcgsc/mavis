import re

from vocab import Vocab

from ..util import MavisNamespace


DEFAULTS = MavisNamespace(
    min_domain_mapping_match=0.9,
    min_orf_size=300,
    max_orf_cap=3,
    annotation_filters='choose_more_annotated,choose_transcripts_by_priority'
)


SPLICE_TYPE = Vocab(
    RETAIN='retained intron',
    SKIP='skipped exon',
    NORMAL='normal',
    MULTI_RETAIN='retained multiple introns',
    MULTI_SKIP='skipped multiple exons',
    COMPLEX='complex'
)
""":class:`Vocab`: holds controlled vocabulary for allowed splice type classification values

- ``RETAIN``: an intron was retained
- ``SKIP``: an exon was skipped
- ``NORMAL``: no exons were skipped and no introns were retained. the normal/expected splicing pattern was followed
- ``MULTI_RETAIN``: multiple introns were retained
- ``MULTI_SKIP``: multiple exons were skipped
- ``COMPLEX``: some combination of exon skipping and intron retention
"""

SPLICE_SITE_TYPE = Vocab(
    DONOR=3,
    ACCEPTOR=5
)

SPLICE_SITE_RADIUS = 2
""":class:`int`: number of bases away from an exon boundary considered to be part of the splice site such that if it were altered
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
    re.compile('([TC]{8}[ATCG]AAG)([GA][ATCG])')
]
