import re

import tab

from ..constants import MavisNamespace, float_fraction
from ..util import WeakMavisNamespace


PASS_FILENAME = 'annotations.tab'

DEFAULTS = WeakMavisNamespace()
"""
- :term:`annotation_filters`
- :term:`max_orf_cap`
- :term:`min_domain_mapping_match`
- :term:`min_orf_size`
"""
DEFAULTS.add(
    'min_domain_mapping_match', 0.9, cast_type=float_fraction,
    defn='a number between 0 and 1 representing the minimum percent match a domain must map to the fusion transcript '
    'to be displayed')
DEFAULTS.add('min_orf_size', 300, defn='the minimum length (in base pairs) to retain a putative open reading frame (ORF)')
DEFAULTS.add('max_orf_cap', 3, defn='the maximum number of ORFs to return (best putative ORFs will be retained)')
DEFAULTS.add(
    'annotation_filters', 'choose_more_annotated,choose_transcripts_by_priority',
    defn='a comma separated list of filters to apply to putative annotations'
)
DEFAULTS.add(
    'draw_fusions_only', True, cast_type=tab.cast_boolean,
    defn='flag to indicate if events which do not produce a fusion transcript should produce illustrations')
DEFAULTS.add(
    'draw_non_synonymous_cdna_only', True, cast_type=tab.cast_boolean,
    defn='flag to indicate if events which are synonymous at the cdna level should produce illustrations')

SPLICE_TYPE = MavisNamespace(
    RETAIN='retained intron',
    SKIP='skipped exon',
    NORMAL='normal',
    MULTI_RETAIN='retained multiple introns',
    MULTI_SKIP='skipped multiple exons',
    COMPLEX='complex'
)
""":class:`MavisNamespace`: holds controlled vocabulary for allowed splice type classification values

- ``RETAIN``: an intron was retained
- ``SKIP``: an exon was skipped
- ``NORMAL``: no exons were skipped and no introns were retained. the normal/expected splicing pattern was followed
- ``MULTI_RETAIN``: multiple introns were retained
- ``MULTI_SKIP``: multiple exons were skipped
- ``COMPLEX``: some combination of exon skipping and intron retention
"""

SPLICE_SITE_TYPE = MavisNamespace(
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
