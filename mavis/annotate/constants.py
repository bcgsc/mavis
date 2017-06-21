from ..util import MavisNamespace
from .variant import choose_more_annotated, choose_transcripts_by_priority

ACCEPTED_FILTERS = {
    'choose_more_annotated': choose_more_annotated,
    'choose_transcripts_by_priority': choose_transcripts_by_priority
}

DEFAULTS = MavisNamespace(
    min_domain_mapping_match=0.9,
    min_orf_size=300,
    max_orf_cap=3,
    annotation_filters='choose_more_annotated,choose_transcripts_by_priority'
)
