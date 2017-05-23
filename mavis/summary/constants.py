from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    filter_min_remapped_reads=5,
    filter_min_spanning_reads=5,
    filter_min_flanking_reads=5,
    filter_min_flanking_only_reads=10,
    filter_min_split_reads=5,
    filter_min_linking_split_reads=1
)
