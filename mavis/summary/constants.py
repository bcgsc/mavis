from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    filter_min_contig_alignment_score=5,
    filter_min_flanking_pairs=4,
    filter_min_linking_split_reads=2,
    filter_min_spanning_reads=3,
    filter_min_split_reads=3
)
