from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    min_contig_alignment_score=5,
    min_flanking_pairs=4,
    min_linking_split_reads=2,
    min_spanning_reads=3,
    min_split_reads=3
)
