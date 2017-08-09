from ..util import MavisNamespace, ChrListString


DEFAULTS = MavisNamespace(
    min_clusters_per_file=50,
    max_files=100,
    cluster_initial_size_limit=25,
    cluster_radius=100,
    max_proximity=5000,
    uninformative_filter=True,
    limit_to_chr=ChrListString(';'.join([str(x) for x in range(1, 23)] + ['X', 'Y']))
)
