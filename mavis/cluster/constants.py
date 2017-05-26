from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    min_clusters_per_file=50,
    max_files=100,
    cluster_clique_size=10,
    cluster_radius=100,
    max_proximity=5000,
    uninformative_filter=True
)
