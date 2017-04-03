from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    min_clusters_per_file=50,
    max_files=10,
    cluster_clique_size=15,
    cluster_radius=20,
    max_proximity=5000,
    uninformative_filter=True
)
