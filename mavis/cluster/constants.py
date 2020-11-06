from ..util import WeakMavisNamespace


DEFAULTS = WeakMavisNamespace()
"""
- [cluster_initial_size_limit](/configuration/settings/#cluster_initial_size_limit)
- [cluster_radius](/configuration/settings/#cluster_radius)
- [limit_to_chr](/configuration/settings/#limit_to_chr)
- [max_files](/configuration/settings/#max_files)
- [max_proximity](/configuration/settings/#max_proximity)
- [min_clusters_per_file](/configuration/settings/#min_clusters_per_file)
- [uninformative_filter](/configuration/settings/#uninformative_filter)
"""
DEFAULTS.add(
    'min_clusters_per_file', 50, defn='the minimum number of breakpoint pairs to output to a file'
)
DEFAULTS.add(
    'max_files', 200, defn='The maximum number of files to output from clustering/splitting'
)
DEFAULTS.add(
    'cluster_initial_size_limit',
    25,
    defn='the maximum cumulative size of both breakpoints for breakpoint pairs to be used in the initial clustering '
    'phase (combining based on overlap)',
)
DEFAULTS.add('cluster_radius', 100, defn='maximum distance allowed between paired breakpoint pairs')
DEFAULTS.add(
    'max_proximity',
    5000,
    defn='the maximum distance away from an annotation before the region in considered to be uninformative',
)
DEFAULTS.add(
    'uninformative_filter',
    False,
    defn='flag that determines if breakpoint pairs which are not within max_proximity to any annotations are filtered '
    'out prior to clustering',
)
DEFAULTS.add(
    'limit_to_chr',
    [str(x) for x in range(1, 23)] + ['X', 'Y'],
    cast_type=str,
    listable=True,
    nullable=True,
    defn='A list of chromosome names to use. BreakpointPairs on other chromosomes will be filtered'
    'out. For example \'1 2 3 4\' would filter out events/breakpoint pairs on any chromosomes but 1, 2, 3, and 4',
)
