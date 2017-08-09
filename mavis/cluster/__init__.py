"""
Sub-package Documentation
==========================

The cluster sub-package is responsible for merging variants coming from different inputs (i.e. different tools).

Types of Output Files
----------------------

+--------------------------------+------------------+----------------------------------------------------------------------+
| expected name/suffix           | file type/format | content                                                              |
+================================+==================+======================================================================+
| ``cluster_assignment.tab``     | text/tabbed      |                                                                      |
+--------------------------------+------------------+----------------------------------------------------------------------+
| ``uninformative_clusters.txt`` | text             | list of cluster ids that were dropped by annotation proximity filter |
+--------------------------------+------------------+----------------------------------------------------------------------+
| ``clusters.bed``               | :term:`bed`      | cluster positions                                                    |
+--------------------------------+------------------+----------------------------------------------------------------------+
| ``cluster-*.tab``              | text/tabbed      | computed clusters                                                    |
+--------------------------------+------------------+----------------------------------------------------------------------+

Algorithm Overview
--------------------

- Collapse any duplicate breakpoint pairs
- Split breakpoint pairs by type
- Cluster breakpoint pairs by distance (within a type)

    - Create a graph representation of the distances between pairs
    - Find cliques up to a given input size (cluster_clique_size)
    - Hierarchically cluster the cliques (allows redundant participation)
    - For each input node/pair pick the best cluster(s)

- Output the clusters and the mapping to the input pairs

"""
from .cluster import merge_breakpoint_pairs
