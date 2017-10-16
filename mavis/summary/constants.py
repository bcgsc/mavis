from ..constants import MavisNamespace
from ..util import WeakMavisNamespace


DEFAULTS = WeakMavisNamespace()
"""
- :term:`filter_cdna_synon`
- :term:`filter_min_flanking_reads`
- :term:`filter_min_linking_split_reads`
- :term:`filter_min_remapped_reads`
- :term:`filter_min_spanning_reads`
- :term:`filter_min_split_reads`
- :term:`filter_protein_synon`
"""
DEFAULTS.add('filter_min_remapped_reads', 5, defn='Minimum number of remapped reads for a call by contig')
DEFAULTS.add('filter_min_spanning_reads', 5, defn='Minimum number of spanning reads for a call by spanning reads')
DEFAULTS.add('filter_min_flanking_reads', 10, defn='Minimum number of flanking pairs for a call by flanking pairs')
DEFAULTS.add('filter_min_split_reads', 5, defn='Minimum number of split reads for a call by split reads')
DEFAULTS.add('filter_min_linking_split_reads', 1, defn='Minimum number of linking split reads for a call by split reads')
DEFAULTS.add('filter_cdna_synon', True, defn='Filter all annotations synonymous at the cdna level')
DEFAULTS.add('filter_protein_synon', True, defn='Filter all annotations synonymous at the protein level')


PAIRING_STATE = MavisNamespace(
    EXP='expressed',
    NO_EXP='not expressed',
    SOMATIC='somatic',
    GERMLINE='germline',
    CO_EXP='co-expressed',
    GERMLINE_EXP='germline expression',
    SOMATIC_EXP='somatic expression',
    MATCH='matched',
    NO_MATCH='not matched',
    GENOMIC='genomic support',
    NO_GENOMIC='no genomic support'
)
