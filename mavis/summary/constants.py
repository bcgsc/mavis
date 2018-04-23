from ..constants import MavisNamespace, float_fraction
from ..util import WeakMavisNamespace


DEFAULTS = WeakMavisNamespace()
HOMOPOLYMER_MIN_LENGTH = 3

"""
- :term:`filter_cdna_synon`
- :term:`filter_min_flanking_reads`
- :term:`filter_min_linking_split_reads`
- :term:`filter_min_remapped_reads`
- :term:`filter_min_spanning_reads`
- :term:`filter_min_split_reads`
- :term:`filter_protein_synon`
- :term:`filter_min_complexity`
- :term:`filter_trans_homopolymers`
"""
DEFAULTS.add('filter_min_remapped_reads', 5, defn='Minimum number of remapped reads for a call by contig')
DEFAULTS.add('filter_min_spanning_reads', 5, defn='Minimum number of spanning reads for a call by spanning reads')
DEFAULTS.add('filter_min_flanking_reads', 10, defn='Minimum number of flanking pairs for a call by flanking pairs')
DEFAULTS.add('filter_min_split_reads', 5, defn='Minimum number of split reads for a call by split reads')
DEFAULTS.add('filter_min_linking_split_reads', 1, defn='Minimum number of linking split reads for a call by split reads')
DEFAULTS.add('filter_cdna_synon', True, defn='Filter all annotations synonymous at the cdna level')
DEFAULTS.add('filter_protein_synon', False, defn='Filter all annotations synonymous at the protein level')
DEFAULTS.add(
    'filter_trans_homopolymers', True,
    defn='Filter all single bp ins/del/dup events that are in a homopolymer region of at least '
         '{} bps and are not paired to a genomic event'.format(HOMOPOLYMER_MIN_LENGTH))
DEFAULTS.add(
    'filter_min_complexity', 0.2, cast_type=float_fraction,
    defn='Filter event calls based on call sequence complexity')


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
