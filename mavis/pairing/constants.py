from ..constants import CALL_METHOD, MavisNamespace
from ..util import WeakMavisNamespace


DEFAULTS = WeakMavisNamespace()
"""
- :term:`contig_call_distance`
- :term:`flanking_call_distance`
- :term:`spanning_call_distance`
- :term:`split_call_distance`
"""
DEFAULTS.add(
    'flanking_call_distance', 50,
    defn='the maximum distance allowed between breakpoint pairs (called by flanking pairs) in order for them to pair')
DEFAULTS.add(
    'split_call_distance', 20,
    defn='the maximum distance allowed between breakpoint pairs (called by split reads) in order for them to pair')
DEFAULTS.add(
    'contig_call_distance', 10,
    defn='the maximum distance allowed between breakpoint pairs (called by contig) in order for them to pair')
DEFAULTS.add(
    'spanning_call_distance', 20,
    defn='the maximum distance allowed between breakpoint pairs (called by spanning reads) in order for them to pair')
DEFAULTS.add(
    'input_call_distance', 20,
    defn='the maximum distance allowed between breakpoint pairs (called by input tools, not validated) in order for them to pair')

PAIRING_DISTANCES = MavisNamespace(**{
    CALL_METHOD.FLANK: DEFAULTS.flanking_call_distance,
    CALL_METHOD.SPAN: DEFAULTS.spanning_call_distance,
    CALL_METHOD.SPLIT: DEFAULTS.split_call_distance,
    CALL_METHOD.CONTIG: DEFAULTS.contig_call_distance,
    CALL_METHOD.INPUT: DEFAULTS.input_call_distance
})
