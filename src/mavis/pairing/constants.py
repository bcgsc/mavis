from typing import Dict

from mavis_config import DEFAULTS

from ..constants import CALL_METHOD

PAIRING_DISTANCES: Dict[str, int] = {
    CALL_METHOD.FLANK: DEFAULTS['pairing.flanking_call_distance'],
    CALL_METHOD.SPAN: DEFAULTS['pairing.spanning_call_distance'],
    CALL_METHOD.SPLIT: DEFAULTS['pairing.split_call_distance'],
    CALL_METHOD.CONTIG: DEFAULTS['pairing.contig_call_distance'],
    CALL_METHOD.INPUT: DEFAULTS['pairing.input_call_distance'],
}
