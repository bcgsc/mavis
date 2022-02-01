"""
Helper classes for type hints
"""

from typing import TYPE_CHECKING, Dict, List, Tuple

from Bio.SeqRecord import SeqRecord

if TYPE_CHECKING:
    from .annotate.genomic import Gene

ReferenceGenome = Dict[str, SeqRecord]
ReferenceAnnotations = Dict[str, List['Gene']]
CigarTuples = List[Tuple[int, int]]
