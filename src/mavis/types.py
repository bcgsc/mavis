"""
Helper classes for type hints
"""

from typing import Dict, List, Tuple

from Bio.SeqRecord import SeqRecord

ReferenceGenome = Dict[str, SeqRecord]

CigarTuples = List[Tuple[int, int]]
