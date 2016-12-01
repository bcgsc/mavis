
from vocab import Vocab
from Bio.Alphabet import Gapped
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq


def reverse_complement(s):
    temp = Seq(s, DNA_ALPHABET)
    return str(temp.reverse_complement())


GAP = '-'

ORIENT = Vocab(LEFT='L', RIGHT='R', NS='?')
"""Vocab: holds controlled vocabulary for allowed orientation values

- LEFT: left wrt to the positive/forward strand
- RIGHT: right wrt to the positive/forward strand
- NS: orientation is not specified
"""
setattr(ORIENT, 'expand', lambda x: [ORIENT.LEFT, ORIENT.RIGHT] if x == ORIENT.NS else [x])
setattr(ORIENT, 'compare', lambda x, y: True if ORIENT.NS in [x, y] else (x == y))

PROTOCOL = Vocab(GENOME='genome', TRANS='transcriptome')
"""Vocab: holds controlled vocabulary for allowed protocol values

- GENOME: genome
- TRANS: transcriptome
"""

STRAND = Vocab(POS='+', NEG='-', NS='?')
"""Vocab: holds controlled vocabulary for allowed strand values

- POS: the positive/forward strand
- NEG: the negative/reverse strand
- NS: strand is not specified
"""
setattr(STRAND, 'expand', lambda x: [STRAND.POS, STRAND.NEG] if x == STRAND.NS else [x])
setattr(STRAND, 'compare', lambda x, y: True if STRAND.NS in [x, y] else (x == y))

SVTYPE = Vocab(
    DEL='deletion',
    TRANS='translocation',
    ITRANS='inverted translocation',
    INV='inversion',
    INS='insertion',
    DUP='duplication'
)
"""Vocab: holds controlled vocabulary for acceptable structural variant classifications

- DEL: deletion
- TRANS: translocation
- ITRANS: inverted translocation
- INV: inversion
- INS: insertion
- DUP: duplication
"""

CIGAR = Vocab(M=0, I=1, D=2, N=3, S=4, H=5, P=6, X=8, EQ=7)
"""Vocab: Enum-like. For readable cigar values

- M: alignment match (can be a sequence match or mismatch)
- I: insertion to the reference
- D: deletion from the reference
- N: skipped region from the reference
- S: soft clipping (clipped sequences present in SEQ)
- H: hard clipping (clipped sequences NOT present in SEQ)
- P: padding (silent deletion from padded reference)
- EQ(=): sequence match
- X: sequence mismatch

note: descriptions are taken from the samfile documentation https://samtools.github.io/hts-specs/SAMv1.pdf
"""


PYSAM_READ_FLAGS = Vocab(
    REVERSE=16,
    MATE_REVERSE=32,
    UNMAPPED=4,
    MATE_UNMAPPED=8,
    FIRST_IN_PAIR=64,
    LAST_IN_PAIR=128,
    SECONDARY=256,
    MULTIMAP=1,
    CUSTOM_REALIGN='cr'
)

"""Vocab: Enum-like. For readable PYSAM flag constants

- MULTIMAP: template having multiple segments in sequencing
- UNMAPPED: segment unmapped
- MATE_UNMAPPED: next segment in the template unmapped
- REVERSE: SEQ being reverse complemented
- MATE_REVERSE: SEQ of the next segment in the template being reverse complemented
- FIRST_IN_PAIR: the first segment in the template
- LAST_IN_PAIR: the last segment in the template
- SECONDARY: secondary alignment

note: descriptions are taken from the samfile documentation https://samtools.github.io/hts-specs/SAMv1.pdf
"""

# read paired, read mapped in proper pair, mate reverse strand, first in pair


def _match_ambiguous_dna(x, y):
    """
    >>> _match_ambiguous_dna('A', 'N')
    True
    >>> _match_ambiguous_dna('A', 'T')
    False
    >>> _match_ambiguous_dna('A', 'A')
    True
    """
    x = x.upper()
    y = y.upper()
    xset = set(ambiguous_dna_values.get(x, x))
    yset = set(ambiguous_dna_values.get(y, y))
    if len(xset.intersection(yset)) == 0:
        return False
    return True

DNA_ALPHABET = alphabet = Gapped(ambiguous_dna, '-')
DNA_ALPHABET.match = lambda x, y: _match_ambiguous_dna(x, y)

NA_MAPPING_QUALITY = 255
PHASE = Vocab(FIRST=0, SECOND=1, LAST=2, NA=-1)

FLAGS = Vocab(LQ='LOWQUAL')

READ_PAIR_TYPE = Vocab(RR='RR', LL='LL', RL='RL', LR='LR')


CALL_METHOD = Vocab(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads', MIXED='split and flanking')
