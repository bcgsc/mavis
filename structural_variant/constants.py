
from vocab import Vocab
from Bio.Alphabet import Gapped
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Alphabet.IUPAC import ambiguous_dna

GAP = '-'

ORIENT = Vocab(LEFT='L', RIGHT='R', NS='?')

PROTOCOL = Vocab(GENOME='genome', TRANS='transcriptome')

STRAND = Vocab(POS='+', NEG='-', NS='?')

SVTYPE = Vocab(
    DEL='deletion',
    TRANS='translocation',
    ITRANS='inverted translocation',
    INV='inversion',
    INS='insertion',
    DUP='duplication'
)

"""
M 0 alignment match (can be a sequence match or mismatch)
I 1 insertion to the reference
D 2 deletion from the reference
N 3 skipped region from the reference
S 4 soft clipping (clipped sequences present in SEQ)
H 5 hard clipping (clipped sequences NOT present in SEQ)
P 6 padding (silent deletion from padded reference)
= 7 sequence match
X 8 sequence mismatch
"""
CIGAR = Vocab(M=0, I=1, D=2, N=3, S=4, H=5, P=6, X=8, EQ=7)

"""
Bit Description
1 0x1 template having multiple segments in sequencing
2 0x2 each segment properly aligned according to the aligner
4 0x4 segment unmapped
8 0x8 next segment in the template unmapped
16 0x10 SEQ being reverse complemented
32 0x20 SEQ of the next segment in the template being reverse complemented
64 0x40 the first segment in the template
128 0x80 the last segment in the template
256 0x100 secondary alignment
512 0x200 not passing filters, such as platform/vendor quality controls
1024 0x400 PCR or optical duplicate
2048 0x800 supplementary alignment
"""
PYSAM_READ_FLAGS = Vocab(
    REVERSE=16,
    MATE_REVERSE=32,
    UNMAPPED=4,
    MATE_UNMAPPED=8,
    FIRST_IN_PAIR=64,
    LAST_IN_PAIR=128
)
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
    xset = set(ambiguous_dna_values.get(x, x))
    yset = set(ambiguous_dna_values.get(y, y))
    if len(xset.intersection(yset)) == 0:
        return False
    return True

DNA_ALPHABET = alphabet = Gapped(ambiguous_dna, '-')
DNA_ALPHABET.match = lambda x, y: _match_ambiguous_dna(x, y)

NA_MAPPING_QUALITY = 255
SUFFIX_DELIM = '--'
PHASE = Vocab(FIRST=0, SECOND=1, LAST=2, NA=-1)

FLAGS = Vocab(LQ='LowQual')

READ_PAIR_TYPE = Vocab(RR='RR', LL='LL', RL='RL', LR='LR')
