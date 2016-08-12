
from vocab import Vocab
from Bio.Alphabet import Gapped
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Alphabet.IUPAC import ambiguous_dna

GAP = '-'

ORIENT = Vocab(LEFT='L', RIGHT='R', NS='?')

STRAND = Vocab(POS='+', NEG='-', NS='?')

SVTYPE = Vocab(DEL='deletion', TRANS='translocation', ITRANS='inverted translocation', 
        INV='inversion', INS='insertion', DUP='duplication', COMPLEX='complex', NS='not specified')

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
