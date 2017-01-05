
from vocab import Vocab
from Bio.Alphabet import Gapped
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq

CODON_SIZE = 3


def reverse_complement(s):
    temp = Seq(str(s), DNA_ALPHABET)
    return str(temp.reverse_complement())


def translate(s, reading_frame=0):
    reading_frame = reading_frame % CODON_SIZE

    temp = s[reading_frame:]
    if len(temp) % 3 == 1:
        temp = temp[:-1]
    elif len(temp) % 3 == 2:
        temp = temp[:-2]
    temp = Seq(temp, DNA_ALPHABET)
    return str(temp.translate())

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

EXON_PHASE = Vocab(MIDDLE_BASE=2, FIRST_BASE=1, LAST_BASE=0, NON_CODING=-1)

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


CALL_METHOD = Vocab(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads')


SPLICE_SITE_RADIUS = 2

PRIME = Vocab(FIVE=5, THREE=3)

START_AA = 'M'
STOP_AA = '*'

# content related to tabbed files for input/output
# ensure that we don't have to change ALL the code when we update column names


class Column:
    def __init__(self, name, defn):
        self.name = name
        self.defn = defn

COLUMNS = Vocab(
    library=Column(
        'library',
        'Identifier for the library/source'),
    cluster_id=Column(
        'cluster_id',
        'Identifier for the merging/clustering step'),
    cluster_size=Column(
        'cluster_size',
        'The number of breakpoint pair calls that were grouped in creating the cluster'),
    validation_id=Column(
        'validation_id',
        'Identifier for the validation step'),
    annotation_id=Column(
        'annotation_id',
        'Identifier for the annotation step'),
    event_type=Column(
        'event_type',
        'The classification of the event'),
    gene1=Column(
        'gene1',
        'Gene for the current annotation at the first breakpoint'),
    gene2=Column(
        'gene2',
        'Gene for the current annotation at the second breakpoint'),
    transcript1=Column(
        'transcript1',
        'Transcript for the current annotation at the first breakpoint'),
    transcript2=Column(
        'transcript2',
        'Transcript for the current annotation at the second breakpoint'),
    genes_encompassed=Column(
        'genes_encompassed',
        'Applies to intrachromosomal events only. List of genes which overlap any region that occurs between both '
        'breakpoints. For example in a deletion event these would be deleted genes.'),
    genes_overlapping_break1=Column(
        'genes_overlapping_break1',
        'List of genes which overlap the first breakpoint'),
    genes_overlapping_break2=Column(
        'genes_overlapping_break2',
        'List of genes which overlap the second breakpoint'),
    genes_proximal_to_break1=Column(
        'genes_proximal_to_break1',
        'List of genes near the breakpoint and the distance away from the breakpoint'),
    genes_proximal_to_break2=Column(
        'genes_proximal_to_break2',
        'List of genes near the breakpoint and the distance away from the breakpoint'),
    break1_chromosome=Column(
        'break1_chromosome',
        'The name of the chromosome on which breakpoint 1 is situated'),
    break1_position_start=Column(
        'break1_position_start',
        'Start (inclusive, 1-based) of the range representing breakpoint 1'),
    break1_position_end=Column(
        'break1_position_end',
        'End (inclusive, 1-based) of the range representing breakpoint 1'),
    break1_orientation=Column(
        'break1_orientation',
        'The side of the breakpoint wrt the positive/forward strand that is retained'),
    break1_strand=Column(
        'break1_strand',
        'The strand wrt to the reference positive/forward strand at this breakpoint'),
    break1_sequence=Column(
        'break1_sequence',
        'The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand'),
    break2_chromosome=Column(
        'break2_chromosome',
        'The name of the chromosome on which breakpoint 2 is situated'),
    break2_position_start=Column(
        'break2_position_start',
        'Start (inclusive, 1-based) of the range representing breakpoint 2'),
    break2_position_end=Column(
        'break2_position_end',
        'End (inclusive, 1-based) of the range representing breakpoint 2'),
    break2_orientation=Column(
        'break2_orientation',
        'The side of the breakpoint wrt the positive/forward strand that is retained'),
    break2_strand=Column(
        'break2_strand',
        'The strand wrt to the reference positive/forward strand at this breakpoint'),
    break2_sequence=Column(
        'break2_sequence',
        'The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand'),
    opposing_strands=Column(
        'opposing_strands',
        'Specifies if breakpoints are on opposite strands wrt to the reference'),
    stranded=Column(
        'stranded',
        'Specifies if the sequencing protocol was strand specific or not'),
    protocol=Column(
        'protocol',
        'Genome or transcriptome'),
    tools=Column(
        'tools',
        'The tools that called the event originally (from the cluster step)'),
    contigs_assembled=Column(
        'contigs_assembled',
        'Number of contigs that were built from split read sequences'),
    contigs_aligned=Column(
        'contigs_aligned',
        'Number of contigs that were able to align'),
    contig_sequence=Column(
        'contig_sequence',
        'Sequence of the current contig (wrt to the positive forward strand if not strand specific)'),
    contig_remap_score=Column(
        'contig_remap_score',
        'Score representing the number of sequences (from the set of sequences given to the assembly algorithm) that '
        'were aligned to the resulting contig with an acceptable scoring (based on user-set thresholds). For any '
        'sequence its contribution to the score is divided by the number of mappings (to give less weight to multimaps)'
    ),
    contig_alignment_score=Column(
        'contig_alignment_score',
        'A rank based on the alignment tool (blat), etc.) of the alignment being used. An average if split alignments '
        'were used. Lower numbers indicate a better alignment. If it was the base alignment possible then this would be'
        'zero'),
    break1_call_method=Column(
        'break1_call_method',
        'The method used to call the first breakpoint'),
    break2_call_method=Column(
        'break2_call_method',
        'The method used to call the second breakpoint'),
    flanking_reads=Column(
        'flanking_reads',
        'Number of read-pairs where one read aligns to the first breakpoint window and the second read aligns to the '
        'other. The count here is based on the number of unique query names'),
    median_insert_size=Column(
        'median_insert_size',
        'The median insert size of the flanking reads being used as evidence'),
    stdev_insert_size=Column(
        'stdev_insert_size',
        'The standard deviation in insert size of the flanking reads being used as evidence'),
    break1_split_reads=Column(
        'break1_split_reads',
        'Number of split reads that call the exact breakpoint given'),
    break1_split_reads_forced=Column(
        'break1_split_reads_forced',
        'Number of split reads which were re-aligned to the opposite breakpoint window'),
    break2_split_reads=Column(
        'break2_split_reads',
        'Number of split reads that call the exact breakpoint given'),
    break2_split_reads_forced=Column(
        'break2_split_reads_forced',
        'Number of split reads which were re-aligned to the opposite breakpoint window'),
    linking_split_reads=Column(
        'linking_split_reads',
        'Number of split reads that align to both breakpoints'),
    untemplated_sequence=Column(
        'untemplated_sequence',
        'The untemplated/novel sequence between the breakpoints'),
    break1_homologous_sequence=Column(
        'break1_homologous_sequence',
        'Sequence in common at the first breakpoint and other side of the second breakpoint'),
    break2_homologous_sequence=Column(
        'break2_homologous_sequence',
        'Sequence in common at the second breakpoint and other side of the first breakpoint'),
    break1_ewindow=Column(
        'break1_ewindow',
        'Window where evidence was gathered for the first breakpoint'),
    break1_ewindow_count=Column(
        'break1_ewindow_count',
        'Number of reads processed/looked-at in the first evidence window'),
    break2_ewindow=Column(
        'break2_ewindow',
        'Window where evidence was gathered for the second breakpoint'),
    break2_ewindow_count=Column(
        'break2_ewindow_count',
        'Number of reads processed/looked-at in the second evidence window'),
    raw_flanking_reads=Column(
        'raw_flanking_reads',
        'Number of flanking reads before calling the breakpoint. The count here is based on the number of unique query '
        'names'),
    raw_break1_split_reads=Column(
        'raw_break1_split_reads',
        'Number of split reads before calling the breakpoint'),
    raw_break2_split_reads=Column(
        'raw_break2_split_reads',
        'Number of split reads before calling the breakpoint')
)


def define_column(col):
    for c in COLUMNS.values():
        if c.name == col:
            return c.defn
    raise KeyError('column name not found', col)


def sort_columns(input_columns):
    order = {}
    for i, col in enumerate(COLUMNS.values()):
        order[col.name] = i

    temp = sorted([c for c in input_columns if c in order], key=lambda x: order[x])
    temp = temp + sorted([c for c in input_columns if c not in order])
    return temp
