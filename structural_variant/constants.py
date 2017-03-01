"""
module responsible for small utility functions and constants used throughout the structural_variant package
"""
from vocab import Vocab
from Bio.Alphabet import Gapped
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq
from datetime import datetime
import random
import math

CODON_SIZE = 3
""":class:`int`: the number of bases making up a codon"""


def reverse_complement(s):
    """
    wrapper for the Bio.Seq reverse_complement method

    Args:
        s (str): the input DNA sequence

    Returns:
        :class:`str`: the reverse complement of the input sequence

    Warning:
        assumes the input is a DNA sequence

    Example:
        >>> reverse_complement('ATCCGGT')
        'ACCGGAT'
    """
    temp = Seq(str(s), DNA_ALPHABET)
    return str(temp.reverse_complement())


def build_batch_id(prefix='', suffix='', size=6):
    date = datetime.now()
    m = int(math.pow(10, size) - 1)
    return 'batch{prefix}{date.year}{date.month:02d}{date.day:02d}r{r:06d}{suffix}'.format(
        prefix=prefix, suffix=suffix, date=date, r=random.randint(1, m))


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def translate(s, reading_frame=0):
    """
    given a DNA sequence, translates it and returns the protein amino acid sequence

    Args:
        s (str): the input DNA sequence
        reading_frame (int): where to start translating the sequence

    Returns:
        str: the amino acid sequence
    """
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
""":class:`Vocab`: holds controlled vocabulary for allowed orientation values

- ``LEFT``: left wrt to the positive/forward strand
- ``RIGHT``: right wrt to the positive/forward strand
- ``NS``: orientation is not specified
"""
setattr(ORIENT, 'expand', lambda x: [ORIENT.LEFT, ORIENT.RIGHT] if x == ORIENT.NS else [x])
setattr(ORIENT, 'compare', lambda x, y: True if ORIENT.NS in [x, y] else (x == y))

PROTOCOL = Vocab(GENOME='genome', TRANS='transcriptome')
""":class:`Vocab`: holds controlled vocabulary for allowed protocol values

- ``GENOME``: genome
- ``TRANS``: transcriptome
"""

STRAND = Vocab(POS='+', NEG='-', NS='?')
""":class:`Vocab`: holds controlled vocabulary for allowed strand values

- ``POS``: the positive/forward strand
- ``NEG``: the negative/reverse strand
- ``NS``: strand is not specified
"""
setattr(STRAND, 'expand', lambda x: [STRAND.POS, STRAND.NEG] if x == STRAND.NS else [x])
setattr(STRAND, 'compare', lambda x, y: True if STRAND.NS in [x, y] else (x == y))

SVTYPE = Vocab(
    DEL='deletion',
    TRANS='translocation',
    ITRANS='inverted translocation',
    INV='inversion',
    INS='insertion',
    DUP='duplication',
    INDEL='indel'
)
""":class:`Vocab`: holds controlled vocabulary for acceptable structural variant classifications

- ``DEL``: deletion
- ``TRANS``: translocation
- ``ITRANS``: inverted translocation
- ``INV``: inversion
- ``INS``: insertion
- ``DUP``: duplication
"""

CIGAR = Vocab(M=0, I=1, D=2, N=3, S=4, H=5, P=6, X=8, EQ=7)
""":class:`Vocab`: Enum-like. For readable cigar values

- ``M``: alignment match (can be a sequence match or mismatch)
- ``I``: insertion to the reference
- ``D``: deletion from the reference
- ``N``: skipped region from the reference
- ``S``: soft clipping (clipped sequences present in SEQ)
- ``H``: hard clipping (clipped sequences NOT present in SEQ)
- ``P``: padding (silent deletion from padded reference)
- ``EQ``: sequence match (=)
- ``X``: sequence mismatch

note: descriptions are taken from the `samfile documentation <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
"""

NA_MAPPING_QUALITY = 255
""":class:`int`: mapping qaulity value to indicate mapping was not performed/calculated"""

PYSAM_READ_FLAGS = Vocab(
    REVERSE=16,
    MATE_REVERSE=32,
    UNMAPPED=4,
    MATE_UNMAPPED=8,
    FIRST_IN_PAIR=64,
    LAST_IN_PAIR=128,
    SECONDARY=256,
    MULTIMAP=1,
    SUPPLEMENTARY=2048,
    TARGETED_ALIGNMENT='ta',
    RECOMPUTED_CIGAR='rc',
    BLAT_RANK='br',
    BLAT_SCORE='bs',
    BLAT_ALIGNMENTS='ba',
    BLAT_PERCENT_IDENTITY='bi',
    BLAT_PMS='bp'
)

""":class:`Vocab`: Enum-like. For readable PYSAM flag constants

- ``MULTIMAP``: template having multiple segments in sequencing
- ``UNMAPPED``: segment unmapped
- ``MATE_UNMAPPED``: next segment in the template unmapped
- ``REVERSE``: SEQ being reverse complemented
- ``MATE_REVERSE``: SEQ of the next segment in the template being reverse complemented
- ``FIRST_IN_PAIR``: the first segment in the template
- ``LAST_IN_PAIR``: the last segment in the template
- ``SECONDARY``: secondary alignment
- ``SUPPLEMENTARY``: supplementary alignment

note: descriptions are taken from the `samfile documentation <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
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

FLAGS = Vocab(LQ='LOWQUAL')

READ_PAIR_TYPE = Vocab(RR='RR', LL='LL', RL='RL', LR='LR')

CALL_METHOD = Vocab(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads')
""":class:`Vocab`: holds controlled vocabulary for allowed call methods

- ``CONTIG``: a contig was assembled and aligned across the breakpoints
- ``SPLIT``: the event was called by split reads
- ``FLANK``: the event was called by flanking reads
"""

GENE_PRODUCT_TYPE = Vocab(SENSE='sense', ANTI_SENSE='anti-sense')
""":class:`Vocab`: controlled vocabulary for gene products

- ``SENSE``: the gene product is a sense fusion
- ``ANTI_SENSE``: the gene product is anti-sense
"""

SPLICE_TYPE = Vocab(
    RETAIN='retained intron',
    SKIP='skipped exon',
    NORMAL='normal',
    MULTI_RETAIN='retained multiple introns',
    MULTI_SKIP='skipped multiple exons',
    COMPLEX='complex'
)
""":class:`Vocab`: holds controlled vocabulary for allowed splice type classification values

- ``RETAIN``: an intron was retained
- ``SKIP``: an exon was skipped
- ``NORMAL``: no exons were skipped and no introns were retained. the normal/expected splicing pattern was followed
- ``MULTI_RETAIN``: multiple introns were retained
- ``MULTI_SKIP``: multiple exons were skipped
- ``COMPLEX``: some combination of exon skipping and intron retention
"""

SPLICE_SITE_RADIUS = 2
""":class:`int`: number of bases away from an exon boundary considered to be part of the splice site such that if it were altered
        the splice site would be considered to be abrogated.
"""

PRIME = Vocab(FIVE=5, THREE=3)
""":class:`Vocab`: holds controlled vocabulary

- ``FIVE``: five prime
- ``THREE``: three prime
"""

START_AA = 'M'
""":class:`str`: The amino acid expected to start translation
"""
STOP_AA = '*'
""":class:`str`: The amino acid expected to end translation
"""

GIESMA_STAIN = Vocab(
    GNEG='gneg',
    GPOS50='gpos50',
    GPOS75='gpos75',
    GPOS25='gpos25',
    GPOS100='gpos100',
    ACEN='acen',
    GVAR='gvar',
    STALK='stalk'
)
""":class:`Vocab`: holds controlled vocabulary relating to stains of chromosome bands"""

# content related to tabbed files for input/output
# ensure that we don't have to change ALL the code when we update column names


COLUMNS = Vocab(
    library='library',
    cluster_id='cluster_id',
    cluster_size='cluster_size',
    validation_id='validation_id',
    annotation_id='annotation_id',
    product_id='product_id',
    event_type='event_type',
    pairing='pairing',
    gene1='gene1',
    gene1_direction='gene1_direction',
    gene2='gene2',
    gene2_direction='gene2_direction',
    gene_product_type='gene_product_type',
    transcript1='transcript1',
    transcript2='transcript2',
    fusion_splicing_pattern='fusion_splicing_pattern',
    fusion_cdna_coding_start='fusion_cdna_coding_start',
    fusion_cdna_coding_end='fusion_cdna_coding_end',
    fusion_mapped_domains='fusion_mapped_domains',
    fusion_sequence_fasta_id='fusion_sequence_fasta_id',
    fusion_sequence_fasta_file='fusion_sequence_fasta_file',
    annotation_figure='annotation_figure',
    annotation_figure_legend='annotation_figure_legend',
    genes_encompassed='genes_encompassed',
    genes_overlapping_break1='genes_overlapping_break1',
    genes_overlapping_break2='genes_overlapping_break2',
    genes_proximal_to_break1='genes_proximal_to_break1',
    genes_proximal_to_break2='genes_proximal_to_break2',
    break1_chromosome='break1_chromosome',
    break1_position_start='break1_position_start',
    break1_position_end='break1_position_end',
    break1_orientation='break1_orientation',
    break1_strand='break1_strand',
    break1_seq='break1_seq',
    break2_chromosome='break2_chromosome',
    break2_position_start='break2_position_start',
    break2_position_end='break2_position_end',
    break2_orientation='break2_orientation',
    break2_strand='break2_strand',
    break2_seq='break2_seq',
    opposing_strands='opposing_strands',
    stranded='stranded',
    protocol='protocol',
    tools='tools',
    break1_call_method='break1_call_method',
    break1_ewindow='break1_ewindow',
    break1_ewindow_count='break1_ewindow_count',
    break1_ewindow_practical_coverage='break1_ewindow_practical_coverage',
    break1_homologous_seq='break1_homologous_seq',
    break1_split_read_names='break1_split_read_names',
    break1_split_reads='break1_split_reads',
    break1_split_reads_forced='break1_split_reads_forced',
    break2_call_method='break2_call_method',
    break2_ewindow='break2_ewindow',
    break2_ewindow_count='break2_ewindow_count',
    break2_ewindow_practical_coverage='break2_ewindow_practical_coverage',
    break2_homologous_seq='break2_homologous_seq',
    break2_split_read_names='break2_split_read_names',
    break2_split_reads='break2_split_reads',
    break2_split_reads_forced='break2_split_reads_forced',
    contig_alignment_score='contig_alignment_score',
    contig_remap_score='contig_remap_score',
    contig_remapped_reads='contig_remapped_reads',
    contig_remapped_read_names='contig_remapped_read_names',
    contig_seq='contig_seq',
    contigs_aligned='contigs_aligned',
    contigs_assembled='contigs_assembled',
    flanking_median_fragment_size='flanking_median_fragment_size',
    flanking_pairs='flanking_pairs',
    flanking_pairs_compatible='flanking_pairs_compatible',
    flanking_pairs_read_names='flanking_pairs_read_names',
    flanking_pairs_compatible_read_names='flanking_pairs_compatible_read_names',
    flanking_stdev_fragment_size='flanking_stdev_fragment_size',
    linking_split_read_names='linking_split_read_names',
    linking_split_reads='linking_split_reads',
    raw_break1_half_mapped_reads='raw_break1_half_mapped_reads',
    raw_break1_split_reads='raw_break1_split_reads',
    raw_break2_half_mapped_reads='raw_break2_half_mapped_reads',
    raw_break2_split_reads='raw_break2_split_reads',
    raw_flanking_pairs='raw_flanking_pairs',
    raw_spanning_reads='raw_spanning_reads',
    untemplated_seq='untemplated_seq',
)
""":class:`Vocab`: Column names for i/o files used throughout the pipeline


.. glossary::
    :sorted:

    library
        Identifier for the library/source

    cluster_id
        Identifier for the merging/clustering step

    cluster_size
        The number of breakpoint pair calls that were grouped in creating the cluster

    validation_id
        Identifier for the validation step

    annotation_id
        Identifier for the annotation step

    product_id
        Unique identifier of the final fusion including splicing and ORF decision from the annotation step

    event_type
        :class:`SVTYPE` - The classification of the event

    pairing
        A semi colon delimited of event identifiers i.e. <annotation_id>_<splicing pattern>_<cds start>_<cds end>

    gene1
        Gene for the current annotation at the first breakpoint

    gene1_direction
        :class:`PRIME` - The direction/prime of the gene

    gene2
        Gene for the current annotation at the second breakpoint

    gene2_direction
        :class:`PRIME` - The direction/prime of the gene. Has the following possible values

    gene_product_type
        :class:`GENE_PRODUCT_TYPE` - Describes if the putative fusion product will be sense or anti-sense

    transcript1
        Transcript for the current annotation at the first breakpoint

    transcript2
        Transcript for the current annotation at the second breakpoint

    fusion_splicing_pattern
        :class:`SPLICE_TYPE` - Type of splicing pattern used to create the fusion cDNA.

    fusion_cdna_coding_start
        Position wrt the 5\ end of the fusion transcript where coding begins first base of the Met amino acid.

    fusion_cdna_coding_end
        Position wrt the 5\ end of the fusion transcript where coding ends last base of the stop codon

    fusion_mapped_domains
        ``JSON`` - List of domains in json format where each domain start and end positions are given wrt to the fusion
        transcript and the mapping quality is the number of matching amino acid positions over the total
        number of amino acids. The sequence is the amino acid sequence of the domain on the reference/original
        transcript

    fusion_sequence_fasta_id
        The sequence identifier for the cdna sequence output fasta file

    fusion_sequence_fasta_file
        Path to the corresponding fasta output file

    annotation_figure
        File path to the svg drawing representing the annotation

    annotation_figure_legend
        ``JSON`` - JSON data for the figure legend

    genes_encompassed
        Applies to intrachromosomal events only. List of genes which overlap any region that occurs between both
        breakpoints. For example in a deletion event these would be deleted genes.

    genes_overlapping_break1
        list of genes which overlap the first breakpoint

    genes_overlapping_break2
        list of genes which overlap the second breakpoint

    genes_proximal_to_break1
        list of genes near the breakpoint and the distance away from the breakpoint

    genes_proximal_to_break2
        list of genes near the breakpoint and the distance away from the breakpoint

    break1_chromosome
        :class:`str` - The name of the chromosome on which breakpoint 1 is situated

    break1_position_start
        :class:`int` - Start integer inclusive 1-based of the range representing breakpoint 1

    break1_position_end
        :class:`int` - End integer inclusive 1-based of the range representing breakpoint 1

    break1_orientation
        :class:`ORIENT` - The side of the breakpoint wrt the positive/forward strand that is retained.

    break1_strand
        :class:`STRAND` - The strand wrt to the reference positive/forward strand at this breakpoint.

    break1_seq
        :class:`str` - The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand

    break2_chromosome
        The name of the chromosome on which breakpoint 2 is situated

    break2_position_start
        :class:`int` - Start integer inclusive 1-based of the range representing breakpoint 2

    break2_position_end
        :class:`int` - End integer inclusive 1-based of the range representing breakpoint 2

    break2_orientation
        :class:`ORIENT` - The side of the breakpoint wrt the positive/forward strand that is retained.

    break2_strand
        :class:`STRAND` - The strand wrt to the reference positive/forward strand at this breakpoint.

    break2_seq
        :class:`str` - The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand

    opposing_strands
        :class:`bool` - Specifies if breakpoints are on opposite strands wrt to the reference. Expects a boolean

    stranded
        :class:`bool` - Specifies if the sequencing protocol was strand specific or not. Expects a boolean

    protocol
        :class:`PROTOCOL` - Specifies the type of library

    tools
        The tools that called the event originally from the cluster step. Should be a semi-colon delimited list of
        <tool name>_<tool version>

    contigs_assembled
        :class:`int` - Number of contigs that were built from split read sequences

    contigs_aligned
        :class:`int` - Number of contigs that were able to align

    contig_seq
        :class:`str` - Sequence of the current contig wrt to the positive forward strand if not strand specific

    contig_remap_score
        :class:`float` - Score representing the number of sequences from the set of sequences given to the assembly
        algorithm that were aligned to the resulting contig with an acceptable scoring based on user-set thresholds.
        For any sequence its contribution to the score is divided by the number of mappings to give less weight to
        multimaps

    contig_remapped_reads
        :class:`int` - the number of reads from the input bam that map to the assembled contig

    contig_remapped_read_names
        read query names for the reads that were remapped. A -1 or -2 has been appended to the end of the name to
        indicate if this is the first or second read in the pair

    contig_alignment_score
        :class:`float` - A rank based on the alignment tool blat etc. of the alignment being used. An average if
        split alignments were used. Lower numbers indicate a better alignment. If it was the best alignment possible
        then this would be zero.

    break1_call_method
        :class:`CALL_METHOD` - The method used to call the first breakpoint

    break2_call_method
        :class:`CALL_METHOD` - The method used to call the second breakpoint

    flanking_pairs
        :class:`int` - Number of read-pairs where one read aligns to the first breakpoint window and the second read
        aligns to the other. The count here is based on the number of unique query names

    flanking_pairs_compatible
        :class:`int` - Number of flanking pairs of a compatible orientation type. This applies to insertions and
        duplications. Flanking pairs supporting an insertion will be compatible to a duplication and flanking pairs
        supporting a duplication will be compatible to an insertion (possibly indicating an internal translocation)

    flanking_median_fragment_size
        :class:`int` - The median fragment size of the flanking reads being used as evidence

    flanking_stdev_fragment_size
        :class:`float` - The standard deviation in fragment size of the flanking reads being used as evidence

    break1_split_reads
        :class:`int` - Number of split reads that call the exact breakpoint given

    break1_split_reads_forced
        :class:`int` - Number of split reads which were aligned to the opposite breakpoint window using a targeted
        alignment

    break2_split_reads
        :class:`int` - Number of split reads that call the exact breakpoint given

    break2_split_reads_forced
        :class:`int` - Number of split reads which were aligned to the opposite breakpoint window using a targeted
        alignment

    linking_split_reads
        :class:`int` - Number of split reads that align to both breakpoints

    untemplated_seq
        :class:`str` - The untemplated/novel sequence between the breakpoints

    break1_homologous_seq
        :class:`str` - Sequence in common at the first breakpoint and other side of the second breakpoint

    break2_homologous_seq
        :class:`str` - Sequence in common at the second breakpoint and other side of the first breakpoint

    break1_ewindow
        Window where evidence was gathered for the first breakpoint

    break1_ewindow_count
        :class:`int` - Number of reads processed/looked-at in the first evidence window

    break1_ewindow_practical_coverage
        :class:`float` - break2_ewindow_practical_coverage, break1_ewindow_count / len(break1_ewindow). Not the actual
        coverage as bins are sampled within and there is a read limit cutoff

    break2_ewindow
        Window where evidence was gathered for the second breakpoint

    break2_ewindow_count
        :class:`int` - Number of reads processed/looked-at in the second evidence window

    break2_ewindow_practical_coverage
        :class:`float` - break2_ewindow_practical_coverage, break2_ewindow_count / len(break2_ewindow). Not the actual
        coverage as bins are sampled within and there is a read limit cutoff

    raw_flanking_pairs
        :class:`int` - Number of flanking reads before calling the breakpoint. The count here is based on the number of
        unique query names

    raw_spanning_reads
        :class:`int` - Number of spanning reads collected during evidence collection before calling the breakpoint

    raw_break1_split_reads
        :class:`int` - Number of split reads before calling the breakpoint

    raw_break2_split_reads
        :class:`int` - Number of split reads before calling the breakpoint

"""

def sort_columns(input_columns):
    order = {}
    for i, col in enumerate(COLUMNS.values()):
        order[col] = i

    temp = sorted([c for c in input_columns if c in order], key=lambda x: order[x])
    temp = temp + sorted([c for c in input_columns if c not in order])
    return temp
