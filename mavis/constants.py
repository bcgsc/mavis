"""
module responsible for small utility functions and constants used throughout the structural_variant package
"""
import argparse
import re

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from tab import cast_boolean


class MavisNamespace(argparse.Namespace):
    """
    Namespace to hold module constants
    """
    reserved_attr = ['_types', '_defns']

    def __init__(self, **kwargs):
        for k in kwargs:
            if k in MavisNamespace.reserved_attr:
                raise AttributeError('reserved attribute {} cannot be used'.format(k))
        self._defns = {}
        self._types = {}
        argparse.Namespace.__init__(self, **kwargs)
        for attr, value in kwargs.items():
            self._set_type(attr, type(value))

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, val):
        if key in MavisNamespace.reserved_attr:
            raise AttributeError('reserved attribute {} cannot be used'.format(key))
        self.__dict__[key] = val

    def flatten(self):
        items = {}
        items.update(self.items())
        return items

    def get(self, key, *pos):
        if len(pos) > 1:
            raise TypeError('too many arguments. get takes a single \'default\' value argument')
        try:
            return self[key]
        except AttributeError as err:
            if len(pos) == 1:
                return pos[0]
            raise err

    def keys(self):
        return [k for k in self.__dict__ if k not in MavisNamespace.reserved_attr]

    def values(self):
        return [self[k] for k in self.keys()]

    def enforce(self, value):
        if value not in self.values():
            raise KeyError('value {0} is not a valid member'.format(value), self.values())
        return value

    def reverse(self, value):
        """for a given value, return the associated key"""
        result = []
        for key in self.keys():
            if self[key] == value:
                result.append(key)
        if len(result) > 1:
            raise KeyError('could not reverse, the mapping is not unique', value, result)
        elif not result:
            raise KeyError('input value is not assigned to a key', value)
        return result[0]

    def __iter__(self):
        return iter(self.keys())

    def _set_type(self, attr, cast_type):
        if cast_type == bool:
            self._types[attr] = cast_boolean
        else:
            self._types[attr] = cast_type

    def type(self, attr):
        return self._types[attr]

    def define(self, attr, *pos):
        if len(pos) > 1:
            raise TypeError('too many arguments. define takes a single \'default\' value argument')
        try:
            return self._defns[attr]
        except KeyError as err:
            if len(pos) == 1:
                return pos[0]
            raise err

    def add(self, attr, *pos, **kwargs):
        """
        Add an attribute to the name space. Optionally include cast_type and definition
        """
        if len(pos) > 1:
            raise TypeError('add() takes 3 positional arguments but more were given')
        elif len(pos) == 1:
            if 'value' in kwargs:
                raise TypeError('add() got multiple values for argument \'value\'')
            kwargs['value'] = pos[0]
        value = kwargs.pop('value', attr)
        self[attr] = value
        if 'cast_type' not in kwargs:
            self._set_type(attr, type(value))
        else:
            self._types[attr] = kwargs.pop('cast_type')
        if 'defn' in kwargs:
            self._defns[attr] = kwargs.pop('defn')
        if kwargs:
            raise TypeError('invalid arguments: {}'.format(kwargs.keys()))


def float_fraction(num):
    """
    cast input to a float

    Args:
        num: input to cast

    Returns:
        float

    Raises:
        TypeError: if the input cannot be cast to a float or the number is not between 0 and 1
    """
    try:
        num = float(num)
    except ValueError:
        raise argparse.ArgumentTypeError('Must be a value between 0 and 1')
    if num < 0 or num > 1:
        raise argparse.ArgumentTypeError('Must be a value between 0 and 1')
    return num


COMPLETE_STAMP = 'MAVIS.COMPLETE'

PIPELINE_STEP = MavisNamespace(
    ANNOTATE='annotate',
    VALIDATE='validate',
    PIPELINE='pipeline',
    CLUSTER='cluster',
    PAIR='pairing',
    SUMMARY='summary',
    CHECKER='checker',
    CONFIG='config',
    CONVERT='convert'
)


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
    input_string = str(s)
    if not re.match('^[A-Za-z]*$', input_string):
        raise ValueError('unexpected sequence format. cannot reverse complement', input_string)
    input_string = Seq(input_string, DNA_ALPHABET)
    return str(input_string.reverse_complement())


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

ORIENT = MavisNamespace(LEFT='L', RIGHT='R', NS='?')
""":class:`MavisNamespace`: holds controlled vocabulary for allowed orientation values

- ``LEFT``: left wrt to the positive/forward strand
- ``RIGHT``: right wrt to the positive/forward strand
- ``NS``: orientation is not specified
"""
setattr(ORIENT, 'expand', lambda x: [ORIENT.LEFT, ORIENT.RIGHT] if x == ORIENT.NS else [x])
setattr(ORIENT, 'compare', lambda x, y: True if ORIENT.NS in [x, y] else (x == y))

PROTOCOL = MavisNamespace(GENOME='genome', TRANS='transcriptome')
""":class:`MavisNamespace`: holds controlled vocabulary for allowed protocol values

- ``GENOME``: genome
- ``TRANS``: transcriptome
"""

DISEASE_STATUS = MavisNamespace(DISEASED='diseased', NORMAL='normal')
""":class:`MavisNamespace`: holds controlled vocabulary for allowed disease status

- ``DISEASED``: diseased
- ``NORMAL``: normal
"""

STRAND = MavisNamespace(POS='+', NEG='-', NS='?')
""":class:`MavisNamespace`: holds controlled vocabulary for allowed strand values

- ``POS``: the positive/forward strand
- ``NEG``: the negative/reverse strand
- ``NS``: strand is not specified
"""
setattr(STRAND, 'expand', lambda x: [STRAND.POS, STRAND.NEG] if x == STRAND.NS else [x])
setattr(STRAND, 'compare', lambda x, y: True if STRAND.NS in [x, y] else (x == y))

SVTYPE = MavisNamespace(
    DEL='deletion',
    TRANS='translocation',
    ITRANS='inverted translocation',
    INV='inversion',
    INS='insertion',
    DUP='duplication'
)
""":class:`MavisNamespace`: holds controlled vocabulary for acceptable structural variant classifications

- ``DEL``: deletion
- ``TRANS``: translocation
- ``ITRANS``: inverted translocation
- ``INV``: inversion
- ``INS``: insertion
- ``DUP``: duplication
"""

CIGAR = MavisNamespace(M=0, I=1, D=2, N=3, S=4, H=5, P=6, X=8, EQ=7)
""":class:`MavisNamespace`: Enum-like. For readable cigar values

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
""":class:`int`: mapping quality value to indicate mapping was not performed/calculated"""

PYSAM_READ_FLAGS = MavisNamespace(
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

""":class:`MavisNamespace`: Enum-like. For readable PYSAM flag constants

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
    if not xset.intersection(yset):
        return False
    return True


DNA_ALPHABET = alphabet = Gapped(ambiguous_dna, '-')
DNA_ALPHABET.match = lambda x, y: _match_ambiguous_dna(x, y)

FLAGS = MavisNamespace(LQ='LOWQUAL')

READ_PAIR_TYPE = MavisNamespace(RR='RR', LL='LL', RL='RL', LR='LR')

CALL_METHOD = MavisNamespace(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads', SPAN='spanning reads')
""":class:`MavisNamespace`: holds controlled vocabulary for allowed call methods

- ``CONTIG``: a contig was assembled and aligned across the breakpoints
- ``SPLIT``: the event was called by :term:`split read`
- ``FLANK``: the event was called by :term:`flanking read pair`
- ``SPAN``: the event was called by :term:`spanning read`
"""

GENE_PRODUCT_TYPE = MavisNamespace(SENSE='sense', ANTI_SENSE='anti-sense')
""":class:`MavisNamespace`: controlled vocabulary for gene products

- ``SENSE``: the gene product is a sense fusion
- ``ANTI_SENSE``: the gene product is anti-sense
"""

PRIME = MavisNamespace(FIVE=5, THREE=3)
""":class:`MavisNamespace`: holds controlled vocabulary

- ``FIVE``: five prime
- ``THREE``: three prime
"""

START_AA = 'M'
""":class:`str`: The amino acid expected to start translation
"""
STOP_AA = '*'
""":class:`str`: The amino acid expected to end translation
"""

GIEMSA_STAIN = MavisNamespace(
    GNEG='gneg',
    GPOS50='gpos50',
    GPOS75='gpos75',
    GPOS25='gpos25',
    GPOS100='gpos100',
    ACEN='acen',
    GVAR='gvar',
    STALK='stalk'
)
""":class:`MavisNamespace`: holds controlled vocabulary relating to stains of chromosome bands"""

# content related to tabbed files for input/output
# ensure that we don't have to change ALL the code when we update column names


COLUMNS = MavisNamespace(
    library='library',
    cluster_id='cluster_id',
    cluster_size='cluster_size',
    validation_id='validation_id',
    annotation_id='annotation_id',
    product_id='product_id',
    event_type='event_type',
    pairing='pairing',
    inferred_pairing='inferred_pairing',
    gene1='gene1',
    gene1_direction='gene1_direction',
    gene2='gene2',
    gene2_direction='gene2_direction',
    gene1_aliases='gene1_aliases',
    gene2_aliases='gene2_aliases',
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
    exon_last_5prime='exon_last_5prime',
    exon_first_3prime='exon_first_3prime',
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
    disease_status='disease_status',
    tools='tools',
    call_method='call_method',
    break1_ewindow='break1_ewindow',
    break1_ewindow_count='break1_ewindow_count',
    break1_ewindow_practical_coverage='break1_ewindow_practical_coverage',
    break1_homologous_seq='break1_homologous_seq',
    break1_split_read_names='break1_split_read_names',
    break1_split_reads='break1_split_reads',
    break1_split_reads_forced='break1_split_reads_forced',
    break2_ewindow='break2_ewindow',
    break2_ewindow_count='break2_ewindow_count',
    break2_ewindow_practical_coverage='break2_ewindow_practical_coverage',
    break2_homologous_seq='break2_homologous_seq',
    break2_split_read_names='break2_split_read_names',
    break2_split_reads='break2_split_reads',
    break2_split_reads_forced='break2_split_reads_forced',
    contig_alignment_query_consumption='contig_alignment_query_consumption',
    contig_alignment_score='contig_alignment_score',
    contig_alignment_query_name='contig_alignment_query_name',
    contig_read_depth='contig_read_depth',
    contig_break1_read_depth='contig_break1_read_depth',
    contig_break2_read_depth='contig_break2_read_depth',
    contig_blat_rank='contig_blat_rank',
    contig_build_score='contig_build_score',
    contig_remap_score='contig_remap_score',
    contig_remap_coverage='contig_remap_coverage',
    contig_remapped_read_names='contig_remapped_read_names',
    contig_remapped_reads='contig_remapped_reads',
    contig_seq='contig_seq',
    contig_strand_specific='contig_strand_specific',
    contigs_aligned='contigs_aligned',
    contigs_assembled='contigs_assembled',
    spanning_reads='spanning_reads',
    spanning_read_names='spanning_read_names',
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
    filter_comment='filter_comment',
    cdna_synon='cdna_synon',
    protein_synon='protein_synon'
)
""":class:`MavisNamespace`: Column names for i/o files used throughout the pipeline



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

    inferred_pairing
        A semi colon delimited of event identifiers i.e. <annotation_id>_<splicing pattern>_<cds start>_<cds end>
        which were paired to the current event based on predicted products

    pairing
        A semi colon delimited of event identifiers i.e. <annotation_id>_<splicing pattern>_<cds start>_<cds end>
        which were paired to the current event based on breakpoint positions

    gene1
        Gene for the current annotation at the first breakpoint

    gene1_direction
        :class:`PRIME` - The direction/prime of the gene

    gene2
        Gene for the current annotation at the second breakpoint

    gene2_direction
        :class:`PRIME` - The direction/prime of the gene. Has the following possible values

    gene1_aliases
        Other gene names associated with the current annotation at the first breakpoint

    gene2_aliases
        Other gene names associated with the current annotation at the second breakpoint

    gene_product_type
        :class:`GENE_PRODUCT_TYPE` - Describes if the putative fusion product will be sense or anti-sense

    fusion_cdna_coding_end
        Position wrt the 5' end of the fusion transcript where coding ends last base of the stop codon

    transcript1
        Transcript for the current annotation at the first breakpoint

    transcript2
        Transcript for the current annotation at the second breakpoint

    fusion_splicing_pattern
        :class:`SPLICE_TYPE` - Type of splicing pattern used to create the fusion cDNA.

    fusion_cdna_coding_start
        Position wrt the 5' end of the fusion transcript where coding begins first base of the Met amino acid.

    fusion_cdna_coding_end
        Position wrt the 5' end of the fusion transcript where coding ends last base of the stop codon

    fusion_mapped_domains
        ``JSON`` - List of domains in :term:`JSON` format where each domain start and end positions are given wrt to the fusion
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
        ``JSON`` - :term:`JSON` data for the figure legend

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

    contig_alignment_query_name
        The query name for the contig alignment. Should match the 'read' name(s) in the .contigs.bam output file

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

    contig_alignment_reference_start
        The reference start(s) <chr>:<position> of the contig alignment. Semi-colon delimited

    contig_alignment_cigar
        The cigar string(s) representing the contig alignment. Semi-colon delimited

    contig_remap_coverage
        :class:`float` - Fraction of the contig sequence which is covered by the remapped reads

    contig_build_score
        :class:`int` - Score representing the edge weights of all edges used in building the sequence

    contig_strand_specific
        :class:`bool` - A flag to indicate if it was possible to resolve the strand for this contig

    spanning_reads
        :class:`int` - the number of spanning reads which support the event

    spanning_read_names
        read query names of the spanning reads which support the current event

    call_method
        :class:`CALL_METHOD` - The method used to call the breakpoints

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

    cdna_synon
        semi-colon delimited list of transcript ids which have an identical cdna sequence to the cdna sequence of the
        current fusion product

    protein_synon
        semi-colon delimited list of transcript ids which produce a translation with an identical amino-acid sequence
        to the current fusion product
"""


def sort_columns(input_columns):
    order = {}
    for i, col in enumerate(COLUMNS.values()):
        order[col] = i
    temp = sorted([c for c in input_columns if c in order], key=lambda x: order[x])
    temp = temp + sorted([c for c in input_columns if c not in order])
    return temp
