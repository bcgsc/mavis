"""
module responsible for small utility functions and constants used throughout the structural_variant package
"""
import argparse
import re

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from tab import cast_boolean, cast_null


class MavisNamespace(argparse.Namespace):
    """
    Namespace to hold module constants

    Example:
        >>> nspace = MavisNamespace(thing=1, otherthing=2)
        >>> nspace.thing
        1
        >>> nspace.otherthing
        2
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
            if not attr.startswith('_'):
                self._set_type(attr, type(value))

    def items(self):
        """
        Example:
            >>> MavisNamespace(thing=1, otherthing=2).items()
            [('thing', 1), ('otherthing', 2)]
        """
        return [(k, self[k]) for k in self.keys()]

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, val):
        if key in MavisNamespace.reserved_attr:
            raise AttributeError('reserved attribute {} cannot be used'.format(key))
        self.__dict__[key] = val

    def flatten(self):
        """
        returns the namespace (minus types and definitions) as a dictionary

        Example:
            >>> MavisNamespace(thing=1, otherthing=2).flatten()
            {'thing': 1, 'otherthing': 2}
        """
        items = {}
        items.update(self.items())
        return items

    def get(self, key, *pos):
        """
        get an attribute, return a default (if given) if the attribute does not exist

        Example:
            >>> nspace = MavisNamespace(thing=1, otherthing=2)
            >>> nspace.get('thing', 2)
            1
            >>> nspace.get('nonexistant_thing', 2)
            2
            >>> nspace.get('nonexistant_thing')
            Traceback (most recent call last):
            ....
        """
        if len(pos) > 1:
            raise TypeError('too many arguments. get takes a single \'default\' value argument')
        try:
            return self[key]
        except AttributeError as err:
            if pos:
                return pos[0]
            raise err

    def keys(self):
        """
        get the attribute keys as a list

        Example:
            >>> MavisNamespace(thing=1, otherthing=2).keys()
            ['thing', 'otherthing']
        """
        return [k for k in self.__dict__ if not k.startswith('_')]

    def values(self):
        """
        get the attribute values as a list

        Example:
            >>> MavisNamespace(thing=1, otherthing=2).values()
            [1, 2]
        """
        return [self[k] for k in self.keys()]

    def enforce(self, value):
        """
        checks that the current namespace has a given value

        Returns:
            the input value

        Raises:
            KeyError: the value did not exist

        Example:
            >>> nspace = MavisNamespace(thing=1, otherthing=2)
            >>> nspace.enforce(1)
            1
            >>> nspace.enforce(3)
            Traceback (most recent call last):
            ....
        """
        if value not in self.values():
            raise KeyError('value {0} is not a valid member'.format(value), self.values())
        return value

    def reverse(self, value):
        """
        for a given value, return the associated key

        Args:
            value: the value to get the key/attribute name for

        Raises:
            KeyError: the value is not unique
            KeyError: the value is not assigned

        Example:
            >>> nspace = MavisNamespace(thing=1, otherthing=2)
            >>> nspace.reverse(1)
            'thing'
        """
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
        """
        returns the type

        Example:
            >>> nspace = MavisNamespace(thing=1, otherthing=2)
            >>> nspace.type('thing')
            <class 'int'>
        """
        return self._types[attr]

    def define(self, attr, *pos):
        """
        Get the definition of a given attribute or return a default (when given) if the attribute does not exist

        Returns:
            str: definition for the attribute

        Raises:
            KeyError: the attribute does not exist and a default was not given

        Example:
            >>> nspace = MavisNamespace()
            >>> nspace.add('thing', 1, defn='I am a thing')
            >>> nspace.add('otherthing', 2)
            >>> nspace.define('thing')
            'I am a thing'
            >>> nspace.define('otherthing')
            Traceback (most recent call last):
            ....
            >>> nspace.define('otherthing', 'I am some other thing')
            'I am some other thing'
        """
        if len(pos) > 1:
            raise TypeError('too many arguments. define takes a single \'default\' value argument')
        try:
            return self._defns[attr]
        except KeyError as err:
            if pos:
                return pos[0]
            raise err

    def add(self, attr, *pos, **kwargs):
        """
        Add an attribute to the name space. Optionally include cast_type and definition

        Example:
            >>> nspace = MavisNamespace()
            >>> nspace.add('thing', 1, int, 'I am a thing')
            >>> nspace = MavisNamespace()
            >>> nspace.add('thing', 1, int)
            >>> nspace = MavisNamespace()
            >>> nspace.add('thing', 1)
            >>> nspace = MavisNamespace()
            >>> nspace.add('thing', value=1, cast_type=int, defn='I am a thing')
        """
        if len(pos) > 1:
            raise TypeError('add() takes 3 positional arguments but more were given')

        for value, argname in zip(pos, ['value', 'cast_type', 'defn'][:len(pos)]):
            if argname in kwargs:
                raise TypeError('add() got multiple values for argument {}'.format(repr(argname)))
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

    def __call__(self, value):
        try:
            return self.enforce(value)
        except KeyError:
            raise TypeError('Invalid value {}. Cannot cast to type {}. Must be a valid member: {}'.format(
                repr(value), self.__class__.__name__, self.values()))


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


def nullable_int(num):
    """
    casts input to an int if not an accepted null value. See :func:tab.tab.cast_null
    """
    try:
        return cast_null(num)
    except TypeError:
        pass
    return int(num)


COMPLETE_STAMP = 'MAVIS.COMPLETE'

SUBCOMMAND = MavisNamespace(
    ANNOTATE='annotate',
    VALIDATE='validate',
    PIPELINE='pipeline',
    CLUSTER='cluster',
    PAIR='pairing',
    SUMMARY='summary',
    CHECKER='checker',
    CONFIG='config',
    CONVERT='convert',
    OVERLAY='overlay'
)
""":class:`MavisNamespace`: holds controlled vocabulary for allowed pipeline stage values

- annotate
- validate
- pipeline
- cluster
- pairing
- summary
- checker
- config
- convert
"""


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

CALL_METHOD = MavisNamespace(CONTIG='contig', SPLIT='split reads', FLANK='flanking reads', SPAN='spanning reads', INPUT='input')
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
    tracking_id='tracking_id',
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
    fusion_protein_hgvs='fusion_protein_hgvs',
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
    protein_synon='protein_synon',
    supplementary_call='supplementary_call',
    net_size='net_size',
    repeat_count='repeat_count',
    assumed_untemplated='assumed_untemplated'
)
""":class:`MavisNamespace`: Column names for i/o files used throughout the pipeline

- :term:`annotation_figure_legend`
- :term:`annotation_figure`
- :term:`annotation_id`
- :term:`break1_chromosome`
- :term:`break1_ewindow_count`
- :term:`break1_ewindow_practical_coverage`
- :term:`break1_ewindow`
- :term:`break1_homologous_seq`
- :term:`break1_orientation`
- :term:`break1_position_end`
- :term:`break1_position_start`
- :term:`break1_seq`
- :term:`break1_split_reads_forced`
- :term:`break1_split_reads`
- :term:`break1_strand`
- :term:`break2_chromosome`
- :term:`break2_ewindow_count`
- :term:`break2_ewindow_practical_coverage`
- :term:`break2_ewindow`
- :term:`break2_homologous_seq`
- :term:`break2_orientation`
- :term:`break2_position_end`
- :term:`break2_position_start`
- :term:`break2_seq`
- :term:`break2_split_reads_forced`
- :term:`break2_split_reads`
- :term:`break2_strand`
- :term:`call_method`
- :term:`cdna_synon`
- :term:`cluster_id`
- :term:`cluster_size`
- :term:`contig_alignment_cigar`
- :term:`contig_alignment_query_name`
- :term:`contig_alignment_reference_start`
- :term:`contig_alignment_score`
- :term:`contig_build_score`
- :term:`contig_remap_coverage`
- :term:`contig_remap_score`
- :term:`contig_remapped_read_names`
- :term:`contig_remapped_reads`
- :term:`contig_seq`
- :term:`contig_strand_specific`
- :term:`contigs_aligned`
- :term:`contigs_assembled`
- :term:`event_type`
- :term:`flanking_median_fragment_size`
- :term:`flanking_pairs_compatible`
- :term:`flanking_pairs`
- :term:`flanking_stdev_fragment_size`
- :term:`fusion_cdna_coding_end`
- :term:`fusion_cdna_coding_end`
- :term:`fusion_cdna_coding_start`
- :term:`fusion_mapped_domains`
- :term:`fusion_sequence_fasta_file`
- :term:`fusion_sequence_fasta_id`
- :term:`fusion_splicing_pattern`
- :term:`fusion_protein_hgvs`
- :term:`gene1_aliases`
- :term:`gene1_direction`
- :term:`gene1`
- :term:`gene2_aliases`
- :term:`gene2_direction`
- :term:`gene2`
- :term:`gene_product_type`
- :term:`genes_encompassed`
- :term:`genes_overlapping_break1`
- :term:`genes_overlapping_break2`
- :term:`genes_proximal_to_break1`
- :term:`genes_proximal_to_break2`
- :term:`inferred_pairing`
- :term:`library`
- :term:`linking_split_reads`
- :term:`net_size`
- :term:`opposing_strands`
- :term:`pairing`
- :term:`product_id`
- :term:`protein_synon`
- :term:`protocol`
- :term:`raw_break1_split_reads`
- :term:`raw_break2_split_reads`
- :term:`raw_flanking_pairs`
- :term:`raw_spanning_reads`
- :term:`spanning_read_names`
- :term:`spanning_reads`
- :term:`stranded`
- :term:`tools`
- :term:`tracking_id`
- :term:`transcript1`
- :term:`transcript2`
- :term:`supplementary_call`
- :term:`untemplated_seq`
- :term:`validation_id`
- :term:`repeat_count`
- :term:`assumed_untemplated`
"""


def sort_columns(input_columns):
    order = {}
    for i, col in enumerate(COLUMNS.values()):
        order[col] = i
    temp = sorted([c for c in input_columns if c in order], key=lambda x: order[x])
    temp = temp + sorted([c for c in input_columns if c not in order])
    return temp
