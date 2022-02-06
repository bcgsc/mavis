"""
module responsible for small utility functions and constants used throughout the structural_variant package
"""
import argparse
import re
from typing import List

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from mavis_config.constants import MavisNamespace

PROGNAME: str = 'mavis'
EXIT_OK: int = 0
EXIT_ERROR: int = 1
EXIT_INCOMPLETE: int = 2


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


class SPLICE_TYPE(MavisNamespace):
    """
    holds controlled vocabulary for allowed splice type classification values

    Attributes:
        RETAIN: an intron was retained
        SKIP: an exon was skipped
        NORMAL: no exons were skipped and no introns were retained. the normal/expected splicing pattern was followed
        MULTI_RETAIN: multiple introns were retained
        MULTI_SKIP: multiple exons were skipped
        COMPLEX: some combination of exon skipping and intron retention
    """

    RETAIN: str = 'retained intron'
    SKIP: str = 'skipped exon'
    NORMAL: str = 'normal'
    MULTI_RETAIN: str = 'retained multiple introns'
    MULTI_SKIP: str = 'skipped multiple exons'
    COMPLEX: str = 'complex'


COMPLETE_STAMP: str = 'MAVIS.COMPLETE'
"""Filename for all complete stamp files"""


CODON_SIZE: int = 3
"""the number of bases making up a codon"""


def reverse_complement(s: str) -> str:
    """
    wrapper for the Bio.Seq reverse_complement method

    Args:
        s: the input DNA sequence

    Returns:
        str: the reverse complement of the input sequence

    Warning:
        assumes the input is a DNA sequence

    Example:
        >>> reverse_complement('ATCCGGT')
        'ACCGGAT'
    """
    input_string = str(s)
    if not re.match('^[A-Za-z]*$', input_string):
        raise ValueError('unexpected sequence format. cannot reverse complement', input_string)
    seq = Seq(input_string, DNA_ALPHABET)
    return str(seq.reverse_complement())


def translate(s: str, reading_frame: int = 0) -> str:
    """
    given a DNA sequence, translates it and returns the protein amino acid sequence

    Args:
        s: the input DNA sequence
        reading_frame: where to start translating the sequence

    Returns:
        the amino acid sequence
    """
    reading_frame = reading_frame % CODON_SIZE

    temp = s[reading_frame:]
    if len(temp) % 3 == 1:
        temp = temp[:-1]
    elif len(temp) % 3 == 2:
        temp = temp[:-2]
    temp = Seq(temp, DNA_ALPHABET)
    return str(temp.translate())  # type: ignore


GAP: str = '-'


class ORIENT(MavisNamespace):
    """
    holds controlled vocabulary for allowed orientation values

    Attributes:
        LEFT: left wrt to the positive/forward strand
        RIGHT: right wrt to the positive/forward strand
        NS: orientation is not specified
    """

    LEFT: str = 'L'
    RIGHT: str = 'R'
    NS: str = '?'

    @classmethod
    def expand(cls, orientation) -> List[str]:
        if orientation == cls.NS:
            return [cls.LEFT, cls.RIGHT]
        return [orientation]

    @classmethod
    def compare(cls, first, second) -> bool:
        if cls.NS in {first, second}:
            return True
        return first == second


class PROTOCOL(MavisNamespace):
    """
    holds controlled vocabulary for allowed protocol values
    """

    GENOME: str = 'genome'
    TRANS: str = 'transcriptome'


class DISEASE_STATUS(MavisNamespace):
    """
    holds controlled vocabulary for allowed disease status
    """

    DISEASED: str = 'diseased'
    NORMAL: str = 'normal'


class STRAND(MavisNamespace):
    """
    holds controlled vocabulary for allowed strand values

    Attributes:
        POS: the positive/forward strand
        NEG: the negative/reverse strand
        NS: strand is not specified
    """

    POS: str = '+'
    NEG: str = '-'
    NS: str = '?'

    @classmethod
    def expand(cls, strand: str) -> List[str]:
        if strand == cls.NS:
            return [cls.POS, cls.NEG]
        return [strand]

    @classmethod
    def compare(cls, first, second) -> bool:
        if cls.NS in {first, second}:
            return True
        return first == second


class SVTYPE(MavisNamespace):
    """
    holds controlled vocabulary for acceptable structural variant classifications
    """

    DEL = 'deletion'
    TRANS = 'translocation'
    ITRANS: str = 'inverted translocation'
    INV: str = 'inversion'
    INS: str = 'insertion'
    DUP: str = 'duplication'


class CIGAR(MavisNamespace):
    """
    Enum-like. For readable cigar values


    Attributes:
        M: alignment match (can be a sequence match or mismatch)
        I: insertion to the reference
        D: deletion from the reference
        N: skipped region from the reference
        S: soft clipping (clipped sequences present in SEQ)
        H: hard clipping (clipped sequences NOT present in SEQ)
        P: padding (silent deletion from padded reference)
        EQ: sequence match (=)
        X: sequence mismatch

    Note:
        descriptions are taken from the `samfile documentation <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
    """

    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    X = 8
    EQ = 7


NA_MAPPING_QUALITY: int = 255
"""mapping quality value to indicate mapping was not performed/calculated"""


class PYSAM_READ_FLAGS(MavisNamespace):
    """
    Enum-like. For readable PYSAM flag constants

    Attributes:
        MULTIMAP: template having multiple segments in sequencing
        UNMAPPED: segment unmapped
        MATE_UNMAPPED: next segment in the template unmapped
        REVERSE: SEQ being reverse complemented
        MATE_REVERSE: SEQ of the next segment in the template being reverse complemented
        FIRST_IN_PAIR: the first segment in the template
        LAST_IN_PAIR: the last segment in the template
        SECONDARY: secondary alignment
        SUPPLEMENTARY: supplementary alignment

    Note:
        descriptions are taken from the `samfile documentation <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
    """

    REVERSE: int = 16
    MATE_REVERSE: int = 32
    UNMAPPED: int = 4
    MATE_UNMAPPED: int = 8
    FIRST_IN_PAIR: int = 64
    LAST_IN_PAIR: int = 128
    SECONDARY: int = 256
    MULTIMAP: int = 1
    SUPPLEMENTARY: int = 2048
    TARGETED_ALIGNMENT: str = 'ta'
    RECOMPUTED_CIGAR: str = 'rc'
    BLAT_RANK: str = 'br'
    BLAT_SCORE: str = 'bs'
    BLAT_ALIGNMENTS: str = 'ba'
    BLAT_PERCENT_IDENTITY: str = 'bi'
    BLAT_PMS: str = 'bp'


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


class FLAGS(MavisNamespace):
    LQ: str = 'LOWQUAL'


class READ_PAIR_TYPE(MavisNamespace):
    RR: str = 'RR'
    LL: str = 'LL'
    RL: str = 'RL'
    LR: str = 'LR'


class CALL_METHOD(MavisNamespace):
    """
    holds controlled vocabulary for allowed call methods

    Attributes:
        CONTIG: a contig was assembled and aligned across the breakpoints
        SPLIT: the event was called by [split read](/glossary/#split-read)
        FLANK: the event was called by [flanking read pair](/glossary/#flanking-read-pair)
        SPAN: the event was called by [spanning read](/glossary/#spanning-read)"""

    CONTIG: str = 'contig'
    SPLIT: str = 'split reads'
    FLANK: str = 'flanking reads'
    SPAN: str = 'spanning reads'
    INPUT: str = 'input'


class GENE_PRODUCT_TYPE(MavisNamespace):
    """
    controlled vocabulary for gene products

    Attributes:
        SENSE: the gene product is a sense fusion
        ANTI_SENSE: the gene product is anti-sense
    """

    SENSE: str = 'sense'
    ANTI_SENSE: str = 'anti-sense'


class PRIME(MavisNamespace):
    """
    Attributes:
        FIVE: five prime
        THREE: three prime
    """

    FIVE: int = 5
    THREE: int = 3


START_AA: str = 'M'
"""The amino acid expected to start translation
"""
STOP_AA: str = '*'
"""The amino acid expected to end translation
"""


class GIEMSA_STAIN(MavisNamespace):
    """
    holds controlled vocabulary relating to stains of chromosome bands
    """

    GNEG: str = 'gneg'
    GPOS33: str = 'gpos33'
    GPOS50: str = 'gpos50'
    GPOS66: str = 'gpos66'
    GPOS75: str = 'gpos75'
    GPOS25: str = 'gpos25'
    GPOS100: str = 'gpos100'
    ACEN: str = 'acen'
    GVAR: str = 'gvar'
    STALK: str = 'stalk'


# content related to tabbed files for input/output
# ensure that we don't have to change ALL the code when we update column names
class COLUMNS(MavisNamespace):
    """
    Column names for i/o files used throughout the pipeline

    see [column descriptions](/outputs/columns)
    """

    tracking_id: str = 'tracking_id'
    library: str = 'library'
    cluster_id: str = 'cluster_id'
    cluster_size: str = 'cluster_size'
    validation_id: str = 'validation_id'
    annotation_id: str = 'annotation_id'
    product_id: str = 'product_id'
    event_type: str = 'event_type'
    pairing: str = 'pairing'
    inferred_pairing: str = 'inferred_pairing'
    gene1: str = 'gene1'
    gene1_direction: str = 'gene1_direction'
    gene2: str = 'gene2'
    gene2_direction: str = 'gene2_direction'
    gene1_aliases: str = 'gene1_aliases'
    gene2_aliases: str = 'gene2_aliases'
    gene_product_type: str = 'gene_product_type'
    transcript1: str = 'transcript1'
    transcript2: str = 'transcript2'
    fusion_splicing_pattern: str = 'fusion_splicing_pattern'
    fusion_cdna_coding_start: str = 'fusion_cdna_coding_start'
    fusion_cdna_coding_end: str = 'fusion_cdna_coding_end'
    fusion_mapped_domains: str = 'fusion_mapped_domains'
    fusion_sequence_fasta_id: str = 'fusion_sequence_fasta_id'
    fusion_sequence_fasta_file: str = 'fusion_sequence_fasta_file'
    fusion_protein_hgvs: str = 'fusion_protein_hgvs'
    annotation_figure: str = 'annotation_figure'
    annotation_figure_legend: str = 'annotation_figure_legend'
    genes_encompassed: str = 'genes_encompassed'
    genes_overlapping_break1: str = 'genes_overlapping_break1'
    genes_overlapping_break2: str = 'genes_overlapping_break2'
    genes_proximal_to_break1: str = 'genes_proximal_to_break1'
    genes_proximal_to_break2: str = 'genes_proximal_to_break2'
    break1_chromosome: str = 'break1_chromosome'
    break1_position_start: str = 'break1_position_start'
    break1_position_end: str = 'break1_position_end'
    break1_orientation: str = 'break1_orientation'
    exon_last_5prime: str = 'exon_last_5prime'
    exon_first_3prime: str = 'exon_first_3prime'
    break1_strand: str = 'break1_strand'
    break1_seq: str = 'break1_seq'
    break2_chromosome: str = 'break2_chromosome'
    break2_position_start: str = 'break2_position_start'
    break2_position_end: str = 'break2_position_end'
    break2_orientation: str = 'break2_orientation'
    break2_strand: str = 'break2_strand'
    break2_seq: str = 'break2_seq'
    opposing_strands: str = 'opposing_strands'
    stranded: str = 'stranded'
    protocol: str = 'protocol'
    disease_status: str = 'disease_status'
    tools: str = 'tools'
    call_method: str = 'call_method'
    break1_ewindow: str = 'break1_ewindow'
    break1_ewindow_count: str = 'break1_ewindow_count'
    break1_homologous_seq: str = 'break1_homologous_seq'
    break1_split_read_names: str = 'break1_split_read_names'
    break1_split_reads: str = 'break1_split_reads'
    break1_split_reads_forced: str = 'break1_split_reads_forced'
    break2_ewindow: str = 'break2_ewindow'
    break2_ewindow_count: str = 'break2_ewindow_count'
    break2_homologous_seq: str = 'break2_homologous_seq'
    break2_split_read_names: str = 'break2_split_read_names'
    break2_split_reads: str = 'break2_split_reads'
    break2_split_reads_forced: str = 'break2_split_reads_forced'
    contig_alignment_query_consumption: str = 'contig_alignment_query_consumption'
    contig_alignment_score: str = 'contig_alignment_score'
    contig_alignment_query_name: str = 'contig_alignment_query_name'
    contig_read_depth: str = 'contig_read_depth'
    contig_break1_read_depth: str = 'contig_break1_read_depth'
    contig_break2_read_depth: str = 'contig_break2_read_depth'
    contig_alignment_rank: str = 'contig_alignment_rank'
    contig_build_score: str = 'contig_build_score'
    contig_remap_score: str = 'contig_remap_score'
    contig_remap_coverage: str = 'contig_remap_coverage'
    contig_remapped_read_names: str = 'contig_remapped_read_names'
    contig_remapped_reads: str = 'contig_remapped_reads'
    contig_seq: str = 'contig_seq'
    contig_strand_specific: str = 'contig_strand_specific'
    contigs_assembled: str = 'contigs_assembled'
    call_sequence_complexity: str = 'call_sequence_complexity'
    spanning_reads: str = 'spanning_reads'
    spanning_read_names: str = 'spanning_read_names'
    flanking_median_fragment_size: str = 'flanking_median_fragment_size'
    flanking_pairs: str = 'flanking_pairs'
    flanking_pairs_compatible: str = 'flanking_pairs_compatible'
    flanking_pairs_read_names: str = 'flanking_pairs_read_names'
    flanking_pairs_compatible_read_names: str = 'flanking_pairs_compatible_read_names'
    flanking_stdev_fragment_size: str = 'flanking_stdev_fragment_size'
    linking_split_read_names: str = 'linking_split_read_names'
    linking_split_reads: str = 'linking_split_reads'
    raw_break1_half_mapped_reads: str = 'raw_break1_half_mapped_reads'
    raw_break1_split_reads: str = 'raw_break1_split_reads'
    raw_break2_half_mapped_reads: str = 'raw_break2_half_mapped_reads'
    raw_break2_split_reads: str = 'raw_break2_split_reads'
    raw_flanking_pairs: str = 'raw_flanking_pairs'
    raw_spanning_reads: str = 'raw_spanning_reads'
    untemplated_seq: str = 'untemplated_seq'
    filter_comment: str = 'filter_comment'
    cdna_synon: str = 'cdna_synon'
    protein_synon: str = 'protein_synon'
    supplementary_call: str = 'supplementary_call'
    net_size: str = 'net_size'
    repeat_count: str = 'repeat_count'
    assumed_untemplated: str = 'assumed_untemplated'


def sort_columns(input_columns):
    order = {}
    for i, col in enumerate(COLUMNS.values()):
        order[col] = i
    temp = sorted([c for c in input_columns if c in order], key=lambda x: order[x])
    temp = temp + sorted([c for c in input_columns if c not in order])
    return temp


INTEGER_COLUMNS = {
    COLUMNS.break1_position_end,
    COLUMNS.break1_position_start,
    COLUMNS.break2_position_end,
    COLUMNS.break2_position_start,
}

FLOAT_COLUMNS = {
    COLUMNS.break1_ewindow_count,
    COLUMNS.break1_split_reads_forced,
    COLUMNS.break1_split_reads,
    COLUMNS.break2_ewindow_count,
    COLUMNS.break2_split_reads_forced,
    COLUMNS.break2_split_reads,
    COLUMNS.cluster_size,
    COLUMNS.contig_alignment_query_consumption,
    COLUMNS.contig_alignment_rank,
    COLUMNS.contig_alignment_score,
    COLUMNS.contig_break1_read_depth,
    COLUMNS.contig_break2_read_depth,
    COLUMNS.contig_build_score,
    COLUMNS.contig_read_depth,
    COLUMNS.contig_remap_score,
    COLUMNS.contig_remapped_reads,
    COLUMNS.contigs_assembled,
    COLUMNS.flanking_pairs_compatible,
    COLUMNS.flanking_pairs,
    COLUMNS.linking_split_reads,
    COLUMNS.raw_break1_half_mapped_reads,
    COLUMNS.raw_break1_split_reads,
    COLUMNS.raw_break2_half_mapped_reads,
    COLUMNS.raw_break2_split_reads,
    COLUMNS.raw_flanking_pairs,
    COLUMNS.raw_spanning_reads,
    COLUMNS.repeat_count,
    COLUMNS.spanning_reads,
}

BOOLEAN_COLUMNS = {COLUMNS.opposing_strands, COLUMNS.stranded, COLUMNS.supplementary_call}

SUMMARY_LIST_COLUMNS = {
    COLUMNS.annotation_figure,
    COLUMNS.annotation_id,
    COLUMNS.break1_split_reads,
    COLUMNS.break2_split_reads,
    COLUMNS.call_method,
    COLUMNS.contig_alignment_score,
    COLUMNS.contig_remapped_reads,
    COLUMNS.contig_seq,
    COLUMNS.event_type,
    COLUMNS.flanking_pairs,
    COLUMNS.pairing,
    COLUMNS.product_id,
    COLUMNS.spanning_reads,
    COLUMNS.tools,
    COLUMNS.tools,
    COLUMNS.tracking_id,
}
