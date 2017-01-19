import TSV
import re
from .genomic import Gene, Transcript, usTranscript, Exon
from .base import BioInterval
from .protein import Domain, Translation
from ..interval import Interval
from ..constants import STRAND
from Bio import SeqIO


def load_masking_regions(filepath):
    """
    reads a file of regions. The expect input format for the file is tab-delimited and

    +---------------+---------------+-----------------------+
    | column name   | value type    | description           |
    +===============+===============+=======================+
    | chr           | string        | the chromosome name   |
    +---------------+---------------+-----------------------+
    | start         | int           | the start position    |
    +---------------+---------------+-----------------------+
    | end           | int           | the end position      |
    +---------------+---------------+-----------------------+
    | name          | string        | label for the region  |
    +---------------+---------------+-----------------------+

    Args:
        filepath (str): path to the input tab-delimited file
    Returns:
        :class:`dict` of :class:`str` and :class:`list` of :class:`BioInterval`:
            a dictionary keyed by chromosome name with values of lists of regions on the chromosome

    Example:
        >>> m = load_masking_regions('filename')
        >>> m['1']
        [BioInterval(), BioInterval(), ...]
    """
    header, rows = TSV.read_file(
        filepath,
        require=['chr', 'start', 'end', 'name'],
        cast={'start': int, 'end': int, 'chr': lambda x: re.sub('^chr', '', x)}
    )
    regions = {}
    for row in rows:
        r = BioInterval(reference_object=row['chr'], start=row['start'], end=row['end'], name=row['name'])
        regions.setdefault(r.reference_object, []).append(r)
    return regions


def load_reference_genes(filepath, verbose=True):
    """
    given a file in the std input format (see below) reads and return a list of genes (and sub-objects)

    +-----------------------+---------------------------+-----------------------------------------------------------+
    | column name           | example                   | description                                               |
    +=======================+===========================+===========================================================+
    | ensembl_transcript_id | ENST000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | ensembl_gene_id       | ENSG000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | strand                | -1                        | positive or negative 1                                    |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_start     | 44                        | where translation begins relative to the start of the cdna|
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_end       | 150                       | where translation terminates                              |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | genomic_exon_ranges   | 100-201;334-412;779-830   | semi-colon demitited exon start/ends                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | AA_domain_ranges      | DBD:220-251,260-271       | semi-colon delimited list of domains                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | hugo_names            | KRAS                      | hugo gene name                                            |
    +-----------------------+---------------------------+-----------------------------------------------------------+

    Args:
        filepath (str): path to the input tab-delimited file

    Returns:
        :class:`dict` of :class:`str` and :class:`list` of :any:`Gene`: a dictionary keyed by chromosome name with
            values of list of genes on the chromosome

    Example:
        >>> ref = load_reference_genes('filename')
        >>> ref['1']
        [Gene(), Gene(), ....]
    """
    def parse_exon_list(row):
        if not row:
            return []
        exons = []
        for temp in re.split('[; ]', row):
            try:
                s, t = temp.split('-')
                exons.append(Exon(int(s), int(t)))
            except Exception as err:
                if verbose:
                    print('exon error:', repr(temp), repr(err))
        return exons

    def parse_domain_list(row):
        if not row:
            return []
        domains = []
        for d in row.split(';'):
            try:
                name, temp = d.rsplit(':')
                temp = temp.split(',')
                temp = [x.split('-') for x in temp]
                temp = [(int(x), int(y)) for x, y in temp]
                temp = Interval.min_nonoverlapping(*temp)
                d = Domain(name, temp)
                domains.append(d)
            except Exception as err:
                if verbose:
                    print('error in domain:', d, row, repr(err))
        return domains

    def nullable_int(row):
        try:
            row = int(row)
        except ValueError:
            row = TSV.null(row)
        return row

    def parse_strand(row):
        if row in ['-1', '-']:
            return STRAND.NEG
        elif row in ['1', '+', '+1']:
            return STRAND.POS
        raise ValueError('cast to strand failed')

    header, rows = TSV.read_file(
        filepath,
        require=[
            'ensembl_gene_id',
            'chr',
            'ensembl_transcript_id'
        ],
        add={
            'cdna_coding_start': 'null',
            'cdna_coding_end': 'null',
            'AA_domain_ranges': '',
            'genomic_exon_ranges': '',
            'hugo_names': '',
            'transcript_genomic_start': 'null',
            'transcript_genomic_end': 'null'
        },
        cast={
            'genomic_exon_ranges': parse_exon_list,
            'AA_domain_ranges': parse_domain_list,
            'cdna_coding_end': nullable_int,
            'cdna_coding_start': nullable_int,
            'transcript_genomic_end': nullable_int,
            'transcript_genomic_start': nullable_int,
            'strand': parse_strand,
            'gene_start': int,
            'gene_end': int
        }
    )
    genes = {}
    for row in rows:
        g = Gene(
            row['chr'],
            row['gene_start'],
            row['gene_end'],
            name=row['ensembl_gene_id'],
            strand=row['strand'],
            aliases=row['hugo_names'].split(';')
        )
        if g.name in genes:
            g = genes[g.name]
        else:
            genes[g.name] = g
        try:
            exons = row['genomic_exon_ranges']
            if len(exons) == 0:
                exons = [(row['transcript_genomic_start'], row['transcript_genomic_end'])]
            ust = usTranscript(
                name=row['ensembl_transcript_id'],
                gene=g,
                exons=exons
            )
            spl_patts = ust.generate_splicing_patterns()
            if 1 != len(spl_patts):
                raise AssertionError('expected 1 splicing pattern for reference loaded transcripts but found', len(spl_patts))
            t = Transcript(ust, spl_patts[0])
            ust.spliced_transcripts.append(t)
            if row['cdna_coding_start'] is not None:
                tx = Translation(
                    row['cdna_coding_start'], row['cdna_coding_end'], transcript=t, domains=row['AA_domain_ranges'])
                t.translations.append(tx)
            g.transcripts.append(ust)
        except Exception as err:
            print('failed loading', row['ensembl_transcript_id'], row['hugo_names'].split(';'), repr(err))

    ref = {}
    for g in genes.values():
        if g.chr not in ref:
            ref[g.chr] = []
        ref[g.chr].append(g)
    return ref


def load_reference_genome(filename):
    """
    Args:
        filename (str): the path to the file containing the input fasta genome

    Returns:
        :class:`dict` of :class:`str` and :class:`Bio.SeqRecord`: a dictionary representing the sequences in the
            fasta file
    """
    HUMAN_REFERENCE_GENOME = None
    with open(filename, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    return HUMAN_REFERENCE_GENOME
