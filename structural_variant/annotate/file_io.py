import TSV
import re
from .genomic import Gene, Transcript, usTranscript, Exon, Template
from .base import BioInterval
from .protein import Domain, Translation
from ..interval import Interval
from ..constants import STRAND, GIESMA_STAIN, CODON_SIZE, translate, START_AA, STOP_AA
from Bio import SeqIO
import json


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


def load_reference_genes(filepath, verbose=True, REFERENCE_GENOME=None, filetype=None):
    if filetype is None:
        m = re.match('.*\.(?P<ext>tsv|tab|json)$', filepath)
        if m:
            filetype = m.group('ext')

    if filetype == 'json':
        return _load_reference_genes_json(filepath, verbose, REFERENCE_GENOME)
    elif filetype == 'tab' or filetype == 'tsv':
        return _load_reference_genes_tabbed(filepath, verbose, REFERENCE_GENOME)
    else:
        raise NotImplementedError('unsupported filetype:', filetype, filepath)


def _load_reference_genes_json(filepath, verbose=True, REFERENCE_GENOME=None):
    genes_by_chr = dict()
    with open(filepath) as fh:
        data = json.load(fh)
        for gene in data['genes']:
            if gene['strand'] == '1':
                gene['strand'] = STRAND.POS
            elif gene['strand'] == '-1':
                gene['strand'] = STRAND.NEG
            else:
                raise AssertionError('input has unexpected form. strand must be 1 or -1 but found', gene['strand'])

            g = Gene(
                chr=gene['chr'],
                start=gene['start'],
                end=gene['end'],
                name=gene['name'],
                aliases=gene['aliases'],
                strand=gene['strand']
            )
            genes_by_chr.setdefault(g.chr, []).append(g)
            print(g, gene['aliases'])

            for transcript in gene['transcripts']:
                if transcript['is_best_transcript'] == 'true':
                    transcript['is_best_transcript'] = True
                elif transcript['is_best_transcript'] == 'false':
                    transcript['is_best_transcript'] = False

                exons = [Exon(**ex) for ex in transcript['exons']]
                if len(exons) == 0:
                    exons[(transcript['start'], transcript['end'])]
                ust = usTranscript(
                    name=transcript['name'],
                    gene=g,
                    exons=exons,
                    is_best_transcript=transcript['is_best_transcript']
                )
                g.transcripts.append(ust)

                if transcript['cdna_coding_end'] is None or transcript['cdna_coding_start'] is None:
                    continue

                for spl_patt in ust.generate_splicing_patterns():
                    # make splice transcripts and translations
                    t = Transcript(ust, spl_patt)
                    ust.spliced_transcripts.append(t)

                    domains = []
                    for dom in transcript['domains']:
                        regions = [Interval(r['start'], r['end']) for r in dom['regions']]
                        regions = Interval.min_nonoverlapping(*regions)
                        domains.append(
                            Domain(name=dom['name'], data={'desc': dom['desc']}, regions=regions)
                        )
                    tx = Translation(
                        transcript['cdna_coding_start'], transcript['cdna_coding_end'], transcript=t, domains=domains
                    )
                    # check that the translation makes sense before including it
                    if len(tx) % CODON_SIZE != 0:
                        continue
                    if REFERENCE_GENOME and g['chr'] in REFERENCE_GENOME:
                        # get the sequence near here to see why these are wrong?
                        s = ust.get_cdna_sequence(t.splicing_pattern, REFERENCE_GENOME)
                        m = s[tx.start - 1:tx.start + 2]
                        stop = s[tx.end - CODON_SIZE: tx.end]
                        if translate(m) != START_AA or translate(stop) != STOP_AA:
                            continue
                    t.translations.append(tx)
    return genes_by_chr


def _load_reference_genes_tabbed(filepath, verbose=True, REFERENCE_GENOME=None):
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

    Warning:
        does not load translations unless then start with 'M', end with '*' and have a length of multiple 3
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
            'transcript_genomic_end': 'null',
            'best_ensembl_transcript_id': 'null'
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
    failures = 0
    genes = {}
    for row in rows:
        g = Gene(
            row['chr'],
            row['gene_start'],
            row['gene_end'],
            name=row['ensembl_gene_id'],
            strand=row['strand'],
            aliases=row['hugo_names'].split(';') if row['hugo_names'] else []
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
            if row['best_ensembl_transcript_id'] == row['ensembl_transcript_id']:
                ust.is_best_transcript = True
            assert(ust.start >= g.start)
            assert(ust.end <= g.end)
            spl_patts = ust.generate_splicing_patterns()
            if 1 != len(spl_patts):
                raise AssertionError('expected 1 splicing pattern for reference loaded transcripts but found', len(spl_patts))
            t = Transcript(ust, spl_patts[0])
            ust.spliced_transcripts.append(t)
            try:
                if row['cdna_coding_start'] is not None:
                    tx = Translation(
                        row['cdna_coding_start'], row['cdna_coding_end'], transcript=t, domains=row['AA_domain_ranges'])

                    if len(tx) % CODON_SIZE != 0:
                        raise AssertionError('the length of the coding sequence is not a multiple of {}. This violates the input assumption.'.format(CODON_SIZE), len(tx) % CODON_SIZE)

                    if REFERENCE_GENOME and row['chr'] in REFERENCE_GENOME:
                        # get the sequence near here to see why these are wrong?
                        s = ust.get_cdna_sequence(t.splicing_pattern, REFERENCE_GENOME)
                        m = s[tx.start - 1:tx.start + 2]
                        stop = s[tx.end - CODON_SIZE: tx.end]
                        if translate(m) != START_AA or translate(stop) != STOP_AA:
                            raise AssertionError('translation sequence check failure: M={}, *={}'.format(translate(m), translate(stop)))

                    t.translations.append(tx)
            except AssertionError:
                pass
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


def load_templates(filename):
    """
    primarily useful if template drawings are required and is not necessary otherwise
    assumes the input file is 0-indexed with [start,end) style

    Args:
        filename (str): the path to the file with the cytoband template information

    Returns:
        :class:`list` of :class:`Template`: list of the templates loaded
    """
    header = ['name', 'start', 'end', 'band_name', 'giesma_stain']
    header, rows = TSV.read_file(
        filename,
        header=header,
        cast={'start': int, 'end': int, 'name': lambda x: re.sub('^chr', '', x)},
        _in={'giesma_stain': GIESMA_STAIN}
    )

    bands_by_template = {}
    for row in rows:
        b = BioInterval(None, row['start'] + 1, row['end'], name=row['band_name'], data=row)
        bands_by_template.setdefault(row['name'], []).append(b)

    templates = dict()
    for tname, bands in bands_by_template.items():
        s = min([b.start for b in bands])
        t = max([b.end for b in bands])
        t = Template(tname, s, t, bands=bands)
        templates[t.name] = t
    return templates
