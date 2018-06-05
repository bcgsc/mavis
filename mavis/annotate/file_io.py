"""
module which holds all functions relating to loading reference files
"""
import json
import re
import warnings
import os

from Bio import SeqIO
import tab

from .base import BioInterval, ReferenceName
from .genomic import Exon, Gene, Template, Transcript, PreTranscript
from .protein import Domain, Translation
from ..constants import CODON_SIZE, GIEMSA_STAIN, START_AA, STOP_AA, STRAND, translate
from ..interval import Interval
from ..util import DEVNULL, LOG, filepath, WeakMavisNamespace


REFERENCE_DEFAULTS = WeakMavisNamespace()
REFERENCE_DEFAULTS.add(
    'template_metadata', [], cast_type=filepath, listable=True,
    defn='file containing the cytoband template information. Used for illustrations only')
REFERENCE_DEFAULTS.add(
    'masking', [], cast_type=filepath, listable=True,
    defn='file containing regions for which input events overlapping them are dropped prior to validation')
REFERENCE_DEFAULTS.add(
    'annotations', [], cast_type=filepath, listable=True,
    defn='path to the reference annotations of genes, transcript, exons, domains, etc')
REFERENCE_DEFAULTS.add(
    'aligner_reference', None, cast_type=filepath, nullable=True,
    defn='path to the aligner reference file used for aligning the contig sequences')
REFERENCE_DEFAULTS.add(
    'dgv_annotation', [], cast_type=filepath, listable=True,
    defn='Path to the dgv reference processed to look like the cytoband file.')
REFERENCE_DEFAULTS.add(
    'reference_genome', [], cast_type=filepath, listable=True,
    defn='Path to the human reference genome fasta file')


def load_masking_regions(*filepaths):
    """
    reads a file of regions. The expect input format for the file is tab-delimited and
    the header should contain the following columns

    - chr: the chromosome
    - start: start of the region, 1-based inclusive
    - end: end of the region, 1-based inclusive
    - name: the name/label of the region

    For example:

    .. code-block:: text

        #chr    start   end     name
        chr20	25600000	27500000	centromere

    Args:
        filepath (str): path to the input tab-delimited file
    Returns:
        :class:`dict` of :class:`list` of :class:`BioInterval` by :class:`str`: a dictionary keyed by chromosome name with values of lists of regions on the chromosome

    Example:
        >>> m = load_masking_regions('filename')
        >>> m['1']
        [BioInterval(), BioInterval(), ...]
    """
    regions = {}
    for filepath in filepaths:
        _, rows = tab.read_file(
            filepath,
            require=['chr', 'start', 'end', 'name'],
            cast={'start': int, 'end': int, 'chr': ReferenceName}
        )
        for row in rows:
            mask_region = BioInterval(reference_object=row['chr'], start=row['start'], end=row['end'], name=row['name'])
            regions.setdefault(mask_region.reference_object, []).append(mask_region)
    return regions


def load_reference_genes(*pos, **kwargs):
    """
    *Deprecated* Use :func:`load_annotations` instead
    """
    warnings.warn('this function has been replaced by load_annotations', DeprecationWarning)
    return load_annotations(*pos, **kwargs)


def load_annotations(*filepaths, warn=DEVNULL, reference_genome=None, best_transcripts_only=False):
    """
    loads gene models from an input file. Expects a tabbed or json file.

    Args:
        filepath (str): path to the input file
        verbose (bool): output extra information to stdout
        reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by
            template/chr name
        filetype (str): json or tab/tsv. only required if the file type can't be interpolated from the path extension

    Returns:
        :class:`dict` of :class:`list` of :class:`~mavis.annotate.genomic.Gene` by :class:`str`: lists of genes keyed by chromosome name
    """
    total_annotations = {}

    for filename in filepaths:
        data = None

        if filename.endswith('.tab') or filename.endswith('.tsv'):
            data = convert_tab_to_json(filename, warn)
        else:
            with open(filename) as fh:
                data = json.load(fh)

        current_annotations = parse_annotations_json(
            data,
            reference_genome=reference_genome,
            best_transcripts_only=best_transcripts_only,
            warn=warn)

        for chrom in current_annotations:
            for gene in current_annotations[chrom]:
                total_annotations.setdefault(chrom, []).append(gene)
    return total_annotations


def parse_annotations_json(data, reference_genome=None, best_transcripts_only=False, warn=DEVNULL):
    """
    parses a json of annotation information into annotation objects
    """
    genes_by_chr = {}
    for gene_dict in data['genes']:
        if gene_dict['strand'] in ['1', '+']:
            gene_dict['strand'] = STRAND.POS
        elif gene_dict['strand'] in ['-1', '-']:
            gene_dict['strand'] = STRAND.NEG
        else:
            raise AssertionError('input has unexpected form. strand must be 1 or -1 but found', gene_dict['strand'])

        gene = Gene(
            chr=gene_dict['chr'],
            start=gene_dict['start'],
            end=gene_dict['end'],
            name=gene_dict['name'],
            aliases=gene_dict['aliases'],
            strand=gene_dict['strand']
        )

        has_best = False
        for transcript in gene_dict['transcripts']:
            transcript['is_best_transcript'] = tab.cast_boolean(transcript['is_best_transcript'])
            transcript.setdefault('exons', [])
            exons = [Exon(strand=gene.strand, **ex) for ex in transcript['exons']]
            if not exons:
                exons = [(transcript['start'], transcript['end'])]
            pre_transcript = PreTranscript(
                name=transcript['name'],
                gene=gene,
                exons=exons,
                is_best_transcript=transcript['is_best_transcript']
            )
            if pre_transcript.is_best_transcript:
                has_best = True
            if best_transcripts_only and not pre_transcript.is_best_transcript:
                continue
            gene.transcripts.append(pre_transcript)

            for spl_patt in pre_transcript.generate_splicing_patterns():
                # make splice transcripts and translations
                spl_tx = Transcript(pre_transcript, spl_patt)
                pre_transcript.spliced_transcripts.append(spl_tx)

                if transcript.get('cdna_coding_end', None) is None or transcript.get('cdna_coding_start', None) is None:
                    continue
                tx_length = transcript['cdna_coding_end'] - transcript['cdna_coding_start'] + 1
                # check that the translation makes sense before including it
                if tx_length % CODON_SIZE != 0:
                    warn('Ignoring translation. The translated region is not a multiple of three')
                    continue
                tx_length = tx_length // CODON_SIZE
                domains = []
                for dom in transcript['domains']:
                    try:
                        regions = [Interval(r['start'], r['end']) for r in dom['regions']]
                        regions = Interval.min_nonoverlapping(*regions)
                        for region in regions:
                            if region.start < 1 or region.end > tx_length:
                                raise AssertionError('region cannot be outside the translated length')
                        domains.append(
                            Domain(name=dom['name'], data={'desc': dom.get('desc', None)}, regions=regions)
                        )
                    except AssertionError as err:
                        warn(repr(err))
                translation = Translation(
                    transcript['cdna_coding_start'], transcript['cdna_coding_end'], transcript=spl_tx, domains=domains
                )
                if reference_genome and gene.chr in reference_genome:
                    # get the sequence near here to see why these are wrong?
                    seq = pre_transcript.get_cdna_seq(spl_tx.splicing_pattern, reference_genome)
                    met = seq[translation.start - 1:translation.start + 2]
                    stop = seq[translation.end - CODON_SIZE: translation.end]
                    if translate(met) != START_AA or translate(stop) != STOP_AA:
                        warn(
                            'Sequence error. The sequence computed from the reference does look like '
                            'a valid translation'
                        )
                        continue
                spl_tx.translations.append(translation)
        if not best_transcripts_only or has_best:
            genes_by_chr.setdefault(gene.chr, []).append(gene)
    return genes_by_chr


def convert_tab_to_json(filepath, warn=DEVNULL):
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
        :class:`dict` of :class:`list` of :any:`Gene` by :class:`str`: a dictionary keyed by chromosome name with values of list of genes on the chromosome

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
                start, end = temp.split('-')
                exons.append({'start': int(start), 'end': int(end)})
            except Exception as err:
                warn('exon error:', repr(temp), repr(err))
        return exons

    def parse_domain_list(row):
        if not row:
            return []
        domains = []
        for domain in row.split(';'):
            try:
                name, temp = domain.rsplit(':')
                temp = temp.split(',')
                temp = [x.split('-') for x in temp]
                regions = [{'start': int(x), 'end': int(y)} for x, y in temp]
                domains.append({'name': name, 'regions': regions})
            except Exception as err:
                warn('error in domain:', domain, row, repr(err))
        return domains

    def nullable_int(row):
        try:
            row = int(row)
        except ValueError:
            row = tab.cast_null(row)
        return row

    _, rows = tab.read_file(
        filepath,
        require=[
            'ensembl_gene_id',
            'chr',
            'ensembl_transcript_id'
        ],
        add_default={
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
            'gene_start': int,
            'gene_end': int
        }
    )
    genes = {}
    for row in rows:
        gene = {
            'chr': row['chr'],
            'start': row['gene_start'],
            'end': row['gene_end'],
            'name': row['ensembl_gene_id'],
            'strand': row['strand'],
            'aliases': row['hugo_names'].split(';') if row['hugo_names'] else [],
            'transcripts': []
        }
        if gene['name'] not in genes:
            genes[gene['name']] = gene
        else:
            gene = genes[gene['name']]

        transcript = {
            'is_best_transcript': row['best_ensembl_transcript_id'] == row['ensembl_transcript_id'],
            'name': row['ensembl_transcript_id'],
            'exons': row['genomic_exon_ranges'],
            'domains': row['AA_domain_ranges'],
            'start': row['transcript_genomic_start'],
            'end': row['transcript_genomic_end'],
            'cdna_coding_start': row['cdna_coding_start'],
            'cdna_coding_end': row['cdna_coding_end'],
            'aliases': []
        }
        gene['transcripts'].append(transcript)

    return {'genes': genes.values()}


def load_reference_genome(*filepaths):
    """
    Args:
        filepaths (list of str): the paths to the files containing the input fasta genomes

    Returns:
        :class:`dict` of :class:`Bio.SeqRecord` by :class:`str`: a dictionary representing the sequences in the fasta file
    """
    reference_genome = {}
    for filename in filepaths:
        with open(filename, 'rU') as fh:
            for chrom, seq in SeqIO.to_dict(SeqIO.parse(fh, 'fasta')).items():
                if chrom in reference_genome:
                    raise KeyError('Duplicate chromosome name', chrom, filename)
                reference_genome[chrom] = seq

    names = list(reference_genome.keys())

    # to fix hg38 issues
    for template_name in names:
        if template_name.startswith('chr'):
            truncated = re.sub('^chr', '', template_name)
            if truncated in reference_genome:
                raise KeyError(
                    'template names {} and {} are considered equal but both have been defined in the reference'
                    'loaded'.format(template_name, truncated))
            reference_genome.setdefault(truncated, reference_genome[template_name].upper())
        else:
            prefixed = 'chr' + template_name
            if prefixed in reference_genome:
                raise KeyError(
                    'template names {} and {} are considered equal but both have been defined in the reference'
                    'loaded'.format(template_name, prefixed))
            reference_genome.setdefault(prefixed, reference_genome[template_name].upper())
        reference_genome[template_name] = reference_genome[template_name].upper()

    return reference_genome


def load_templates(*filepaths):
    """
    primarily useful if template drawings are required and is not necessary otherwise
    assumes the input file is 0-indexed with [start,end) style. Columns are expected in
    the following order, tab-delimited. A header should not be given

    1. name
    2. start
    3. end
    4. band_name
    5. giemsa_stain

    for example

    .. code-block:: text

        chr1	0	2300000	p36.33	gneg
        chr1	2300000	5400000	p36.32	gpos25

    Args:
        filename (str): the path to the file with the cytoband template information

    Returns:
        :class:`list` of :class:`Template`: list of the templates loaded

    """
    header = ['name', 'start', 'end', 'band_name', 'giemsa_stain']
    templates = {}

    for filename in filepaths:
        header, rows = tab.read_file(
            filename,
            header=header,
            cast={'start': int, 'end': int},
            in_={'giemsa_stain': GIEMSA_STAIN.values()}
        )

        bands_by_template = {}
        for row in rows:
            band = BioInterval(None, row['start'] + 1, row['end'], name=row['band_name'], data=row)
            bands_by_template.setdefault(row['name'], []).append(band)

        for tname, bands in bands_by_template.items():
            start = min([b.start for b in bands])
            end = max([b.end for b in bands])
            end = Template(tname, start, end, bands=bands)
            templates[end.name] = end
    return templates


class ReferenceFile:

    CACHE = {}  # store loaded file to avoid re-loading

    LOAD_FUNCTIONS = {
        'annotations': load_annotations,
        'reference_genome': load_reference_genome,
        'masking': load_masking_regions,
        'template_metadata': load_templates,
        'dgv_annotation': load_masking_regions,
        'aligner_reference': None
    }
    """:class:`dict`: Mapping of file types (based on ENV name) to load functions"""

    def __init__(self, file_type, *filepaths, eager_load=False, assert_exists=False, **opt):
        """
        Args:
            *filepaths (str): list of paths to load
            file_type (str): Type of file to load
            eager_load (bool=False): load the files immeadiately
            assert_exists (bool=False): check that all files exist
            **opt: key word arguments to be passed to the load function and used as part of the file cache key

        Raises
            FileNotFoundError: when assert_exists and an input does not exist
        """
        self.name = sorted(filepaths)
        self.file_type = file_type
        self.key = tuple(self.name + sorted(list(opt.items())))  # freeze the input state so we know when to reload
        self.content = None
        self.opt = opt
        self.loader = self.LOAD_FUNCTIONS[self.file_type]
        if assert_exists:
            self.files_exist()
        if eager_load:
            self.load()

    def __repr__(self):
        cls = self.__class__.__name__
        return '{}(file_type={}, files={}, loaded={}, content={})'.format(cls, self.file_type, self.name, self.content is not None, object.__repr__(self.content))

    def files_exist(self, not_empty=False):
        if not_empty and not self.name:
            raise FileNotFoundError('expected files but given an empty list', self)
        for filename in self.name:
            if not os.path.exists(filename):
                raise FileNotFoundError('Missing file', filename, self)

    def __iter__(self):
        return iter(self.name)

    def is_empty(self):
        return not self.name

    def is_loaded(self):
        return False if self.content is None else True

    def load(self, ignore_cache=False, verbose=True):
        """
        load (or return) the contents of a reference file and add it to the cache if enabled
        """
        if self.content is not None:
            return self
        if self.key in ReferenceFile.CACHE and not ignore_cache:
            if verbose:
                LOG('cached content:', self.name)
            self.content = ReferenceFile.CACHE[self.key].content
            return self
        self.files_exist()
        try:
            LOG('loading:', self.name, time_stamp=True)
            self.content = self.loader(*self.name, **self.opt)
            ReferenceFile.CACHE[self.key] = self
        except Exception as err:
            message = 'Error in loading files: {}. {}'.format(', '.join(self.name), err)
            raise err.__class__(message)
        return self
