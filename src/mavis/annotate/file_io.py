"""
module which holds all functions relating to loading reference files
"""
import json
import os
import re
from typing import Callable, Dict, List, Optional

import pandas as pd
from Bio import SeqIO
from snakemake.utils import validate as snakemake_validate

from ..constants import CODON_SIZE, GIEMSA_STAIN, START_AA, STOP_AA, STRAND, translate
from ..interval import Interval
from ..types import ReferenceAnnotations, ReferenceGenome
from ..util import logger
from .base import BioInterval, ReferenceName
from .genomic import Exon, Gene, PreTranscript, Template, Transcript
from .protein import Domain, Translation


def load_masking_regions(*filepaths: str) -> Dict[str, List[BioInterval]]:
    """
    reads a file of regions. The expect input format for the file is tab-delimited and
    the header should contain the following columns

    - chr: the chromosome
    - start: start of the region, 1-based inclusive
    - end: end of the region, 1-based inclusive
    - name: the name/label of the region

    For example:

    .. code-block:: text

        #chr    start       end         name
        chr20   25600000    27500000    centromere

    Args:
        filepath: path to the input tab-delimited file
    Returns:
        a dictionary keyed by chromosome name with values of lists of regions on the chromosome
    """
    regions: Dict[str, List[BioInterval]] = {}
    for filepath in filepaths:
        df = pd.read_csv(
            filepath, sep='\t', dtype={'chr': str, 'start': int, 'end': int, 'name': str}
        )
        for col in ['chr', 'start', 'end', 'name']:
            if col not in df:
                raise KeyError(f'missing required column ({col})')
        df['chr'] = df['chr'].apply(lambda c: ReferenceName(c))
        for row in df.to_dict('records'):
            mask_region = BioInterval(
                reference_object=row['chr'], start=row['start'], end=row['end'], name=row['name']
            )
            regions.setdefault(mask_region.reference_object, []).append(mask_region)
    return regions


def load_annotations(
    *filepaths: str,
    reference_genome: Optional[ReferenceGenome] = None,
    best_transcripts_only: bool = False,
) -> Dict[str, List[Gene]]:
    """
    loads gene models from an input file. Expects a tabbed or json file.

    Args:
        filepath: path to the input file
        reference_genome: dict of reference sequence by template/chr name

    Returns:
        lists of genes keyed by chromosome name
    """
    total_annotations: Dict[str, List[Gene]] = {}

    for filename in filepaths:
        data = None

        with open(filename) as fh:
            data = json.load(fh)

        current_annotations = parse_annotations_json(
            data,
            reference_genome=reference_genome,
            best_transcripts_only=best_transcripts_only,
        )

        for chrom in current_annotations:
            for gene in current_annotations[chrom]:
                total_annotations.setdefault(chrom, []).append(gene)
    return total_annotations


def parse_annotations_json(
    data,
    reference_genome: Optional[ReferenceGenome] = None,
    best_transcripts_only=False,
) -> ReferenceAnnotations:
    """
    parses a json of annotation information into annotation objects
    """
    try:
        snakemake_validate(
            data,
            os.path.join(os.path.dirname(__file__), 'annotations_schema.json'),
        )
    except Exception as err:
        short_msg = '. '.join(
            [line for line in str(err).split('\n') if line.strip()][:3]
        )  # these can get super long
        raise AssertionError(short_msg)

    genes_by_chr: ReferenceAnnotations = {}
    tx_skipped = 0
    domain_errors = 0

    for gene_dict in data['genes']:

        gene = Gene(
            chr=gene_dict['chr'],
            start=gene_dict['start'],
            end=gene_dict['end'],
            name=gene_dict['name'],
            aliases=gene_dict['aliases'],
            strand=gene_dict['strand'],
        )

        has_best = False
        for transcript in gene_dict['transcripts']:
            exons = []
            for ex in transcript.get('exons', []):
                exons.append(
                    Exon(strand=gene.strand, start=ex['start'], end=ex['end'], name=ex.get('name'))
                )
            if not exons:
                exons = [Exon(transcript['start'], transcript['end'], strand=gene.strand)]
            pre_transcript = PreTranscript(
                name=transcript['name'],
                gene=gene,
                exons=exons,
                is_best_transcript=transcript['is_best_transcript'],
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

                for translation in transcript.get('translations', []):
                    try:
                        if (
                            'cdna_coding_end' not in translation
                            or 'cdna_coding_start' not in translation
                        ):
                            if 'cdna_coding_end' not in translation:
                                translation['cdna_coding_end'] = spl_tx.convert_genomic_to_cdna(
                                    translation['end']
                                )
                            if 'cdna_coding_start' not in translation:
                                translation['cdna_coding_start'] = spl_tx.convert_genomic_to_cdna(
                                    translation['start']
                                )

                            if gene.strand == STRAND.NEG:
                                translation['cdna_coding_start'], translation['cdna_coding_end'] = (
                                    translation['cdna_coding_end'],
                                    translation['cdna_coding_start'],
                                )

                    except IndexError as err:
                        raise IndexError(
                            f'Invalid specification of CDS ({translation["name"]}: {translation["start"]}-{translation["end"]}) '
                            f'region on transcript ({transcript["name"]}: {transcript["start"]}-{transcript["end"]}): {err}'
                        )

                    tx_length = (
                        translation['cdna_coding_end'] - translation['cdna_coding_start'] + 1
                    )
                    # check that the translation makes sense before including it
                    if tx_length % CODON_SIZE != 0:
                        tx_skipped += 1
                        logger.debug(
                            f'Ignoring translation ({translation.get("name")}). The translated region is not a multiple of three (length={tx_length})'
                        )
                        continue
                    tx_length = tx_length // CODON_SIZE
                    domains = []
                    for dom in translation.get('domains', []):
                        try:
                            regions = [Interval(r['start'], r['end']) for r in dom['regions']]
                            regions = Interval.min_nonoverlapping(*regions)
                            for region in regions:
                                if region.start < 1 or region.end > tx_length:
                                    raise AssertionError(
                                        f'region ({dom["name"]}:{region.start}-{region.end}) cannot be outside the translated length ({tx_length})'
                                    )
                            domains.append(
                                Domain(
                                    name=dom['name'],
                                    data={'desc': dom.get('desc', None)},
                                    regions=regions,
                                )
                            )
                        except AssertionError as err:
                            domain_errors += 1
                            logger.debug(repr(err))
                    translation = Translation(
                        start=translation['cdna_coding_start'],
                        end=translation['cdna_coding_end'],
                        transcript=spl_tx,
                        domains=domains,
                        name=translation.get('name'),
                    )
                    if reference_genome and gene.chr in reference_genome:
                        # get the sequence near here to see why these are wrong?
                        seq = pre_transcript.get_cdna_seq(spl_tx.splicing_pattern, reference_genome)
                        met = seq[translation.start - 1 : translation.start + 2]
                        stop = seq[translation.end - CODON_SIZE : translation.end]
                        if translate(met) != START_AA or translate(stop) != STOP_AA:
                            logger.warning(
                                'Sequence error. The sequence computed from the reference does look like a valid translation'
                            )
                            continue
                    spl_tx.translations.append(translation)
        if not best_transcripts_only or has_best:
            genes_by_chr.setdefault(gene.chr, []).append(gene)
    if tx_skipped:
        logger.warning(
            f'Skipped {tx_skipped} translations where the CDS length was not a multiple of 3'
        )
    if domain_errors:
        logger.warning(
            f'Skipped {domain_errors} domains due to errors (coordinates defined outside the translated region)'
        )
    return genes_by_chr


def load_reference_genome(*filepaths: str) -> ReferenceGenome:
    """
    Args:
        filepaths: the paths to the files containing the input fasta genomes

    Returns:
        a dictionary representing the sequences in the fasta file
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
                    'loaded'.format(template_name, truncated)
                )
            reference_genome.setdefault(truncated, reference_genome[template_name].upper())
        else:
            prefixed = 'chr' + template_name
            if prefixed in reference_genome:
                raise KeyError(
                    'template names {} and {} are considered equal but both have been defined in the reference'
                    'loaded'.format(template_name, prefixed)
                )
            reference_genome.setdefault(prefixed, reference_genome[template_name].upper())
        reference_genome[template_name] = reference_genome[template_name].upper()

    return reference_genome


def load_templates(*filepaths: str) -> Dict[str, Template]:
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

        chr1    0   2300000 p36.33  gneg
        chr1    2300000 5400000 p36.32  gpos25

    Returns:
        templates loaded
    """
    header = ['name', 'start', 'end', 'band_name', 'giemsa_stain']
    templates: Dict[str, Template] = {}

    for filename in filepaths:
        df = pd.read_csv(
            filename,
            sep='\t',
            dtype={
                'start': int,
                'end': int,
                'name': str,
                'band_name': str,
                'giemsa_stain': str,
            },
            names=header,
            comment='#',
        )
        df['giemsa_stain'].apply(lambda v: GIEMSA_STAIN.enforce(v))

        bands_by_template: Dict[str, List[BioInterval]] = {}
        for row in df.to_dict('records'):
            band = BioInterval(None, row['start'] + 1, row['end'], name=row['band_name'], data=row)
            bands_by_template.setdefault(row['name'], []).append(band)

        for tname, bands in bands_by_template.items():
            start = min([b.start for b in bands])
            end = max([b.end for b in bands])
            end = Template(tname, start, end, bands=bands)
            templates[end.name] = end
    return templates


class ReferenceFile:
    # store loaded file to avoid re-loading
    CACHE = {}  # type: ignore

    LOAD_FUNCTIONS: Dict[str, Optional[Callable]] = {
        'annotations': load_annotations,
        'reference_genome': load_reference_genome,
        'masking': load_masking_regions,
        'template_metadata': load_templates,
        'dgv_annotation': load_masking_regions,
        'aligner_reference': None,
    }
    """dict: Mapping of file types (based on ENV name) to load functions"""

    def __init__(
        self,
        file_type: str,
        *filepaths: str,
        eager_load: bool = False,
        assert_exists: bool = False,
        **opt,
    ):
        """
        Args:
            *filepaths: list of paths to load
            file_type: Type of file to load
            eager_load: load the files immediately
            assert_exists: check that all files exist
            **opt: key word arguments to be passed to the load function and used as part of the file cache key

        Raises
            FileNotFoundError: when assert_exists and an input does not exist
        """
        self.name = sorted(filepaths)
        self.file_type = file_type
        self.key = tuple(
            self.name + sorted(list(opt.items()))
        )  # freeze the input state so we know when to reload
        self.content = None
        self.opt = opt
        self.loader = self.LOAD_FUNCTIONS[self.file_type]
        if assert_exists:
            self.files_exist()
        if eager_load:
            self.load()

    def __repr__(self):
        cls = self.__class__.__name__
        return '{}(file_type={}, files={}, loaded={}, content={})'.format(
            cls, self.file_type, self.name, self.content is not None, object.__repr__(self.content)
        )

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
                logger.info(f'cached content: {self.name}')
            self.content = ReferenceFile.CACHE[self.key].content
            return self
        self.files_exist()
        try:
            logger.info(f'loading: {self.name}')
            self.content = self.loader(*self.name, **self.opt)
            ReferenceFile.CACHE[self.key] = self
        except Exception as err:
            message = 'Error in loading files: {}. {}'.format(', '.join(self.name), err)
            raise err.__class__(message)
        return self

    @classmethod
    def load_from_config(cls, config, file_type: str, **kwargs):
        return ReferenceFile(file_type, *config.get(f'reference.{file_type}', []), **kwargs)
