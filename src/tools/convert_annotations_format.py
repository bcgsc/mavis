import argparse
import json
import logging
import re
from typing import Dict, Tuple

import pandas as pd
from mavis.annotate.file_io import parse_annotations_json

PANDAS_DEFAULT_NA_VALUES = [
    '-1.#IND',
    '1.#QNAN',
    '1.#IND',
    '-1.#QNAN',
    '#N/A',
    'N/A',
    'NA',
    '#NA',
    'NULL',
    'NaN',
    '-NaN',
    'nan',
    '-nan',
]


GFF_GENELIKE_FEATURES = {
    'gene',
    'ncRNA_gene',
    'biological_region',
    'pseudogene',
    'enhancer',
    'promoter',
    'region',
    'protein_binding_site',
}
GFF_RNALIKE_FEATURES = {
    'rna',
    'mRNA',
    'lncRNA',
    'transcript',
    'lnc_RNA',
    'pseudogenic_transcript',
    'snRNA',
    'miRNA',
    'unconfirmed_transcript',
    'ncRNA',
    'snoRNA',
    'scRNA',
}
GFF_ALL_FEATURES = GFF_GENELIKE_FEATURES | GFF_RNALIKE_FEATURES | {'CDS', 'exon'}
GFF_ID_DELIMITER = '_'
GFF_ATTRS = [
    'Alias',
    'bound_moiety',
    'DBxref',
    'Derives_from',
    'exon_id',
    'exon_number',
    'exon_version',
    'function',
    'Gap',
    'gene_id',
    'gene_name',
    'gene_version',
    'ID',
    'Name',
    'Note',
    'old-name',
    'Ontology_term',
    'Parent',
    'product',
    'protein_id',
    'protein_version',
    'rank',
    'standard_name',
    'Target',
    'transcript_id',
    'transcript_name',
    'transcript_version',
    'version',
]
GFF_KEY_COLS = ['feature_id', 'type', 'seqid', 'strand']


def agg_strings_unique(series):
    series = series.fillna('')
    return ';'.join([s for s in series.astype(str).unique()])


def strip_empty_fields(input_obj):
    """Remove all empty string or null fields from some dictionary object to reduce the size"""

    if isinstance(input_obj, dict):
        result = {}
        for k, v in input_obj.items():
            if v == '' or v is None or (isinstance(v, list) and not len(v)):
                continue
            result[k] = strip_empty_fields(v)
        return result
    elif isinstance(input_obj, list):
        return [strip_empty_fields(v) for v in input_obj]
    return input_obj


def coerce_number_types(input_obj, fields=['start', 'end', 'coding_cdna_start', 'coding_cdna_end']):
    if isinstance(input_obj, dict):
        result = {}
        for k, v in input_obj.items():
            if k in fields and isinstance(v, str):
                if v.lower() in {'', 'null', 'none'}:
                    continue
                result[k] = int(v)
            else:
                result[k] = coerce_number_types(v)
        return result
    elif isinstance(input_obj, list):
        return [coerce_number_types(v) for v in input_obj]
    return input_obj


def convert_tab_to_json(filepath: str) -> Dict:
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
    | genomic_exon_ranges   | 100-201;334-412;779-830   | semi-colon delimited exon start/ends                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | AA_domain_ranges      | DBD:220-251,260-271       | semi-colon delimited list of domains                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | hugo_names            | KRAS                      | hugo gene name                                            |
    +-----------------------+---------------------------+-----------------------------------------------------------+

    Args:
        filepath (str): path to the input tab-delimited file

    Returns:
        Dict[str,List[Gene]]: a dictionary keyed by chromosome name with values of list of genes on the chromosome

    Warning:
        does not load translations unless then start with 'M', end with '*' and have a length of multiple 3
    """

    def parse_exon_list(row):
        if pd.isnull(row):
            return []
        exons = []
        for temp in re.split('[; ]', row):
            try:
                start, end = temp.split('-')
                exons.append({'start': int(start), 'end': int(end)})
            except Exception as err:
                logging.warning(f'exon error: {repr(temp)}, {repr(err)}')
        return exons

    def parse_domain_list(row):
        if pd.isnull(row):
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
                logging.warning(f'error in domain: {domain}, {row}, {repr(err)}')
        return domains

    skip_lines = 0
    with open(filepath, 'r') as fh:
        lines = fh.readlines()
        skip_lines = len([line for line in lines if line.startswith('##')])

    df = pd.read_csv(
        filepath,
        skiprows=skip_lines,
        dtype={
            '#ensembl_gene_id': str,
            'ensembl_transcript_id': str,
            'chr': str,
            'cdna_coding_start': pd.Int64Dtype(),
            'cdna_coding_end': pd.Int64Dtype(),
            'AA_domain_ranges': str,
            'genomic_exon_ranges': str,
            'hugo_names': str,
            'transcript_genomic_start': pd.Int64Dtype(),
            'transcript_genomic_end': pd.Int64Dtype(),
            'best_ensembl_transcript_id': str,
            'gene_start': int,
            'gene_end': int,
        },
        sep='\t',
    ).rename(columns={'#ensembl_gene_id': 'ensembl_gene_id'})

    for col in ['ensembl_gene_id', 'chr', 'ensembl_transcript_id', 'gene_start', 'gene_end']:
        if col not in df:
            raise KeyError(f'missing required column: {col}')

    for col, parser in [
        ('genomic_exon_ranges', parse_exon_list),
        ('AA_domain_ranges', parse_domain_list),
    ]:
        if col in df:
            df[col] = df[col].apply(parser)

    genes = {}
    rows = df.where(df.notnull(), None).to_dict('records')

    for row in rows:
        gene = {
            'chr': row['chr'],
            'start': int(row['gene_start']),
            'end': int(row['gene_end']),
            'name': row['ensembl_gene_id'],
            'strand': row['strand'],
            'aliases': row['hugo_names'].split(';') if row.get('hugo_names') else [],
            'transcripts': [],
        }
        if gene['strand'] in {'true', '1', '+', '+1', 'True', 1, True}:
            gene['strand'] = '+'
        elif gene['strand'] in {'false', '-1', '-', 'False', -1, False}:
            gene['strand'] = '-'
        if gene['name'] not in genes:
            genes[gene['name']] = gene
        else:
            gene = genes[gene['name']]
        is_best_transcript = (
            row.get('best_ensembl_transcript_id', row['ensembl_transcript_id'])
            == row['ensembl_transcript_id']
        )
        transcript = {
            'is_best_transcript': is_best_transcript,
            'name': row['ensembl_transcript_id'],
            'exons': row.get('genomic_exon_ranges', []),
            'start': row.get('transcript_genomic_start'),
            'end': row.get('transcript_genomic_end'),
            'aliases': [],
            'translations': [
                {
                    'domains': row.get('AA_domain_ranges', []),
                    'cdna_coding_start': row.get('cdna_coding_start'),
                    'cdna_coding_end': row.get('cdna_coding_end'),
                }
            ],
        }
        gene['transcripts'].append(transcript)

    return coerce_number_types({'genes': list(genes.values())})


def strip_id_field(feature_id) -> Tuple[str, str]:
    """
    Remove type prefix from ID if applicable
    """
    prefix_map = {k: k for k in ['gene', 'transcript', 'cds', 'exon']}
    prefix_map.update({k: 'gene' for k in GFF_GENELIKE_FEATURES})
    prefix_map.update({k: 'transcript' for k in GFF_RNALIKE_FEATURES})
    if feature_id:
        for prefix in prefix_map:
            if feature_id.lower().startswith(prefix):
                return prefix_map.get(prefix, prefix), feature_id[len(prefix) + 1 :]
    return '', feature_id


def parse_gff_id(row):
    """
    Get the unique ID of the current row/feature
    """
    _, feature_id = strip_id_field(row.ID if 'ID' in row else '')

    if not feature_id:
        if row.type == 'exon' and 'exon_id' in row:
            return row.exon_id
        elif row.type == 'gene' and 'gene_id' in row:
            return row.gene_id
        elif row.type == 'transcript' and 'transcript_id' in row:
            return row.transcript_id
        elif row.type.lower() == 'cds' and 'protein_id' in row:
            return row.protein_id
    return feature_id


def pull_alias_terms(row):
    aliases = []
    for field in ['Name', 'standard_name', 'old-name']:
        if row[field] and not pd.isnull(row[field]):
            aliases.extend(row[field].split(';'))
    if row.Alias and not pd.isnull(row.Alias):
        aliases.extend(row.Alias.split(','))
    return [a for a in aliases if a != row.feature_id]


class NumberedFeatureGenerator:
    def __init__(self):
        self.counter = 0

    def __call__(self, features, parent_id, prefix='-T'):
        result = f'{parent_id}{prefix}{self.counter}'
        while result in features:
            self.counter += 1
            result = f'{parent_id}{prefix}{self.counter}'
        return result


def split_col_into_rows(df, col, delimiter=',', new_col=None):
    """
    Given some string column in a dataframe, split the column by the delimiter and for each resulting value duplicate the existing row
    """
    if not new_col:
        new_col = col
    new_df = df.copy().reset_index()

    s = new_df[col].str.split(delimiter).apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = new_col

    if new_col == col:
        new_df = new_df.drop(columns=[new_col])
    return new_df.merge(s, left_index=True, right_index=True)


def fix_dangling_parent_reference(nodes_df, links_df):
    """
    Insert a pseudo element for any parents referenced by an element that do not already have their own line/definition

    Returns the elements to be added to the node definitions
    """
    dangling_refs = links_df.rename(
        {
            'parent_id': 'feature_id',
            'parent_type': 'type',
            'feature_id': 'child_id',
            'type': 'child_type',
        }
    ).merge(nodes_df[GFF_KEY_COLS], how='left', indicator=True)
    dangling_refs = dangling_refs[dangling_refs._merge == 'left_only']
    # now join back to its children to create coordinates that are the interval covering all connected children
    dangling_refs = dangling_refs.merge(
        nodes_df[GFF_KEY_COLS + ['start', 'end', 'row_index']].rename(
            columns={'feature_id': 'child_id', 'type': 'child_type'}
        )
    )
    dangling_refs = (
        dangling_refs.groupby(GFF_KEY_COLS)
        .agg(
            {
                'start': 'min',
                'end': 'max',
                'row_index': agg_strings_unique,
            }
        )
        .reset_index()
    )
    if dangling_refs.shape[0]:
        logging.warning(f'Inserting {dangling_refs.shape[0]} missing parent element definitions')

    return pd.concat([nodes_df, dangling_refs]).reset_index(drop=True), links_df


def fix_orphan_elements(nodes_df, links_df):
    """
    When there are non-gene elements that do not have a parent assigned to them, connect them to a
    inserted 'mock' gene instead
    """
    links_df = links_df.copy()

    links_df['_orphan'] = False
    links_df.loc[
        (links_df.parent_id == '') & (links_df.type.isin({'CDS', 'exon'})), '_orphan'
    ] = True
    links_df.loc[links_df._orphan, 'parent_id'] = 'G' + GFF_ID_DELIMITER + links_df.feature_id
    links_df.loc[links_df._orphan, 'parent_type'] = 'gene'

    new_genes_df = (
        links_df[links_df._orphan]
        .merge(nodes_df[GFF_KEY_COLS + ['start', 'end']])
        .rename(
            columns={
                'feature_id': 'child_id',
                'type': 'child_type',
                'parent_id': 'feature_id',
                'parent_type': 'type',
            }
        )
    )
    new_genes_df = (
        new_genes_df.groupby(GFF_KEY_COLS)
        .agg({'start': 'min', 'end': 'max', 'row_index': agg_strings_unique})
        .reset_index()
    )

    links_df = links_df.drop(columns=['_orphan'])
    if new_genes_df.shape[0]:
        logging.warning(
            f'Inserting {new_genes_df.shape[0]} new genes to connect to orphan elements'
        )
    return pd.concat([nodes_df, new_genes_df]).reset_index(drop=True), links_df


def insert_missing_transcripts(nodes_df, links_df):
    """
    For any cds elements with a direct parent gene, create a transcript and link them through that instead
    """
    direct_links_df = links_df[(links_df.parent_type == 'gene') & (links_df.type != 'transcript')]
    rest_links_df = links_df[(links_df.parent_type != 'gene') | (links_df.type == 'transcript')]

    src_transcript_df = direct_links_df.copy()
    src_transcript_df['feature_id'] = src_transcript_df.parent_id + GFF_ID_DELIMITER + 'T'
    src_transcript_df['type'] = 'transcript'

    tgt_transcript_df = direct_links_df.copy()
    tgt_transcript_df['parent_id'] = tgt_transcript_df.parent_id + GFF_ID_DELIMITER + 'T'
    tgt_transcript_df['parent_type'] = 'transcript'

    links_df = pd.concat([rest_links_df, src_transcript_df, tgt_transcript_df]).reset_index(
        drop=True
    )

    if direct_links_df.shape[0]:
        logging.warning(
            f'Inserting {direct_links_df.shape[0]} transcripts between lower element to gene connections'
        )

    return fix_dangling_parent_reference(nodes_df, links_df)


def validate_gff_coordinates(nodes_df, links_df):
    """
    Check that all child elements have coordinates within the coordinates of their parent elements
    """
    df = links_df.merge(nodes_df[GFF_KEY_COLS + ['start', 'end']]).merge(
        nodes_df[GFF_KEY_COLS + ['start', 'end']].rename(
            columns={
                'feature_id': 'parent_id',
                'type': 'parent_type',
                'start': 'parent_start',
                'end': 'parent_end',
            }
        )
    )
    df['error'] = False
    df.loc[(df.parent_start > df.start) | (df.parent_end < df.end), 'error'] = True

    errors = df[df.error]
    if errors.shape[0]:
        for _, row in errors.iterrows():
            logging.debug(
                f'{row.feature_id} ({row.start}-{row.end}) is not within its parent element {row.parent_id} ({row.parent_start}-{row.parent_end})'
            )
        raise ValueError(f'{errors.shape[0]} entries with impossible coordinates')


def enforce_uniq_transcript_ids(input_df) -> pd.DataFrame:
    df = input_df.copy()
    duplicates = df[df.type == 'transcript'].drop_duplicates(['seqid', 'parent_id', 'feature_id'])

    if duplicates.shape[0] == duplicates.feature_id.nunique():
        return df

    # there are some non-unique transcript IDs, make them all pre-pend the seqid
    # do not change ensembl transcript IDs since they should already be unique
    df.loc[(df.type == 'transcript') & (~df.feature_id.str.startswith('ENST')), 'feature_id'] = (
        df.seqid + GFF_ID_DELIMITER + df.feature_id
    )
    df.loc[
        (df.parent_type == 'transcript') & (~df.parent_id.str.startswith('ENST')), 'parent_id'
    ] = (df.seqid + GFF_ID_DELIMITER + df.parent_id)
    duplicates = df[df.type == 'transcript'].drop_duplicates(['seqid', 'parent_id', 'feature_id'])

    if duplicates.shape[0] == duplicates.feature_id.nunique():
        return df.copy()

    raise ValueError(
        f'Unable to enforce unique transcript IDs: ({duplicates.shape[0]},{duplicates.feature_id.nunique()})'
    )


def convert_pandas_gff_to_mavis(df) -> Dict:
    df['error'] = ''
    df.loc[~df.type.isin(GFF_ALL_FEATURES), 'error'] = 'unrecognized type ' + df.type
    df = split_col_into_rows(df, 'Parent', ',')
    # simplify the type
    df['biotype'] = df.type.fillna('')

    def simplify_type(t):
        if t in GFF_GENELIKE_FEATURES:
            return 'gene'
        elif t in GFF_RNALIKE_FEATURES:
            return 'transcript'
        return t

    df['type'] = df.type.apply(simplify_type).fillna('')
    df['parent_type'] = (
        df.Parent.apply(lambda x: strip_id_field(x)[0]).fillna('').apply(simplify_type)
    )
    df['parent_id'] = df.Parent.apply(lambda x: strip_id_field(x)[1]).fillna('')
    df.loc[df.type == 'gene', 'parent_type'] = 'seq'
    df.loc[df.type == 'gene', 'parent_id'] = df.seqid

    if df[df.error != ''].shape[0]:
        logging.warning(
            f'dropping {df[df.error != ""].shape[0]} features that did not match an expected type: {df[df.error != ""].type.unique()}'
        )
    df = df[df.error == '']

    if df[df.feature_id == ''].shape[0]:
        logging.warning(f'dropping {df[df.feature_id == ""].shape[0]} rows for missing ID')
    df = df[df.feature_id != '']
    df['regions'] = df.start.astype(str) + '-' + df.end.astype(str)

    # use the feature key to group elements that are discontinuous
    links_df = (
        df.sort_values(['seqid', 'start'])
        .groupby(GFF_KEY_COLS + ['parent_type', 'parent_id'])
        .agg({'row_index': agg_strings_unique})
        .reset_index()
    )
    nodes_df = (
        df.sort_values(['seqid', 'start'])
        .groupby(GFF_KEY_COLS)
        .agg(
            {
                'start': 'min',
                'end': 'max',
                'regions': agg_strings_unique,
                'version': agg_strings_unique,
                'Note': agg_strings_unique,
                'Name': agg_strings_unique,
                'Alias': agg_strings_unique,
                'biotype': agg_strings_unique,
                'exon_number': agg_strings_unique,
                'row_index': agg_strings_unique,
                'source': agg_strings_unique,
                'standard_name': agg_strings_unique,
                'old-name': agg_strings_unique,
            }
        )
        .reset_index()
    )
    nodes_df, links_df = fix_dangling_parent_reference(nodes_df, links_df)
    nodes_df, links_df = fix_orphan_elements(nodes_df, links_df)
    nodes_df, links_df = insert_missing_transcripts(nodes_df, links_df)
    validate_gff_coordinates(nodes_df, links_df)
    df = nodes_df.merge(
        links_df[GFF_KEY_COLS + ['parent_type', 'parent_id']].drop_duplicates(),
        how='outer',
        on=GFF_KEY_COLS,
    ).fillna('')

    df = enforce_uniq_transcript_ids(df)

    def feature_key(row, parent=False):
        if not parent:
            return tuple([row[c] for c in ['feature_id', 'type', 'seqid', 'strand']])
        else:
            return tuple([row[c] for c in ['parent_id', 'parent_type', 'seqid', 'strand']])

    genes_by_id = {}
    for _, row in df[df.type == 'gene'].iterrows():
        genes_by_id[feature_key(row)] = {
            'start': row.start,
            'end': row.end,
            'chr': row.seqid,
            'aliases': pull_alias_terms(row),
            'strand': row.strand,
            'transcripts': [],
            'name': row.feature_id,
            'version': row.version,
            'biotype': row.biotype,
            'note': row.Note,
        }
    logging.info(f'loaded {len(genes_by_id)} genes')

    transcripts_by_id = {}
    df = df.fillna('')

    for _, row in df[df.type == 'transcript'].iterrows():
        parent_key = feature_key(row, True)
        if parent_key not in genes_by_id:
            raise KeyError(
                f'cannot find gene ({row.parent_id}) skipping feature ({row.feature_id}) on line ({row.row_index})'
            )
        feature_id = row.feature_id
        transcript = {
            'name': feature_id,
            'start': row.start,
            'end': row.end,
            'aliases': pull_alias_terms(row),
            'domains': [],
            'exons': [],
            'version': row.version,
            'note': row.Note,
            'biotype': row.biotype,
        }
        genes_by_id[parent_key]['transcripts'].append(transcript)
        transcripts_by_id[feature_key(row)] = transcript

    # now cds
    cds_by_id = {}
    for _, row in df[df.type == 'CDS'].iterrows():
        parent_key = feature_key(row, True)
        if parent_key not in transcripts_by_id:
            print(row)
            raise KeyError(
                f'failed to find parent transcript ({row.parent_id}) skipping cds ({row.feature_id}) on line ({row.row_index})'
            )
        parent = transcripts_by_id[parent_key]
        parent.setdefault('translations', [])
        cds = {
            'start': row.start,
            'end': row.end,
            'name': row.feature_id,
            'aliases': pull_alias_terms(row),
            'version': row.version,
            'note': row.Note,
            'biotype': row.biotype,
        }
        parent['translations'].append(cds)
        cds_by_id[feature_key(row)] = cds

    logging.info(f'loaded {len(transcripts_by_id)} transcripts')
    logging.info(f'loaded {len(cds_by_id)} cds regions')
    # exons
    exons_by_id = {}

    for _, row in df[df.type == 'exon'].iterrows():
        parent_key = feature_key(row, True)
        if parent_key not in transcripts_by_id:
            raise KeyError(
                f'failed to find parent transcript ({row.parent_id}) skipping exon ({row["feature_id"]}) index={row["row_index"]}'
            )
        exon = {
            'start': row.start,
            'end': row.end,
            'name': row.feature_id,
            'version': row.version,
            'number': row.exon_number,
        }
        transcripts_by_id[parent_key]['exons'].append(exon)
        exons_by_id[feature_key(row)] = exon

    logging.info(f'loaded {len(exons_by_id)} exons')

    ignored_df = df[~df.type.isin({'exon', 'CDS', 'transcript', 'gene'})]
    if ignored_df.shape[0]:
        logging.warning(
            f'Ignored {ignored_df.shape[0]} rows that did not match the expected types: {ignored_df.type.unique()}'
        )

    result = strip_empty_fields({'genes': list(genes_by_id.values())})

    try:
        parse_annotations_json(result)
    except Exception as err:
        short_msg = '. '.join(
            [line for line in str(err).split('\n') if line.strip()][:3]
        )  # these can get super long
        raise AssertionError(short_msg)
    # re-strip (mavis adds defaults)
    result = strip_empty_fields({'genes': list(genes_by_id.values())})
    return result


def convert_gff3_to_mavis(filename: str, no_alt=False) -> Dict:
    """
    Convert an input gff3 file to the JSON format accepted by MAVIS
    """
    logging.info(f'reading: {filename}')
    df = pd.read_csv(
        filename,
        sep='\t',
        dtype={
            'seqid': str,
            'source': str,
            'type': str,
            'start': int,
            'end': int,
            'score': str,
            'strand': str,
            'phase': str,
            'attributes': str,
        },
        index_col=False,
        header=None,
        comment='#',
        na_values=['.'] + PANDAS_DEFAULT_NA_VALUES,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
    )
    df['row_index'] = df.index
    if no_alt:
        df = df[~df.seqid.str.startswith('GL')]
        df = df[~df.seqid.str.startswith('KI')]

    skip_types = {
        'five_prime_UTR',
        'five_prime_UTR',
    }
    df = df[~df.type.isin(skip_types)]

    def split_attributes(row):
        result = {}
        for attr in row.attributes.split(';'):
            name, value = attr.split('=')
            result[name] = value
        return [row.row_index] + [result.get(c, '') for c in GFF_ATTRS]

    prev_size = df.shape[0]
    attrs_df = pd.DataFrame(
        df.apply(split_attributes, axis=1).tolist(),
        columns=['row_index'] + GFF_ATTRS,
    )
    assert prev_size == attrs_df.shape[0]
    df = df.merge(attrs_df, on=['row_index'])
    df = df.drop(columns=['attributes'])

    assert prev_size == df.shape[0]

    df['feature_id'] = df.apply(parse_gff_id, axis=1)
    df.loc[(df.feature_id == '') & (df.type == 'exon'), 'feature_id'] = df.exon_id
    df = df[df.feature_id != '']
    df['strand'] = df.strand.fillna('')
    return convert_pandas_gff_to_mavis(df)


def convert_gff2_to_mavis(filename: str, no_alt=False) -> Dict:
    """
    Convert an input gff2/gtf file to the JSON format accepted by MAVIS
    """
    logging.info(f'reading: {filename}')
    df = pd.read_csv(
        filename,
        sep='\t',
        dtype={
            'seqname': str,
            'source': str,
            'feature': str,
            'start': int,
            'end': int,
            'score': str,
            'strand': str,
            'frame': str,
            'attribute': str,
        },
        index_col=False,
        header=None,
        comment='#',
        na_values=['.'] + PANDAS_DEFAULT_NA_VALUES,
        names=[
            'seqname',
            'source',
            'feature',
            'start',
            'end',
            'score',
            'strand',
            'frame',
            'attribute',
        ],
    ).rename(
        columns={'feature': 'type', 'seqname': 'seqid', 'frame': 'phase', 'attribute': 'attributes'}
    )  # match gff3 names
    df['row_index'] = df.index

    if no_alt:
        df = df[~df.seqid.str.startswith('GL')]
        df = df[~df.seqid.str.startswith('KI')]

    skip_types = {
        'five_prime_utr',
        'five_prime_utr',
    }
    df = df[~df.type.isin(skip_types)]

    def split_attributes(row):
        result = {}
        for attr in row.attributes.split('";'):
            if not attr.strip():
                continue
            m = re.match(r'^\s*([^"]+)\s+"(.*)"?$', attr)
            if not m:
                raise KeyError(f'attributes do not follow expected pattern: {attr}')
            result[m.group(1)] = m.group(2)
        return [row.row_index] + [result.get(c, '') for c in GFF_ATTRS]

    prev_size = df.shape[0]
    attrs_df = pd.DataFrame(
        df.apply(split_attributes, axis=1).tolist(),
        columns=['row_index'] + GFF_ATTRS,
    )
    assert prev_size == attrs_df.shape[0]
    df = df.merge(attrs_df, on=['row_index'])
    assert prev_size == df.shape[0]
    df = df.drop(columns=['attributes'])

    df['Alias'] = ''
    df['feature_id'] = df.apply(parse_gff_id, axis=1)

    df['Name'] = ''
    df.loc[df.type == 'gene', 'Name'] = df.gene_name
    df.loc[df.type == 'transcript', 'Name'] = df.transcript_name
    df['strand'] = df.strand.fillna('')
    df['gene_id'] = df.gene_id.astype(str)
    df.loc[df.gene_id.str.startswith('unassigned_gene_'), 'gene_id'] = ''

    df['Parent'] = ''
    df.loc[(df.type == 'transcript') & (df.gene_id != ''), 'Parent'] = 'gene:' + df.gene_id
    df.loc[(df.type == 'exon') & (df.transcript_id != ''), 'Parent'] = (
        'transcript:' + df.transcript_id
    )
    df.loc[(df.type == 'CDS') & (df.transcript_id != ''), 'Parent'] = (
        'transcript:' + df.transcript_id
    )
    df.loc[
        (df.type == 'CDS') & df.Parent.str.startswith('transcript:unassigned_transcript_'), 'Parent'
    ] = ''
    df.loc[(df.type == 'CDS') & (df.Parent == '') & (df.gene_id != ''), 'Parent'] = (
        'gene:' + df.gene_id
    )

    df['version'] = ''
    df.loc[df.type == 'transcript', 'version'] = df.transcript_version
    df.loc[df.type == 'exon', 'version'] = df.exon_version
    df.loc[df.type == 'gene', 'version'] = df.gene_version
    df.loc[df.type == 'CDS', 'version'] = df.protein_version

    df['strand'] = df.strand.fillna('')
    return convert_pandas_gff_to_mavis(df)


def convert_mavis_json_2to3(filename):
    logging.info(f'loading: {filename}')
    with open(filename, 'r') as fh:
        content = json.load(fh)

    # move translations into sep object
    skipped_tx = 0
    total_tx = 0
    for gene in content['genes']:
        if str(gene['strand']) == '1':
            gene['strand'] = '+'
        elif str(gene['strand']) == '-1':
            gene['strand'] = '-'
        for transcript in gene.get('transcripts', []):
            if all(transcript.get(k) for k in ['cdna_coding_start', 'cdna_coding_end']):
                total_tx += 1
                translation = {
                    'cdna_coding_start': transcript['cdna_coding_start'],
                    'cdna_coding_end': transcript['cdna_coding_end'],
                    'domains': transcript.get('domains', []),
                }
                translated_length = (
                    1 + transcript['cdna_coding_end'] - transcript['cdna_coding_start']
                )

                if 'domains' in transcript:
                    del transcript['domains']

                del transcript['cdna_coding_start']
                del transcript['cdna_coding_end']

                if translated_length % 3 != 0:
                    skipped_tx += 1
                    logging.debug(
                        f'Ignoring translation ({transcript.get("name")}). The translated region is not a multiple of three (length={translated_length})'
                    )
                    continue
                transcript['translations'] = [translation]
    if skipped_tx:
        logging.warning(
            f'dropped {skipped_tx} / {total_tx} translations for lengths that were not a multiple of 3'
        )
    content = coerce_number_types(content)
    content = strip_empty_fields(content)
    logging.info('testing new JSON with MAVIS loader')
    parse_annotations_json(content)
    logging.info('removing unnecessary empty fields')
    content = strip_empty_fields(content)
    return content


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input', help='path to the tab-delimated mavis v2 style reference annotations file'
    )
    parser.add_argument('--input_type', default='v2', choices=['v2-tab', 'v2-json', 'gff3', 'gtf'])
    parser.add_argument('output', help='path to the JSON output file')
    parser.add_argument(
        '--filter_alt',
        help='filter out chromosome/seqid names starting with GL or KI',
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--log_level', choices=['INFO', 'DEBUG', 'WARNING', 'ERROR'], default='INFO'
    )

    args = parser.parse_args()
    logging.basicConfig(
        format='{asctime} [{levelname}] {message}',
        style='{',
        level=logging.getLevelName(args.log_level),
    )

    if args.input_type == 'v2-tab':
        annotations = convert_tab_to_json(args.input)
    elif args.input_type == 'v2-json':
        annotations = convert_mavis_json_2to3(args.input)
    elif args.input_type == 'gtf':
        annotations = convert_gff2_to_mavis(args.input, args.filter_alt)
    else:
        annotations = convert_gff3_to_mavis(args.input, args.filter_alt)

    logging.info(f'writing: {args.output}')
    with open(args.output, 'w') as fh:
        fh.write(json.dumps(annotations, sort_keys=True))


if __name__ == '__main__':
    main()
