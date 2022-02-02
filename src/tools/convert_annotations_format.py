import argparse
import json
import logging
import re
from typing import Dict

import pandas as pd
import pkg_resources
from snakemake.utils import validate as snakemake_validate

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
    | genomic_exon_ranges   | 100-201;334-412;779-830   | semi-colon demitited exon start/ends                      |
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

    df = pd.read_csv(
        filepath,
        dtype={
            'ensembl_gene_id': str,
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
        comment='#',
    )

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
            'domains': row.get('AA_domain_ranges', []),
            'start': row.get('transcript_genomic_start'),
            'end': row.get('transcript_genomic_end'),
            'cdna_coding_start': row.get('cdna_coding_start'),
            'cdna_coding_end': row.get('cdna_coding_end'),
            'aliases': [],
        }
        for int_value in ['start', 'end', 'cdna_coding_start', 'cdna_coding_end']:
            if transcript.get(int_value) is not None:
                transcript[int_value] = int(transcript[int_value])
        gene['transcripts'].append(transcript)

    return {'genes': list(genes.values())}


def convert_pandas_gff_to_mavis(df) -> Dict:
    df['parent_type'] = df.Parent.str.split(':').str[0]
    genelike_features = {'gene', 'ncRNA_gene', 'biological_region', 'pseudogene'}
    consumed = set()

    def pull_alias_terms(row):
        aliases = []
        if row['Name']:
            aliases.append(row['Name'])
        if row['Alias']:
            aliases.extend(row['Alias'].split(','))
        return aliases

    genes_by_id = {}
    for row in df[df.type.isin(genelike_features)].to_dict('records'):
        genes_by_id[row['feature_id']] = {
            'start': row['start'],
            'end': row['end'],
            'chr': row['seqid'],
            'aliases': pull_alias_terms(row),
            'strand': row['strand'],
            'transcripts': [],
            'name': row['feature_id'] + '.' + row['version'],
        }
        consumed.add(row['row_index'])
    logging.info(f'loaded {len(genes_by_id)} genes')

    transcripts_by_id = {}

    for row in df[df.parent_type == 'gene'].to_dict('records'):
        for parent in row['Parent'].split(','):
            gene_id = parent.split(':')[1]
            if gene_id not in genes_by_id:
                raise KeyError(
                    f'cannot find gene ({gene_id}) skipping transcript ({row["feature_id"]})'
                )
            feature_id = row['feature_id']
            transcript = {
                'name': feature_id + '.' + row['version'],
                'start': row['start'],
                'end': row['end'],
                'aliases': pull_alias_terms(row),
                'domains': [],
                'exons': [],
                'cdna_coding_start': None,
                'cdna_coding_end': None,
            }
            genes_by_id[gene_id]['transcripts'].append(transcript)
            transcripts_by_id[feature_id] = transcript
            consumed.add(row['row_index'])

    logging.info(f'loaded {len(transcripts_by_id)} transcripts')
    # now cds
    cds_count = 0
    for row in df[df.type == 'CDS'].to_dict('records'):
        for parent in row['Parent'].split(','):
            transcript_id = parent.split(':')[1]
            if transcript_id not in transcripts_by_id:
                raise KeyError(
                    f'failed to find parent transcript ({transcript_id}) skipping cds on line ({row["row_index"] + 1})'
                )
            transcripts_by_id[transcript_id].update(
                {'cdna_coding_start': row['start'], 'cdna_coding_end': row['end']}
            )
            cds_count += 1
            consumed.add(row['row_index'])
    logging.info(f'loaded {cds_count} cds regions')
    # exons
    exons_count = 0
    for row in df[df.type == 'exon'].to_dict('records'):
        for parent in row['Parent'].split(','):
            transcript_id = parent.split(':')[1]
            if transcript_id not in transcripts_by_id:
                raise KeyError(
                    f'failed to find parent transcript ({transcript_id}) skipping exon ({row["feature_id"]}) on line {row["row_index"] + 1}'
                )
            transcripts_by_id[transcript_id]['exons'].append(
                {
                    'start': row['start'],
                    'end': row['end'],
                    'name': row['feature_id'] + '.' + row['version'],
                }
            )
            exons_count += 1
            consumed.add(row['row_index'])

    logging.info(f'loaded {exons_count} exons')

    ignored_df = df[~df.row_index.isin(consumed)]
    if ignored_df.shape[0]:
        logging.warning(
            f'Ignored {ignored_df.shape[0]} rows that did not match the expected types: {ignored_df.type.unique()}'
        )

    result = {'genes': list(genes_by_id.values())}
    try:
        snakemake_validate(
            result, pkg_resources.resource_filename('mavis.annotate', 'annotations_schema.json')
        )
    except Exception as err:
        short_msg = '. '.join(
            [line for line in str(err).split('\n') if line.strip()][:3]
        )  # these can get super long
        raise AssertionError(short_msg)
    return result


def convert_gff3_to_mavis(filename: str, no_alt) -> Dict:
    """
    Convert an input gff3 file to the JSON format accepted by MAVIS
    """
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

    attribute_columns = [
        'ID',
        'Name',
        'Alias',
        'Parent',
        'Target',
        'Gap',
        'Derives_from',
        'Note',
        'DBxref',
        'Ontology_term',
        'rank',
        'version',
        'exon_id',
    ]

    def split_attributes(row):
        result = {}
        for attr in row.attributes.split(';'):
            name, value = attr.split('=')
            result[name] = value
        return [row.row_index] + [result.get(c, '') for c in attribute_columns]

    prev_size = df.shape[0]
    attrs_df = pd.DataFrame(
        df.apply(split_attributes, axis=1).tolist(),
        columns=['row_index'] + attribute_columns,
    )
    assert prev_size == attrs_df.shape[0]
    df = df.merge(attrs_df, on=['row_index'])

    assert prev_size == df.shape[0]

    df['feature_id'] = df['ID'].apply(lambda id: id.split(':')[1] if ':' in id else '')
    df.loc[(df.feature_id == '') & (df.type == 'exon'), 'feature_id'] = df.exon_id
    df = df[df.feature_id != '']
    df['strand'] = df.strand.fillna('?')
    return convert_pandas_gff_to_mavis(df)


def convert_gff2_to_mavis(filename: str, no_alt) -> Dict:
    """
    Convert an input gff2/gtf file to the JSON format accepted by MAVIS
    """
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

    attribute_columns = [
        'gene_id',
        'gene_version',
        'gene_name',
        'transcript_id',
        'transcript_version',
        'transcript_name',
        'exon_id',
        'exon_version',
    ]

    def split_attributes(row):
        result = {}
        for attr in row.attributes.split('";'):
            if not attr:
                continue
            m = re.match(r'^\s*([^"]+)\s+"(.*)"?$', attr)
            if not m:
                raise KeyError(f'attributes do not follow expected pattern: {attr}')
            result[m.group(1)] = m.group(2)
        return [row.row_index] + [result.get(c, '') for c in attribute_columns]

    prev_size = df.shape[0]
    attrs_df = pd.DataFrame(
        df.apply(split_attributes, axis=1).tolist(),
        columns=['row_index'] + attribute_columns,
    )
    assert prev_size == attrs_df.shape[0]
    df = df.merge(attrs_df, on=['row_index'])
    assert prev_size == df.shape[0]

    df['Alias'] = ''
    df['feature_id'] = ''
    df.loc[df.type == 'exon', 'feature_id'] = df.exon_id
    df.loc[df.type == 'gene', 'feature_id'] = df.gene_id
    df.loc[df.type == 'transcript', 'feature_id'] = df.transcript_id

    df['Name'] = ''
    df.loc[df.type == 'gene', 'Name'] = df.gene_name
    df.loc[df.type == 'transcript', 'Name'] = df.transcript_name
    df['strand'] = df.strand.fillna('?')

    df['Parent'] = ''
    df.loc[(df.type == 'transcript') & (df.gene_id != ''), 'Parent'] = 'gene:' + df.gene_id
    df.loc[(df.type == 'exon') & (df.transcript_id != ''), 'Parent'] = (
        'transcript:' + df.transcript_id
    )
    df.loc[(df.type == 'CDS') & (df.transcript_id != ''), 'Parent'] = (
        'transcript:' + df.transcript_id
    )

    df['version'] = ''
    df.loc[df.type == 'transcript', 'version'] = df.transcript_version
    df.loc[df.type == 'exon', 'version'] = df.exon_version
    df.loc[df.type == 'gene', 'version'] = df.gene_version

    df['strand'] = df.strand.fillna('?')
    return convert_pandas_gff_to_mavis(df)


if __name__ == '__main__':
    logging.basicConfig(format='{message}', style='{', level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input', help='path to the tab-delimated mavis v2 style reference annotations file'
    )
    parser.add_argument('--input_type', default='v2', choices=['v2', 'gff3', 'gtf'])
    parser.add_argument('output', help='path to the JSON output file')
    parser.add_argument(
        '--keep_alt',
        help='do not filter out chromosome/seqid names starting with GL or KI',
        action='store_true',
        default=False,
    )

    args = parser.parse_args()

    if args.input_type == 'v2':
        annotations = convert_tab_to_json(args.input)
    elif args.input_type == 'gtf':
        annotations = convert_gff2_to_mavis(args.input, not args.keep_alt)
    else:
        annotations = convert_gff3_to_mavis(args.input, not args.keep_alt)

    logging.info(f'writing: {args.output}')
    with open(args.output, 'w') as fh:
        fh.write(json.dumps(annotations, sort_keys=True))
