import argparse
import json
import logging
import re
from typing import Dict

import pandas as pd


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


if __name__ == '__main__':
    logging.basicConfig(**{'format': '{message}', 'style': '{', 'level': logging.INFO})
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input', help='path to the tab-delimated mavis v2 style reference annotations file'
    )
    parser.add_argument('output', help='path to the JSON output file')

    args = parser.parse_args()

    annotations = convert_tab_to_json(args.input)

    rows = []
    logging.info(f'writing: {args.output}')
    if args.output_format == 'jsonl':
        with open(args.output, 'w') as fh:
            for gene in annotations['genes']:
                fh.write(json.dumps(gene, sort_keys=True) + '\n')
    elif args.output_format == 'json':
        with open(args.output, 'w') as fh:
            fh.write(json.dumps(annotations, sort_keys=True))
    else:
        transcripts = []

        for gene in annotations['genes']:
            meta = {**gene}
            del meta['transcripts']
            if gene['transcripts']:
                for transcript in gene['transcripts']:
                    transcripts.append(
                        {**meta, **{f'transcript.{k}': v for k, v in transcript.items()}}
                    )
            else:
                transcripts.append(meta)
        df = pd.json_normalize(transcripts, max_level=1)
        json_cols = [
            'aliases',
            'transcript.aliases',
            'transcript.exons',
            'transcript.domains',
        ]
        for col in json_cols:
            df[col] = df[col].apply(json.dumps)
        df.to_csv(args.output, index=False, sep='\t')
