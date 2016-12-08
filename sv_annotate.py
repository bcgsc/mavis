
"""
collapses and annotates breakpoint pairs from the validation step
"""
import argparse
from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.annotate import gather_annotations, load_reference_genes
from structural_variant import __version__
import TSV
from structural_variant.constants import PROTOCOL, CALL_METHOD


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        '-f', '--overwrite',
        action='store_true', default=False,
        help='set flag to overwrite existing reviewed files'
    )
    # parser.add_argument(
    #     '-o', '--output',
    #     help='path to the output directory', required=True
    # )
    parser.add_argument(
        '-n', '--input', action='append',
        help='path to the input file', required=True
    )
    # parser.add_argument(
    #     '-l', '--library',
    #     help='library id', required=True
    # )
    parser.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.',
    )
    args = parser.parse_args()
    return args


def compare_validation_row(current, other):
    if current['call_method'] != other['call_method']:
        if current['call_method'] == CALL_METHOD.CONTIG:
            return current
        elif other['call_method'] == CALL_METHOD.CONTIG:
            return other
        elif current['call_method'] == CALL_METHOD.SPLIT:
            return current
        elif other['call_method'] == CALL_METHOD.SPLIT:
            return other
        elif current['call_method'] == CALL_METHOD.MIXED:
            return current
        elif other['call_method'] == CALL_METHOD.MIXED:
            return other

def collapse_dictionaries(dicts, columns_to_collapse=[]):
    pass

def main():
    # load the file
    args = parse_arguments()

    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations)

    bpps = []
    for f in args.input:
        temp = read_bpp_from_input_file(
            f,
            require=[
                'tools',
                'contig_sequence',
                # 'break1_homologous_sequence',
                # 'break2_homologous_sequence',
                'break1_evidence_window',
                'break2_evidence_window'
            ],
            cast={
                'cluster_id': lambda x: x.split(';'),
                'stranded': TSV.bool,
                'contigs_assembled': int,
                'contigs_aligned': int,
                'contig_remap_score': float,
                'contig_alignment_score': float,
                'flanking_reads': int,
                'median_insert_size': int,
                'stdev_insert_size': float,
                'break1_split_reads': int,
                'break2_split_reads': int,
                'linking_split_reads': int,
            },
            _in={
                'protocol': PROTOCOL,
                'call_method': CALL_METHOD,
                'event_type': SVTYPE
            })
        bpps.extend(temp)
    print('read {} breakpoint pairs'.format(len(bpps)))

    annotations = []
    for bpp in bpps:
        ann = gather_annotations(REFERENCE_ANNOTATIONS, bpp, event_type=bpp.data['event_type'])
        annotations.extend(ann)
        print(bpp)
        for a in ann:
            print('transcript1', a.transcript1)
            print('transcript2', a.transcript2)
            print('encompassed genes', a.encompassed_genes)
            print('overlap1', a.genes_at_break1)
            print('overlap2', a.genes_at_break2)
            print('nearest_gene_break1', a.nearest_gene_break1)
            print('nearest_gene_break2', a.nearest_gene_break2)

    print('generated {} annotations'.format(len(annotations)))


if __name__ == '__main__':
    main()
