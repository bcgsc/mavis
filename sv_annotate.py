
"""
annotates breakpoint pairs and draws visualizations
"""
import argparse
from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.annotate import gather_annotations, load_reference_genes
from structural_variant import __version__
import TSV
from structural_variant.constants import PROTOCOL, SVTYPE
import types
import re
import datetime


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
    parser.add_argument(
        '--no_draw', default=True, action='store_false',
        help='set flag to suppress svg drawings of putative annotations')
    parser.add_argument(
        '-o', '--output',
        help='path to the output directory', required=True
    )
    parser.add_argument(
        '-n', '--input', action='append',
        help='path to the input file', required=True
    )
    parser.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.',
    )
    args = parser.parse_args()
    return args


def main():
    # load the file
    args = parse_arguments()

    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations)

    bpps = []
    for f in args.input:
        temp = read_bpp_from_input_file(
            f,
            require=['cluster_id', 'validation_id'],
            cast={
                'stranded': TSV.bool
            },
            _in={
                'protocol': PROTOCOL,
                'event_type': SVTYPE
            },
            simplify=False)
        bpps.extend(temp)
    print('read {} breakpoint pairs'.format(len(bpps)))
    
    with open(args.output, 'w') as fh:
        annotations = []
        header = set()
        for bpp in bpps:
            ann = gather_annotations(REFERENCE_ANNOTATIONS, bpp, event_type=bpp.data['event_type'])
            annotations.extend(ann)
            header.update(ann.data.keys())
        temp = header
        header = [
            'cluster_id',
            'validation_id',
            'annotation_id',
            'break1_chromosome',
            'break1_position_start',
            'break1_position_end',
            'break1_orientation',
            'break1_strand',
            'break2_chromosome',
            'break2_position_start',
            'break2_position_end',
            'break2_orientation',
            'break2_strand',
            'opposing_strands',
            'stranded',
            'untemplated_sequence',
            'event_type',
            'transcript1',
            'transcript2',
            'genes_encompassed',
            'genes_overlapping_break1',
            'genes_overlapping_break2',
            'genes_proximal_to_break1',
            'genes_proximal_to_break2'
        ]
        for col in sorted(list(temp)):
            if col not in header:
                header.append(col)
        fh.write('\t'.join(header) + '\n')

        id_prefix = 'annotation_{}-'.format(re.sub(' ', '_', str(datetime.now())))

        for i, ann in enumerate(annotations):
            row = ann.flatten()
            row['annotation_id'] = id_prefix + str(i + 1)
            temp = []
            for col in header:
                if not isinstance(row[col], types.StringTypes):
                    try:
                        temp.append(';'.join([str(x) for x in row[col]]))
                        continue
                    except TypeError:
                        pass
                temp.append(str(row[col]))
            fh.write('\t'.join(temp) + '\n')

        print('generated {} annotations'.format(len(annotations)))


if __name__ == '__main__':
    main()
