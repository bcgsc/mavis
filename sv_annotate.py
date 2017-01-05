
"""
annotates breakpoint pairs and draws visualizations
"""
import argparse
from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.annotate import gather_annotations, load_reference_genes
from structural_variant import __version__
import TSV
from structural_variant.constants import PROTOCOL, SVTYPE, COLUMNS, sort_columns
import types
import re
import datetime


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


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
            require=[COLUMNS.cluster_id.name, COLUMNS.validation_id.name],
            cast={
                COLUMNS.stranded.name: TSV.bool
            },
            _in={
                COLUMNS.protocol.name: PROTOCOL,
                COLUMNS.event_type.name: SVTYPE
            },
            simplify=False)
        bpps.extend(temp)
    log('read {} breakpoint pairs'.format(len(bpps)))
    
    with open(args.output, 'w') as fh:
        annotations = []
        header = set()
        for bpp in bpps:
            ann = gather_annotations(REFERENCE_ANNOTATIONS, bpp, event_type=bpp.data[COLUMNS.event_type.name])
            annotations.extend(ann)
            header.update(ann.data.keys())
        header = sort_columns(header)
        fh.write('\t'.join(header) + '\n')

        id_prefix = 'annotation_{}-'.format(re.sub(' ', '_', str(datetime.now())))

        for i, ann in enumerate(annotations):
            row = ann.flatten()
            row[COLUMNS.annotation_id.name] = id_prefix + str(i + 1)
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

        log('generated {} annotations'.format(len(annotations)))


if __name__ == '__main__':
    main()
