
"""
About
------

This is the fourth and final step in the svmerge pipeline. It is responsible for pairing calls between libraries
::

    <output_dir_name>/
    |-- clustering/
    |-- validation/
    |-- annotation/
    |-- pairing/
    |   `--<library 1>_<protocol 1>_and_<library 2>_<protocol 2>/
    |       |-- <library 1>_<protocol 1>.paired.tab
    |       `-- <library 2>_<protocol 2>.paired.tab
    `-- summary/

General Process
----------------


"""
import argparse
from structural_variant.breakpoint import read_bpp_from_input_file
from structural_variant.annotate import load_reference_genes, load_reference_genome, load_templates
from structural_variant.interval import Interval
from structural_variant import __version__
import TSV
from structural_variant.constants import PROTOCOL, SVTYPE, COLUMNS, SPLICE_TYPE, CALL_METHOD, STRAND, ORIENT
from datetime import datetime


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s version ' + __version__,
        help='Outputs the version number'
    )
    parser.add_argument(
        '-o', '--output',
        help='path to the output directory', required=True
    )
    parser.add_argument(
        '-n', '--inputs',
        help='path to the input file(s)', required=True, nargs='+'
    )
    parser.add_argument(
        '-f', '--fasta', help='path to the sequence files', required=True, nargs='+'
    )

    parser.add_argument(
        '-s', '--split_call_distance', default=10, type=int,
        help='distance allowed between breakpoint calls when pairing from split read (and higher) resolution calls'
    )
    parser.add_argument(
        '-c', '--contig_call_distance', default=0, type=int,
        help='distance allowed between breakpoint calls when pairing from contig (and higher) resolution calls'
    )
    parser.add_argument(
        '--flanking_call_distance', default=0, type=int,
        help='distance allowed between breakpoint calls when pairing from contig (and higher) resolution calls'
    )

    args = parser.parse_args()
    return args


def compare_trans_events(ev1, ev2):
    """
    checks if two events, at least one of which is transcriptomic, produce the same product
    """
    pass


def compare_genome_events(ev1, ev2, contig_distance=0, split_distance=10, flank_distance=0):
    """
    checks if 2 events are equivalent. Uses genomic breakpoints and call method
    """
    if ev1.break1.chr != ev2.break1.chr or ev1.break2.chr != ev2.break2.chr or \
            ev1.data[COLUMNS.event_type] != ev2.data[COLUMNS.event_type] or \
            len(set([STRAND.NS, ev1.break1.strand, ev2.break1.strand])) > 2 or \
            len(set([STRAND.NS, ev1.break2.strand, ev2.break2.strand])) > 2 or \
            len(set([ORIENT.NS, ev1.break1.orient, ev2.break1.orient])) > 2 or \
            len(set([ORIENT.NS, ev1.break2.orient, ev2.break2.orient])) > 2 or \
            ev1.opposing_strands != ev2.opposing_strands:
        return False

    methods = set([
        ev1.data[COLUMNS.break1_call_method],
        ev1.data[COLUMNS.break2_call_method],
        ev2.data[COLUMNS.break1_call_method],
        ev2.data[COLUMNS.break2_call_method]
    ])

    d = abs(Interval.dist(ev1.break1, ev2.break1)) + abs(Interval.dist(ev1.break2, ev2.break2))
    ref_d = contig_distance

    if CALL_METHOD.FLANK in methods:  # lowest level
        ref_d = flank_distance
    elif CALL_METHOD.SPLIT in methods:
        ref_d = split_distance
    else:  # highest level of confidence
        assert({CALL_METHOD.CONTIG} == methods)

    if d > ref_d * 2:
        return False
    return True


def main():
    # load the file
    args = parse_arguments()

    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)

    bpps = []
    for f in args.input:
        log('loading:', f)
        bpps.extend(
            read_bpp_from_input_file(
                f,
                require=[
                    COLUMNS.cluster_id,
                    COLUMNS.validation_id,
                    COLUMNS.annotation_id,
                    COLUMNS.library,
                    COLUMNS.fusion_cdna_coding_start,
                    COLUMNS.fusion_cdna_coding_end
                ],
                cast={
                    COLUMNS.stranded.name: TSV.tsv_boolean
                },
                _in={
                    COLUMNS.protocol: PROTOCOL,
                    COLUMNS.event_type: SVTYPE,
                    COLUMNS.fusion_splicing_pattern: SPLICE_TYPE,
                    COLUMNS.break1_call_method: CALL_METHOD,
                    COLUMNS.break2_call_method: CALL_METHOD
                },
                simplify=False
            ))
    log('read {} breakpoint pairs'.format(len(bpps)))

    log('loading:', args.reference_genome)
    REFERENCE_GENOME = load_reference_genome(args.reference_genome)

    log('loading:', args.template_metadata)
    TEMPLATES = load_templates(args.template_metadata)

    log('loading:', args.annotations)
    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations, REFERENCE_GENOME=REFERENCE_GENOME)


if __name__ == '__main__':
    main()
