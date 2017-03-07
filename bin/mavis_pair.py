
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
    |       |-- <library 2>_<protocol 2>.paired.tab
    |       `-- edges.tab
    `-- summary/

General Process
----------------


"""
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.breakpoint import read_bpp_from_input_file
from mavis.annotate import load_reference_genes
from mavis.pairing import equivalent_events
from mavis import __version__
from Bio import SeqIO
import TSV
from mavis.constants import PROTOCOL, SVTYPE, COLUMNS, SPLICE_TYPE, CALL_METHOD, log
import networkx as nx
import itertools


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
    parser.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )

    args = parser.parse_args()
    return args


def main():
    # load the file
    args = parse_arguments()
    DISTANCES = {
        CALL_METHOD.FLANK: args.flanking_call_distance,
        CALL_METHOD.SPLIT: args.split_call_distance,
        CALL_METHOD.CONTIG: args.contig_call_distance
    }
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
                    COLUMNS.fusion_cdna_coding_end,
                    COLUMNS.fusion_sequence_fasta_id,
                    COLUMNS.fusion_sequence_fasta_file,
                    COLUMNS.transcript1,
                    COLUMNS.transcript2
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
    libraries = set()

    SEQUENCES = dict()
    sequence_files = set()
    for bpp in bpps:
        if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
            sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
        libraries.add(bpp.data[COLUMNS.library])
    for f in sorted(list(sequence_files)):
        log('loading:', f)
        with open(f, 'rU') as fh:
            temp = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
            for k in temp:
                if k in SEQUENCES:
                    raise AssertionError('sequence identifiers are not unique', k)
            SEQUENCES.update(temp)

    log('loading:', args.annotations)
    REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations)

    TRANSCRIPTS = dict()

    for chr in REFERENCE_ANNOTATIONS:
        for gene in REFERENCE_ANNOTATIONS[chr]:
            for t in gene.transcripts:
                k = (gene.name, t.name)
                if k in TRANSCRIPTS:
                    raise AssertionError('gene + transcript is not unique', gene, t)
                TRANSCRIPTS[k] = t

    # now try comparing breakpoints between libraries
    calls_by_lib = dict()

    def pair_key(bpp):
        key = [
            bpp.data[COLUMNS.library],
            bpp.data[COLUMNS.protocol],
            bpp.data[COLUMNS.annotation_id],
            bpp.data[COLUMNS.fusion_splicing_pattern],
            bpp.data[COLUMNS.fusion_cdna_coding_start],
            bpp.data[COLUMNS.fusion_cdna_coding_end]
        ]
        return '_'.join([str(k) for k in key])

    all_bpp = dict()
    for bpp in bpps:
        key = (bpp.data[COLUMNS.library], bpp.data[COLUMNS.protocol])
        if key not in calls_by_lib:
            calls_by_lib[key] = dict()

        k = pair_key(bpp)
        if k in calls_by_lib[key]:
            raise KeyError('duplicate bpp is not unique within lib', k, bpp, bpp.data)
        calls_by_lib[key][k] = bpp
        all_bpp[k] = bpp

    pairing = nx.Graph()

    # pairwise comparison of breakpoints between all libraries
    for l1, l2 in itertools.combinations(calls_by_lib.keys(), 2):
        # for each two libraries pair all calls
        print(len(calls_by_lib[l1]) * len(calls_by_lib[l2]), 'comparisons between', l1, 'and', l2)
        for bpp1, bpp2 in itertools.product(calls_by_lib[l1], calls_by_lib[l2]):
            if equivalent_events(
                calls_by_lib[l1][bpp1],
                calls_by_lib[l2][bpp2],
                DISTANCES=DISTANCES,
                TRANSCRIPTS=TRANSCRIPTS,
                SEQUENCES=SEQUENCES
            ):
                pairing.add_edge(bpp1, bpp2)

    OUTPUT_DIR = os.path.join(args.output, '_'.join(sorted(list(libraries))))
    of = os.path.join(OUTPUT_DIR, 'edges.tab')
    with open(of, 'w') as fh:
        log('writing:', of)
        fh.write('source\ttarget\t\n')
        for src, tgt in pairing.edges():
            fh.write('{}\t{}\n'.format(src, tgt))

    for lib, protocol in calls_by_lib:
        for k, bpp in calls_by_lib[(lib, protocol)].items():
            # for each pair list all the pairings to other libraries
            paired_to = set()
            for node in nx.all_neighbors(pairing, k):
                p = all_bpp[node]
                temp = '{}_{}'.format(p.data[COLUMNS.library], p.data[COLUMNS.protocol])
                paired_to.add(temp)


if __name__ == '__main__':
    main()
