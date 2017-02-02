
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
from structural_variant.annotate import load_reference_genes
from structural_variant.interval import Interval
from structural_variant.annotate.variant import predict_transcriptome_breakpoint
from structural_variant import __version__
from Bio import SeqIO
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


def equivalent_events(ev1, ev2, TRANSCRIPTS, DISTANCES=None, SEQUENCES=None):
    if DISTANCES is None:
        DISTANCES = {CALL_METHOD.CONTIG: 0, CALL_METHOD.SPLIT: 10, CALL_METHOD.FLANK: 0}
    SEQUENCES = dict() if SEQUENCES is None else SEQUENCES
    
    # basic checks
    if ev1.break1.chr != ev2.break1.chr or ev1.break2.chr != ev2.break2.chr or \
            len(set([STRAND.NS, ev1.break1.strand, ev2.break1.strand])) > 2 or \
            len(set([STRAND.NS, ev1.break2.strand, ev2.break2.strand])) > 2 or \
            len(set([ORIENT.NS, ev1.break1.orient, ev2.break1.orient])) > 2 or \
            len(set([ORIENT.NS, ev1.break2.orient, ev2.break2.orient])) > 2 or \
            ev1.opposing_strands != ev2.opposing_strands:
        print('basic diff')
        return False
    
    methods = set([
        ev1.data[COLUMNS.break1_call_method],
        ev1.data[COLUMNS.break2_call_method],
        ev2.data[COLUMNS.break1_call_method],
        ev2.data[COLUMNS.break2_call_method]
    ])

    call_method = CALL_METHOD.CONTIG

    if CALL_METHOD.FLANK in methods:  # lowest level
        call_method = CALL_METHOD.FLANK
    elif CALL_METHOD.SPLIT in methods:
        call_method = CALL_METHOD.SPLIT
    else:  # highest level of confidence
        assert({CALL_METHOD.CONTIG} == methods)
    
    max_distance = DISTANCES[call_method]
    
    fusion1 = SEQUENCES.get(ev1.data[COLUMNS.fusion_sequence_fasta_id], None)
    fusion2 = SEQUENCES.get(ev2.data[COLUMNS.fusion_sequence_fasta_id], None)
    
    break1_match = False
    break2_match = False

    if ev1.data[COLUMNS.protocol] != ev2.data[COLUMNS.protocol]:  # mixed
        if fusion1 and fusion2:
            # compare product
            if fusion1 != fusion2:
                return False
            for col in [COLUMNS.fusion_cdna_coding_start, COLUMNS.fusion_cdna_coding_end]:
                if ev1.data[col] != ev2.data[col]:
                    return False
            return True

        # predict genome breakpoints to compare by location
        if ev1.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            t1 = TRANSCRIPTS.get(ev1.data[COLUMNS.transcript1], None)
            if t1:
                pbreaks = predict_transcriptome_breakpoint(ev1.break1, t1)
                for b in pbreaks:
                    if Interval.dist(b, ev2.break1) <= max_distance:
                        break1_match = True
                        break
            t2 = TRANSCRIPTS.get(ev1.data[COLUMNS.transcript2], None)
            if t2:
                pbreaks = predict_transcriptome_breakpoint(ev1.break2, t1)
                for b in pbreaks:
                    if Interval.dist(b, ev2.break2) <= max_distance:
                        break2_match = True
                        break
        else:
            t1 = TRANSCRIPTS.get(ev2.data[COLUMNS.transcript1], None)
            if t1:
                pbreaks = predict_transcriptome_breakpoint(ev2.break1, t1)
                for b in pbreaks:
                    if Interval.dist(b, ev1.break1) <= max_distance:
                        break1_match = True
                        break
            t2 = TRANSCRIPTS.get(ev2.data[COLUMNS.transcript2], None)
            if t2:
                pbreaks = predict_transcriptome_breakpoint(ev2.break2, t1)
                for b in pbreaks:
                    if Interval.dist(b, ev1.break2) <= max_distance:
                        break2_match = True
                        break
    elif ev1.data[COLUMNS.event_type] != ev2.data[COLUMNS.event_type]:
        print('diff events')
        return False

    # location comparison
    if Interval.dist(ev1.break1, ev2.break1) <= max_distance:
        break1_match = True
    print('break1 dist:', Interval.dist(ev1.break1, ev2.break1), ev1.break1, ev2.break1)
    if Interval.dist(ev1.break2, ev2.break2) <= max_distance:
        break2_match = True
    print('break2 dist:', Interval.dist(ev1.break2, ev2.break2), ev1.break2, ev2.break2)
    print('break1_match', break1_match, 'break2_match', break2_match)
    return break1_match and break2_match


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
                    COLUMNS.fusion_cdna_coding_end,
                    COLUMNS.fusion_sequence_fasta_id,
                    COLUMNS.fusion_sequence_fasta_file
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

    SEQUENCES = dict()
    sequence_files = set()
    for bpp in bpps:
        if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
            sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
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



if __name__ == '__main__':
    main()
