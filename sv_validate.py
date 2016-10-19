#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import TSV
import argparse
import os
import re
from structural_variant import __version__
from structural_variant.constants import *
from structural_variant.align import *
from structural_variant.error import *
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.validate import Evidence
from structural_variant.blat import blat_contigs
from structural_variant.interval import Interval
import structural_variant.annotate as ann
from difflib import SequenceMatcher
from datetime import datetime

__prog__ = os.path.basename(os.path.realpath(__file__))

MAX_POOL_SIZE = 4
HRG = '/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'


class Profile:
    """
    this is only being used while building the application
    to get a sense of how long the different steps are taking
    """
    steps = []

    @classmethod
    def mark_step(cls, name):
        cls.steps.append((datetime.now(), name))

    @classmethod
    def print(cls):
        for step, time in cls.steps:
            print(step, time)


def compare_strings(first, second):
    return SequenceMatcher(None, first, second).ratio()


def mkdirp(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


def read_cluster_file(name, bamcache, HUMAN_REFERENCE_GENOME):
    header, rows = TSV.read_file(
        name,
        retain=[
            'cluster_id',
            'cluster_size',
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
            'protocol',
            'tools'
        ],
        cast={
            'cluster_size': 'int',
            'break1_position_start': 'int',
            'break1_position_end': 'int',
            'break2_position_start': 'int',
            'break2_position_end': 'int',
            'opposing_strands': 'bool'
        }
    )
    evidence = []
    for row in rows:
        bpp = BreakpointPair(
            Breakpoint(
                row['break1_chromosome'],
                row['break1_position_start'],
                row['break1_position_end'],
                strand=row['break1_strand'],
                orient=row['break1_orientation']
            ),
            Breakpoint(
                row['break2_chromosome'],
                row['break2_position_start'],
                row['break2_position_end'],
                strand=row['break2_strand'],
                orient=row['break2_orientation']
            ),
            opposing_strands=row['opposing_strands']
        )
        e = Evidence(bpp, bamcache, HUMAN_REFERENCE_GENOME, labels=row, protocol=row['protocol'])
        bpp.label = row['cluster_id']
        evidence.append(e)
    return evidence


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                        help='Outputs the version number')
    parser.add_argument('-f', '--overwrite', action='store_true', default=False,
                        help='set flag to overwrite existing reviewed files')
    parser.add_argument(
        '-o', '--output', help='path to the output directory', required=True)
    parser.add_argument(
        '-n', '--input', help='path to the input file', required=True)
    parser.add_argument('--max_pools', '-p', default=MAX_POOL_SIZE, type=int,
                        help='defines the maximum number of processes', dest='MAX_POOL_SIZE')
    parser.add_argument('-b', '--bamfile',
                        help='path to the input bam file', required=True)
    parser.add_argument('-l', '--library', help='library id', required=True)
    parser.add_argument('-r', '--reference_genome', default=HRG,
                        help='path to the reference genome fasta file')
    args = parser.parse_args()
    if args.MAX_POOL_SIZE < 1:
        print('\nerror: MAX_POOL_SIZE must be a positive integer')
        parser.print_help()
        exit(1)
    return args


def main():
    """
    - read the evidence
    - assemble contigs from the split reads
    - blat the contigs
    - pair the blatted contigs (where appropriate)
    - TODO: call the breakpoints and summarize the evidence
    """
    args = parse_arguments()

    FILENAME_PREFIX = re.sub('\.(txt|tsv)$', '', os.path.basename(args.input))
    EVIDENCE_BAM = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bam')
    CONTIG_BAM = os.path.join(args.output, FILENAME_PREFIX + '.contigs.bam')
    EVIDENCE_BED = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bed')
    MIN_EXTEND_OVERLAP = 6  # on each end
    MIN_CONTIG_READ_REMAP = 3
    MIN_BREAKPOINT_RESOLUTION = 3
    INPUT_BAM_CACHE = BamCache(args.bamfile)
    # load the reference genome
    print('loading the reference genome', args.reference_genome)
    Profile.mark_step('start loading reference genome')
    HUMAN_REFERENCE_GENOME = ann.load_reference_genome(args.reference_genome)
    Profile.mark_step('finished loading reference genome')
    print('loading complete')

    evidence_reads = set()

    split_read_contigs = set()
    chr_to_index = {}

    evidence = []
    Profile.mark_step('load bam reads and assemble')
    for e in read_cluster_file(args.input, INPUT_BAM_CACHE, HUMAN_REFERENCE_GENOME):
        if e.protocol == PROTOCOL.GENOME:
            # gather evidence
            # resolve the breakpoints
            # classify the evidence resolved breakpoints
            # annotate
            # write ouputs
            print('======================================================== loading evidence for',
                  e.breakpoint_pair, BreakpointPair.classify(e.breakpoint_pair))
            e.load_evidence()
            print('\ninitial pair', e.breakpoint_pair, e.linking_evidence(), e.classification,
                  'flanking reads:', len(e.flanking_reads),
                  'split reads:', [len(a) for a in e.split_reads.values()],
                  )
            print('try assembling split reads')
            assembly_sequences = {}
            for r in itertools.chain.from_iterable(e.split_reads.values()):
                s = r.query_sequence
                temp = Seq(s, DNA_ALPHABET)
                temp = str(temp.reverse_complement())
                assembly_sequences[s] = assembly_sequences.get(s, set())
                assembly_sequences[s].add(r)
                assembly_sequences[temp] = assembly_sequences.get(temp, set())
                assembly_sequences[temp].add(r)
                for m in INPUT_BAM_CACHE.get_mate(r):
                    s = m.query_sequence
                    temp = Seq(s, DNA_ALPHABET)
                    temp = str(temp.reverse_complement())
                    assembly_sequences[s] = assembly_sequences.get(s, set())
                    assembly_sequences[s].add(m)
                    assembly_sequences[temp] = assembly_sequences.get(temp, set())
                    assembly_sequences[temp].add(m)

            contigs = assemble(assembly_sequences)
            filtered_contigs = []
            for contig in contigs:
                if contig.remap_score() < MIN_CONTIG_READ_REMAP:
                    continue
                print('contig{', 'score =', contig.score,
                      'adjusted-score =', round(contig.score / len(contig.seq), 2),
                      'read-remap =', contig.remap_score(),
                      'read-total =', len(contig.remapped_reads.keys()),
                      'seq =', contig.seq,
                      '}')
                filtered_contigs.append(contig)
            contigs = filtered_contigs if len(filtered_contigs) > 0 else contigs
            split_read_contigs.update(set([c.seq for c in contigs]))
            evidence_reads.update(e.supporting_reads())
            evidence.append((e, contigs, assembly_sequences))
        else:
            raise NotImplementedError('currently only genome protocols are supported')
    Profile.mark_step('start blat')
    blat_contig_alignments = blat_contigs(
        split_read_contigs,
        INPUT_BAM_CACHE
    )
    Profile.mark_step('blat complete')
    print('blatted sequences')
    for seq in blat_contig_alignments:
        print('>', seq)
    print()
    # go through the alignments for each blat result and compute a score
    # based on the mate-pairs of the reads the realigned to the contigs

    final_contigs = []
    output_rows = []

    for e, contigs, assembly_sequences in evidence:
        print('\n===================== blatted', len(contigs), 'contigs for', e)
        for c in contigs:
            print('>', c.seq)
        print('\ninitial pair', e.breakpoint_pair, e.linking_evidence(), e.classification,
                  'flanking reads:', len(e.flanking_reads),
                  'split reads:', [len(a) for a in e.split_reads.values()],
                  )
        # TODO score low if it is not the best alignment
        putative_alignments = []
        putative_alignment_pairs = []

        # gather a list of the contigs alignments that apply to this evidence object/breakpoint pair
        for contig in contigs:
            if e.break1.chr == e.break2.chr:
                for read in blat_contig_alignments[contig.seq]:
                    # if it covers both breakpoints add to putative alignments
                    temp = Interval(read.reference_start, read.reference_end - 1)
                    if INPUT_BAM_CACHE.chr(read) == e.break1.chr \
                            and e.window1.overlaps(temp) \
                            and e.window2.overlaps(temp):
                        putative_alignments.append(read)

            for a1, a2 in itertools.combinations(
                    [x for x in blat_contig_alignments[contig.seq] if x not in putative_alignments], 2):
                # do they overlap both breakpoints
                # do their mate pairs make sense
                union = Interval.union(a1.query_coverage_interval(),
                                       a2.query_coverage_interval())
                if len(union) - len(a1.query_coverage_interval()) < MIN_EXTEND_OVERLAP \
                        or len(union) - len(a2.query_coverage_interval()) < MIN_EXTEND_OVERLAP:
                    continue
                if INPUT_BAM_CACHE.chr(a1) == e.break1.chr \
                        and e.window1.overlaps((a1.reference_start, a1.reference_end - 1)) \
                        and INPUT_BAM_CACHE.chr(a2) == e.break2.chr \
                        and e.window2.overlaps((a2.reference_start, a2.reference_end - 1)):
                    putative_alignment_pairs.append((a1, a2))
                elif INPUT_BAM_CACHE.chr(a2) == e.break1.chr \
                        and e.window1.overlaps((a2.reference_start, a2.reference_end - 1)) \
                        and INPUT_BAM_CACHE.chr(a1) == e.break2.chr \
                        and e.window2.overlaps((a1.reference_start, a1.reference_end - 1)):
                    putative_alignment_pairs.append((a2, a1))

        print('single alignments', len(putative_alignments))
        for a in putative_alignments:
            print(a)
        print('paired alignments', len(putative_alignment_pairs))
        for a1, a2 in sorted(putative_alignment_pairs, key=lambda x: x[0].blat_score() + x[1].blat_score()):
            print('pair score', a1.blat_score() + a2.blat_score())
            print(a1)
            print(a2)
        print()
        final_contigs.extend(putative_alignments)
        final_contigs.extend([x for x, y in putative_alignment_pairs])
        final_contigs.extend([y for x, y in putative_alignment_pairs])

        row = {
            'cluster_id': e.labels['cluster_id'],
            'cluster_size': e.labels['cluster_size'],
            'event_type': '?',
            'protocol': e.protocol,
            'tools': e.labels['tools'],
            'contigs_assembled': len(contigs),
            'contigs_aligned': len(blat_contig_alignments[contig.seq]),
            'contigs_accepted': len(putative_alignments) + len(putative_alignment_pairs)
        }
        if len(putative_alignments) + len(putative_alignment_pairs) > 0:  # contigs aligned
            print('METHOD: contigs')
            continue
        print('METHOD: resolve by softclipped reads')
        if len(blat_contig_alignments[contig.seq]) > 0:
            print('had', len(blat_contig_alignments[contig.seq]), 'aligned contigs but none survived pairing')
            for c in sorted(blat_contig_alignments[contig.seq], key=lambda x: x.get_tag('br'))[0:5]:
                print(repr(c))
                print(c)

        else:
            print('no contigs assembled')
        for temp in e.call_by_softclipping_resolution():
            print(
                temp, temp.breakpoint_pair, 'splitreads',
                sum([len(x) for x in temp.split_reads.values()]),
                'linking reads:',
                temp.linking_evidence())
        print(row)

    # write the bam file for guiding the user in debugging why a particular cluster may not
    # pass validation
    # with open(EVIDENCE_BED, 'w') as fh:
    #     print('writing:', EVIDENCE_BED)
    #     for s in bedfile_regions:
    #         fh.write(s + '\n')

    with pysam.AlignmentFile(CONTIG_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        print('writing:', CONTIG_BAM)
        for read in final_contigs:
            read.cigar = CigarTools.convert_for_igv(read.cigar)
            fh.write(read)

    # write the evidence
    with pysam.AlignmentFile(EVIDENCE_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        print('writing:', EVIDENCE_BAM)
        for read in evidence_reads:
            read.cigar = CigarTools.convert_for_igv(read.cigar)
            fh.write(read)

    INPUT_BAM_CACHE.close()

    # write the output validated clusters (split by type and contig)
    header = [
        'cluster_id',
        'cluster_size',
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
        'event_type',
        'opposing_strands',
        'protocol',
        'tools',
        'flanking_reads',
        'spanning_reads',
        'split_reads',
        'contigs_assembled',
        'contig_sequence',
        'contigs_aligned',
        'call_method'  # contig alignment, split-read consensus, flanking read interval, input
    ]
    print(header)

try:
    main()
finally:
    print()
    Profile.print()
