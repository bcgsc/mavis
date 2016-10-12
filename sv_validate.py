#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
import TSV
import argparse
import os
from structural_variant import __version__
from structural_variant.constants import *
from structural_variant.align import *
from structural_variant.error import *
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.validate import Evidence
import structural_variant.annotate as ann
from difflib import SequenceMatcher

__prog__ = os.path.basename(os.path.realpath(__file__))

MAX_POOL_SIZE = 4
HRG = '/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa'

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

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + __version__,
                        help='Outputs the version number')
parser.add_argument('-f', '--overwrite', action='store_true', default=False,
                        help='set flag to overwrite existing reviewed files')
parser.add_argument('-o', '--output', help='path to the output directory', required=True)
parser.add_argument('-n', '--input', help='path to the input file', required=True)
parser.add_argument('--max_pools', '-p', default=MAX_POOL_SIZE, type=int,
        help='defines the maximum number of processes', dest='MAX_POOL_SIZE')
parser.add_argument('-b', '--bamfile', help='path to the input bam file', required=True)
parser.add_argument('-l', '--library', help='library id', required=True)
parser.add_argument('-r', '--reference_genome', default=HRG,
        help='path to the reference genome fasta file')
args = parser.parse_args()

if args.MAX_POOL_SIZE < 1:
    print('\nerror: MAX_POOL_SIZE must be a positive integer')
    parser.print_help()
    exit(1)

MAX_POOL_SIZE = args.MAX_POOL_SIZE

header, rows = TSV.read_file(args.input,
        retain = [
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
        cast = {
            'cluster_size': 'int', 
            'break1_position_start': 'int',
            'break1_position_end': 'int',
            'break2_position_start': 'int',
            'break2_position_end': 'int',
            'opposing_strands': 'bool'
            }
        )

# load the reference genome
print('loading the reference genome', args.reference_genome)
ann.load_reference_genome(args.reference_genome)
print('loading complete')

input_bamfile = pysam.AlignmentFile(args.bamfile, 'rb')
new_bam_name = args.input + '.bam'
new_contig_bam = args.input + '.contigs.bam'
new_bamfile = pysam.AlignmentFile(args.input + '.bam', 'wb', template=input_bamfile)
evidence_reads = set()
print('writing:', new_bamfile)
output_rows = set()

split_read_contigs = set()
chr_to_index = {}

evidence_obj_list = []

for row in rows:
    bpp = BreakpointPair(
            Breakpoint(
                row['break1_chromosome'], 
                row['break1_position_start'],
                row['break1_position_end'],
                strand = row['break1_strand'],
                orient = row['break1_orientation']
                ),
            Breakpoint(
                row['break2_chromosome'], 
                row['break2_position_start'],
                row['break2_position_end'],
                strand = row['break2_strand'],
                orient = row['break2_orientation']
                ),
            opposing_strands = row['opposing_strands']
            )
    e = Evidence(bpp, bamfile = args.bamfile, labels=row)
    bpp.label = row['cluster_id']
    
    if row['protocol'] == 'genome':
        # gather evidence
        # resolve the breakpoints
        # classify the evidence resolved breakpoints
        # annotate
        # write ouputs
        print('======================================================== loading evidence for', 
                e.breakpoint_pair, BreakpointPair.classify(e.breakpoint_pair))
        e.load_evidence(input_bamfile)
        print('\ninitial pair', e.breakpoint_pair, e.linking_evidence(), e.classification,
                    'flanking reads:', [len(a) for a in e.flanking_reads.values()],
                    'split reads:', [len(a) for a in e.split_reads.values()],
                        )
        print('try assembling split reads')
        assembly_sequences = {}
        first = True
        for reads in e.split_reads.values():
            for r in reads:
                s = r.query_sequence
                if not first and e.breakpoint_pair.opposing_strands:
                    temp = Seq(s, DNA_ALPHABET)
                    temp = temp.reverse_complement()
                    s = str(temp)
                assembly_sequences[s] = assembly_sequences.get(s, set())
                assembly_sequences[s].add(r)
            first = False
        contigs = assemble(assembly_sequences)
        filtered_contigs = []
        for contig in contigs:
            if contig.remap_score() < 3:
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
        chr_to_index.update(e.settings.convert_chr_to_index)
        evidence_reads.update(e.supporting_reads())
        evidence_obj_list.append((e, contigs, assembly_sequences))
    else:
        raise NotImplementedError('currently only genome protocols are supported')

results = blat_contigs(split_read_contigs, chr_to_index = chr_to_index)

# go through the alignments for each blat result and compute a score
# based on the mate-pairs of the reads the realigned to the contigs

MIN_EXTEND_OVERLAP = 6 # on each end

bedfile_regions = set()

for e, contigs, assembly_sequences in evidence_obj_list:
    print('\nblatted', len(contigs), 'contigs for', e, e.breakpoint_pair, BreakpointPair.classify(e.breakpoint_pair))
    if e.breakpoint_pair.break1.chr == e.breakpoint_pair.break2.chr:
        bedfile_regions.add('{0}\t{1}\t{2}\t{3}'.format(
            e.breakpoint_pair.break1.chr, 
            e._window(e.breakpoint_pair.break1)[0],
            e._window(e.breakpoint_pair.break2)[1],
            str(e.breakpoint_pair) + str(BreakpointPair.classify(e.breakpoint_pair))
            ))
    else:
        bedfile_regions.add('{0}\t{1}\t{2}\t{3}'.format(
            e.breakpoint_pair.break1.chr, 
            e._window(e.breakpoint_pair.break1)[0],
            e._window(e.breakpoint_pair.break1)[1],
            str(e.breakpoint_pair) + str(BreakpointPair.classify(e.breakpoint_pair))
            ))
        bedfile_regions.add('{0}\t{1}\t{2}\t{3}'.format(
            e.breakpoint_pair.break2.chr, 
            e._window(e.breakpoint_pair.break2)[0],
            e._window(e.breakpoint_pair.break2)[1],
            str(e.breakpoint_pair) + str(BreakpointPair.classify(e.breakpoint_pair))
            ))
    # for each evidence contig 
    #   1. choose the alignment that aligns to both breakpoint windows first
    #   2. choose the pair of alignments that make up both breakpoints
    #   3. choose the best scoring alignment that lands in one of the windows
    # score low if it is not the best alignment
    aligned_contigs = []
    for contig in contigs:
        aligned_contigs.extend(results[contig.seq])
        print('\n>', contig.seq)
        for aln in sorted(results[contig.seq], key=lambda x: x.get_tag('br')):
            print('{0}:{1}'.format(aln.reference_id, aln.reference_start), aln.cigar, aln.get_tags())
    
    putative_alignments = []
    putative_alignment_pairs = []
    
    if e.break1.chr == e.break2.chr:
        for read in aligned_contigs:
            # if it covers both breakpoints add to putative alignments
            temp = Interval(read.reference_start, read.reference_end - 1)
            if e.settings.convert_index_to_chr[read.reference_id] == e.break1.chr \
                    and e.window1.overlap(temp) \
                    and e.window2.overlap(temp):
                putative_alignments.append(read)

    for a1, a2 in itertools.combinations([x for x in aligned_contigs if x not in putative_alignments], 2):
        # do they overlap both breakpoints
        # do their mate pairs make sense
        union = Interval.union(a1.query_coverage_interval(), a2.query_coverage_interval())
        if len(union) - len(a1.query_coverage_interval()) < MIN_EXTEND_OVERLAP \
                or len(union) - len(a2.query_coverage_interval()) < MIN_EXTEND_OVERLAP:
                    continue
        if e.settings.convert_index_to_chr[a1.reference_id] == e.break1.chr \
                and e.window1.overlap((a1.reference_start, a1.reference_end - 1)) \
                and e.settings.convert_index_to_chr[a2.reference_id] == e.break2.chr \
                and e.window2.overlap((a2.reference_start, a2.reference_end - 1)):
            putative_alignment_pairs.append((a1, a2))
        elif e.settings.convert_index_to_chr[a2.reference_id] == e.break1.chr \
                and e.window1.overlap((a2.reference_start, a2.reference_end - 1)) \
                and e.settings.convert_index_to_chr[a1.reference_id] == e.break2.chr \
                and e.window2.overlap((a1.reference_start, a1.reference_end - 1)):
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




temp = args.input + '.bed'
with  open(temp, 'w') as fh:
    print('writing regions to:', temp)
    for s in bedfile_regions:
        fh.write(s + '\n')

new_contig_bamfile = pysam.AlignmentFile(new_contig_bam, 'wb', template=input_bamfile)

for contig, read_set in results.items():
    # take the top 2 alignments or the top alignment
    reads = sorted(read_set, key = lambda x: x.get_tag('bs'), reverse=True) # decreasing order
    unique_scores = []
    filtered = []
    for read in reads:
        if read.get_tag('bs') not in unique_scores:
            if len(unique_scores) > 1:
                break
            else:
                unique_scores.append(read.get_tag('bs'))
        filtered.append(read)
    for read in filtered:
        read.cigar = CigarTools.convert_for_igv(read.cigar)
        new_contig_bamfile.write(read)
print('writing:', new_contig_bam)

new_contig_bamfile.close()
input_bamfile.close()

for read in evidence_reads:
    read.cigar = CigarTools.convert_for_igv(read.cigar)
    new_bamfile.write(read)

new_bamfile.close()

