#!/projects/tumour_char/analysis_scripts/python/centos06/anaconda3_v2.3.0/bin/python
"""
About
---------

This is the second step in the svmerge pipeline. This is the step responsible for validating
the events/clusters from the first/merge step. The putative breakpoint pairs are investigated
in the respective bam file. The evidence is collected and summarized. Outputs are written to
the validation subfolder in the pattern as follows

::

    <output_dir_name>/
    |-- clustering/
    |-- validation/
    |   `-- <library>_<protocol>/
    |       |-- qsub.sh
    |       |-- log/
    |       |-- clusterset-#.igv.batch
    |       |-- clusterset-#.validation-failed.tab
    |       |-- clusterset-#.validation-passed.tab
    |       |-- clusterset-#.validation-passed.bed
    |       |-- clusterset-#.contigs.bam
    |       |-- clusterset-#.contigs.tab
    |       |-- clusterset-#.contigs.sorted.bam
    |       |-- clusterset-#.contigs.sorted.bam.bai
    |       |-- clusterset-#.evidence.bam
    |       |-- clusterset-#.evidence.sorted.bam
    |       `-- clusterset-#.evidence.sorted.bam.bai
    |-- annotation/
    |-- pairing/
    `-- summary/
"""
import subprocess
import TSV
import argparse
import os
import re
from structural_variant import __version__
from structural_variant.constants import *
from structural_variant.error import *
from structural_variant.read_tools import CigarTools
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.read_tools import BamCache
from structural_variant.validate import Evidence, EvidenceSettings
from structural_variant.blat import blat_contigs
from structural_variant.interval import Interval
from structural_variant.annotate import load_masking_regions, load_reference_genome, load_reference_genes
import itertools
from datetime import datetime
import pysam

try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser

__prog__ = os.path.basename(os.path.realpath(__file__))

INPUT_BAM_CACHE = None
REFERENCE_ANNOTATIONS = None
HUMAN_REFERENCE_GENOME = None
MASKED_REGIONS = None
EVIDENCE_SETTINGS = None


def log(*pos, time_stamp=True):
    if time_stamp:
        print('[{}]'.format(datetime.now()), *pos)
    else:
        print(' ' * 28, *pos)


def mkdirp(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise


def read_cluster_file(name, is_stranded):
    header, rows = TSV.read_file(
        name,
        require=[
            COLUMNS.cluster_id,
            COLUMNS.cluster_size,
            COLUMNS.break1_chromosome,
            COLUMNS.break1_position_start,
            COLUMNS.break1_position_end,
            COLUMNS.break1_orientation,
            COLUMNS.break1_strand,
            COLUMNS.break2_chromosome,
            COLUMNS.break2_position_start,
            COLUMNS.break2_position_end,
            COLUMNS.break2_orientation,
            COLUMNS.break2_strand,
            COLUMNS.opposing_strands
        ],
        cast={
            COLUMNS.cluster_size: int,
            COLUMNS.break1_position_start: int,
            COLUMNS.break1_position_end: int,
            COLUMNS.break2_position_start: int,
            COLUMNS.break2_position_end: int,
            COLUMNS.opposing_strands: TSV.tsv_boolean
        }
    )
    evidence = []
    for row in rows:
        strands = [(row[COLUMNS.break1_strand], row[COLUMNS.break2_strand])]
        if is_stranded:
            strands = itertools.product(
                STRAND.expand(row[COLUMNS.break1_strand]),
                STRAND.expand(row[COLUMNS.break2_strand])
            )
            strands = [(s1, s2) for s1, s2 in strands if row[COLUMNS.opposing_strands] == (s1 != s2)]
        if len(strands) == 0:
            raise UserWarning('error in reading input file. could not resolve strands', row)

        for s1, s2 in strands:
            bpp = BreakpointPair(
                Breakpoint(
                    row[COLUMNS.break1_chromosome],
                    row[COLUMNS.break1_position_start],
                    row[COLUMNS.break1_position_end],
                    strand=s1,
                    orient=row[COLUMNS.break1_orientation]
                ),
                Breakpoint(
                    row[COLUMNS.break2_chromosome],
                    row[COLUMNS.break2_position_start],
                    row[COLUMNS.break2_position_end],
                    strand=s2,
                    orient=row[COLUMNS.break2_orientation]
                ),
                opposing_strands=row[COLUMNS.opposing_strands]
            )
            try:
                e = Evidence(
                    bpp,
                    INPUT_BAM_CACHE,
                    HUMAN_REFERENCE_GENOME,
                    annotations=REFERENCE_ANNOTATIONS,
                    protocol=row[COLUMNS.protocol],
                    data=row
                )
                evidence.append(e)
            except UserWarning as e:
                warnings.warn('failed to read cluster {}'.format(repr(e)))
    return evidence

def add_evidence_args_to_parser(parse_group):
    default = EvidenceSettings()
    parse_group.add_argument(
        '--read_length', default=default.read_length, type=int, help='length of reads in the bam file')
    parse_group.add_argument(
        '--stdev_isize', default=default.stdev_isize, type=int,
        help='the standard deviation in insert sizes of paired end reads')
    parse_group.add_argument(
        '--median_insert_size', default=default.median_insert_size, type=int,
        help='the median insert size of paired end reads')
    parse_group.add_argument(
        '--stdev_count_abnormal', default=default.stdev_count_abnormal, type=int,
        help='the number of standard deviations away from the normal considered expected and therefore not '
        'qualifying as flanking reads')
    parse_group.add_argument(
        '--call_error', default=default.call_error, type=int, help='buffer zone for the evidence window')
    parse_group.add_argument(
        '--min_splits_reads_resolution', default=default.min_splits_reads_resolution, type=int,
        help='minimum number of split reads required to call a breakpoint by split reads')
    parse_group.add_argument(
        '--min_anchor_exact', default=default.min_anchor_exact, type=int)
    parse_group.add_argument(
        '--min_anchor_fuzzy', default=default.min_anchor_fuzzy, type=int)
    parse_group.add_argument(
        '--min_anchor_match', default=default.min_anchor_match, type=float)
    parse_group.add_argument(
        '--min_mapping_quality', default=default.min_mapping_quality, type=int,
        help='the minimum mapping quality of reads to be used as evidence')
    parse_group.add_argument(
        '--fetch_reads_limit', default=default.fetch_reads_limit, type=int,
        help='maximum number of reads, cap, to loop over for any given evidence window')
    parse_group.add_argument(
        '--fetch_reads_bins', default=default.fetch_reads_bins, type=int,
        help='number of bins to split an evidence window into to ensure more even sampling of high coverage '
        'regions')
    parse_group.add_argument(
        '--filter_secondary_alignments', default=default.filter_secondary_alignments, type=bool,
        help='filter secondary alignments when gathering read evidence')
    parse_group.add_argument(
        '--consensus_req', default=default.consensus_req, type=int)
    parse_group.add_argument(
        '--sc_extension_stop', default=default.sc_extension_stop, type=int)
    parse_group.add_argument(
        '--assembly_min_edge_weight', default=default.assembly_min_edge_weight, type=int,
        help='when building the initial deBruijn graph edge weights are determined by the frequency of the kmer '
        'they represent in all the input sequences. The parameter here discards edges to simply the graph if they '
        'have a weight less than specified')
    parse_group.add_argument(
        '--assembly_min_remap', default=default.assembly_min_remap, type=int,
        help='The minimum input sequences that must remap for an assembly to be used')
    parse_group.add_argument(
        '--min_non_target_aligned_split_reads', default=default.min_non_target_aligned_split_reads, type=int,
        help='The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local '
        'alignment to the target region to call a breakpoint by split read evidence')
    parse_group.add_argument(
        '--min_linking_split_reads', default=default.min_linking_split_reads, type=int,
        help='The minimum number of split reads which aligned to both breakpoints')
    parse_group.add_argument(
        '--min_flanking_reads_resolution', default=default.min_flanking_reads_resolution, type=int,
        help='the minimum number of flanking reads required to call a breakpoint by flanking evidence')


def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
        '-o', '--output',
        help='path to the output directory', required=True
    )
    parser.add_argument(
        '-n', '--input',
        help='path to the input file', required=True
    )
    parser.add_argument(
        '-b', '--bamfile',
        help='path to the input bam file', required=True
    )
    parser.add_argument(
        '--stranded', default=False,
        help='indicates that the input bam file is strand specific'
    )
    parser.add_argument(
        '-l', '--library',
        help='library id', required=True
    )
    g = parser.add_argument_group('reference files')
    g.add_argument(
        '-m', '--masking_file',
        default='/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv',
        help='path to the masking regions file'
    )
    g.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )
    g.add_argument(
        '-r', '--reference_genome',
        default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
        help='path to the human reference genome in fa format'
    )
    parser.add_argument('--igv_genome', help='the genome name to use for the igv batch file output', default='hg19')
    parser.add_argument(
        '-p', '--protocol',
        default=PROTOCOL.GENOME,
        choices=[PROTOCOL.GENOME, PROTOCOL.TRANS]
    )
    g = parser.add_argument_group('evidence settings')
    add_evidence_args_to_parser(g)
    args = parser.parse_args()
    return args


def gather_evidence_from_bam(clusters):
    evidence = []

    for i, e in enumerate(clusters):
        if e.protocol == PROTOCOL.GENOME:
            print()
            log(
                '({} of {})'.format(i + 1, len(clusters)),
                'gathering evidence for:',
                e.breakpoint_pair
            )
            log('possible event type(s):', BreakpointPair.classify(e.breakpoint_pair), time_stamp=False)
            try:
                e.load_evidence()
            except NotImplementedError as err:
                log(repr(err), time_stamp=False)
                continue
            log(
                'flanking reads:', [len(a) for a in e.flanking_reads],
                'split reads:', [len(a) for a in e.split_reads],
                time_stamp=False
            )
            e.assemble_split_reads()
            log('assembled {} contigs'.format(len(e.contigs)), time_stamp=False)
            evidence.append(e)
        else:
            raise NotImplementedError('currently only genome protocols are supported')
    return evidence


def main():
    global INPUT_BAM_CACHE, REFERENCE_ANNOTATIONS, MASKED_REGIONS, HUMAN_REFERENCE_GENOME, EVIDENCE_SETTINGS
    """
    - read the evidence
    - assemble contigs from the split reads
    - blat the contigs
    - pair the blatted contigs (where appropriate)
    - TODO: call the breakpoints and summarize the evidence
    """
    args = parse_arguments()
    EVIDENCE_SETTINGS = EvidenceSettings.parse_args(args)
    FILENAME_PREFIX = re.sub('\.(txt|tsv|tab)$', '', os.path.basename(args.input))
    EVIDENCE_BAM = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bam')
    CONTIG_BAM = os.path.join(args.output, FILENAME_PREFIX + '.contigs.bam')
    EVIDENCE_BED = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bed')

    PASSED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-passed.tab')
    PASSED_BED_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-passed.bed')
    FAILED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-failed.tab')
    CONTIG_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.contigs.tab')
    IGV_BATCH_FILE = os.path.join(args.output, FILENAME_PREFIX + '.igv.batch')
    INPUT_BAM_CACHE = BamCache(args.bamfile)

    log('input arguments listed below')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', val, time_stamp=False)
    log('loading the masking regions:', args.masking_file)
    MASKED_REGIONS = load_masking_regions(args.masking_file)
    for chr in MASKED_REGIONS:
        for m in MASKED_REGIONS[chr]:
            if m.name == 'nspan':
                m.position.start -= EVIDENCE_SETTINGS.read_length
                m.position.end += EVIDENCE_SETTINGS.read_length

    # load the reference genome
    log('loading the reference genome:', args.reference_genome)
    HUMAN_REFERENCE_GENOME = load_reference_genome(args.reference_genome)
    if args.protocol == PROTOCOL.TRANS:
        log('loading the reference annotations:', args.annotations)
        REFERENCE_ANNOTATIONS = load_reference_genes(args.annotations)
    log('loading complete')

    evidence_reads = set()

    split_read_contigs = set()
    chr_to_index = {}

    clusters = read_cluster_file(args.input, args.stranded)
    failed_cluster_rows = []
    filtered_clusters = []
    for cluster in clusters:
        overlaps_mask = None
        for mask in MASKED_REGIONS.get(cluster.break1.chr, []):
            if Interval.overlaps(cluster.window1, mask):
                overlaps_mask = mask
                break
        for mask in MASKED_REGIONS.get(cluster.break2.chr, []):
            if Interval.overlaps(cluster.window2, mask):
                overlaps_mask = mask
                break
        if overlaps_mask is None:
            filtered_clusters.append(cluster)
        else:
            log('dropping cluster {} overlapping mask {}:{}-{}'.format(
                cluster.breakpoint_pair, mask.reference_object, mask.start, mask.end))
            row = {}
            row.update(cluster.data)
            row.update(cluster.breakpoint_pair.flatten())
            fl = set([r.query_name for r in cluster.flanking_reads[0]]) | \
                set([r.query_name for r in cluster.flanking_reads[1]])
            row[COLUMNS.raw_flanking_reads] = len(fl)
            row[COLUMNS.raw_break1_split_reads] = len(cluster.split_reads[0])
            row[COLUMNS.raw_break2_split_reads] = len(cluster.split_reads[1])
            row['failure_comment'] = 'dropped b/c overlapped a masked region {}:{}-{}'.format(
                mask.reference_object, mask.start, mask.end
            )
            failed_cluster_rows.append(row)

    evidence = gather_evidence_from_bam(filtered_clusters)

    # output all the assemblies to a file of contigs
    with open(CONTIG_OUTPUT_FILE, 'w') as fh:
        log('writing the contigs to an output file:', CONTIG_OUTPUT_FILE)
        fh.write('#{}\t{}\t{}\n'.format(COLUMNS.cluster_id, COLUMNS.contig_sequence, COLUMNS.contig_remap_score))
        for ev in evidence:
            for c in ev.contigs:
                fh.write('{}\t{}\t{}\n'.format(ev.data[COLUMNS.cluster_id], c.seq, c.remap_score()))

    blat_sequences = set()
    for e in evidence:
        for c in e.contigs:
            blat_sequences.add(c.seq)
    print()
    log('aligning {} contig sequences'.format(len(blat_sequences)))
    blat_contig_alignments = blat_contigs(
        evidence,
        INPUT_BAM_CACHE,
        REFERENCE_GENOME=HUMAN_REFERENCE_GENOME
    )
    log('alignment complete')
    event_calls = []
    with open(EVIDENCE_BED, 'w') as fh:
        for e in evidence:
            fh.write('{}\t{}\t{}\t{}\n'.format(
                e.break1.chr, e.window1.start, e.window1.end, e.data[COLUMNS.cluster_id]))
            fh.write('{}\t{}\t{}\t{}\n'.format(
                e.break2.chr, e.window2.start, e.window2.end, e.data[COLUMNS.cluster_id]))
            print()
            log('calling events for', e.breakpoint_pair)
            calls = []
            failure_comment = None
            try:
                calls = e.call_events()
                event_calls.extend(calls)
            except UserWarning as err:
                log('warning: error in calling events', repr(err), time_stamp=False)
                failure_comment = str(err)

            if failure_comment:
                row = {}
                row.update(e.data)
                row.update(e.breakpoint_pair.flatten())
                fl = set([r.query_name for r in e.flanking_reads[0]]) | set([r.query_name for r in e.flanking_reads[1]])
                row[COLUMNS.raw_flanking_reads] = len(fl)
                row[COLUMNS.raw_break1_split_reads] = len(e.split_reads[0])
                row[COLUMNS.raw_break2_split_reads] = len(e.split_reads[1])
                row['failure_comment'] = failure_comment
                failed_cluster_rows.append(row)

            log('called {} event(s)'.format(len(calls)))

    # write the output validated clusters (split by type and contig)

    id_prefix = re.sub(' ', '_', str(datetime.now()))
    id = 1
    with open(PASSED_OUTPUT_FILE, 'w') as fh:
        print()
        log('writing:', PASSED_OUTPUT_FILE)
        rows = []
        header = set()
        for ec in event_calls:
            flank_count, flank_median, flank_stdev = ec.count_flanking_support()
            b1_count, b1_custom, b2_count, b2_custom, link_count = ec.count_split_read_support()
            b1_homseq = None
            b2_homseq = None
            try:
                b1_homseq, b2_homseq = ec.breakpoint_sequence_homology(HUMAN_REFERENCE_GENOME)
            except AttributeError:
                pass
            row = {
                COLUMNS.cluster_id: ec.data[COLUMNS.cluster_id],
                COLUMNS.validation_id: 'validation_{}-{}'.format(id_prefix, id),
                COLUMNS.break1_chromosome: ec.break1.chr,
                COLUMNS.break1_position_start: ec.break1.start,
                COLUMNS.break1_position_end: ec.break1.end,
                COLUMNS.break1_strand: STRAND.NS,
                COLUMNS.break1_orientation: ec.break1.orient,
                COLUMNS.break1_sequence: ec.break1.seq,
                COLUMNS.break2_chromosome: ec.break2.chr,
                COLUMNS.break2_position_start: ec.break2.start,
                COLUMNS.break2_position_end: ec.break2.end,
                COLUMNS.break2_strand: STRAND.NS,
                COLUMNS.break2_orientation: ec.break2.orient,
                COLUMNS.break2_sequence: ec.break2.seq,
                COLUMNS.event_type: ec.classification,
                COLUMNS.opposing_strands: ec.opposing_strands,
                COLUMNS.stranded: ec.stranded,
                COLUMNS.protocol: ec.evidence.protocol,
                COLUMNS.tools: ec.data[COLUMNS.tools],
                COLUMNS.contigs_assembled: len(ec.evidence.contigs),
                COLUMNS.contigs_aligned: sum([len(c.alignments) for c in ec.evidence.contigs]),
                COLUMNS.contig_sequence: None,
                COLUMNS.contig_remap_score: None,
                COLUMNS.contig_alignment_score: None,
                COLUMNS.break1_call_method: ec.call_method[0],
                COLUMNS.break2_call_method: ec.call_method[1],
                COLUMNS.flanking_reads: flank_count,
                COLUMNS.median_insert_size: round(flank_median, 0) if flank_median is not None else None,
                COLUMNS.stdev_insert_size: round(flank_stdev, 0) if flank_stdev is not None else None,
                COLUMNS.break1_split_reads: b1_count,
                COLUMNS.break1_split_reads_forced: b1_custom,
                COLUMNS.break2_split_reads: b2_count,
                COLUMNS.break2_split_reads_forced: b2_custom,
                COLUMNS.linking_split_reads: link_count,
                COLUMNS.untemplated_sequence: None,
                COLUMNS.break1_homologous_sequence: b1_homseq,
                COLUMNS.break2_homologous_sequence: b2_homseq,
                COLUMNS.break1_ewindow: '{}-{}'.format(*ec.evidence.window1),
                COLUMNS.break2_ewindow: '{}-{}'.format(*ec.evidence.window2),
                COLUMNS.break1_ewindow_count: ec.evidence.counts[0],
                COLUMNS.break2_ewindow_count: ec.evidence.counts[1]
            }
            if ec.contig:
                row[COLUMNS.contig_sequence] = ec.contig.seq
                row[COLUMNS.contig_remap_score] = ec.contig.remap_score()
                if ec.break1.strand == STRAND.NEG and not ec.stranded:
                    row[COLUMNS.contig_sequence] = reverse_complement(row[COLUMNS.contig_sequence])
            if ec.alignment:
                r1, r2 = ec.alignment
                if r2 is None:
                    row[COLUMNS.contig_alignment_score] = r1.get_tag('br')
                else:
                    row[COLUMNS.contig_alignment_score] = int(round((r1.get_tag('br') + r2.get_tag('br')) / 2, 0))
            if ec.untemplated_sequence is not None:
                row[COLUMNS.untemplated_sequence] = ec.untemplated_sequence
            if ec.stranded:
                row[COLUMNS.break1_strand] = ec.break1.strand
                row[COLUMNS.break2_strand] = ec.break2.strand
            rows.append(row)
            header.update(row.keys())
            id += 1
        header = sort_columns(header)
        fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row[col]) for col in header]) + '\n')

    with open(FAILED_OUTPUT_FILE, 'w') as fh:
        log('writing:', FAILED_OUTPUT_FILE)
        rows = []
        header = set()
        for row in failed_cluster_rows:
            header.update(row.keys())

        header = sort_columns(header)
        fh.write('#' + '\t'.join([str(c) for c in header]) + '\n')
        for row in failed_cluster_rows:
            fh.write('\t'.join([str(row.get(col, None)) for col in header]) + '\n')

    with pysam.AlignmentFile(CONTIG_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', CONTIG_BAM)
        for ev in evidence:
            for c in ev.contigs:
                for read1, read2 in c.alignments:
                    read1.cigar = CigarTools.convert_for_igv(read1.cigar)
                    fh.write(read1)
                    if read2:
                        read2.cigar = CigarTools.convert_for_igv(read2.cigar)
                        fh.write(read2)

    # write the evidence
    with pysam.AlignmentFile(EVIDENCE_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', EVIDENCE_BAM)
        reads = set()
        for ev in evidence:
            temp = ev.supporting_reads()
            reads.update(temp)
        for read in reads:
            read.cigar = CigarTools.convert_for_igv(read.cigar)
            fh.write(read)
    # now sort the contig bam
    sort = re.sub('.bam$', '.sorted', CONTIG_BAM)
    log('sorting the bam file:', CONTIG_BAM)
    subprocess.call(['samtools', 'sort', CONTIG_BAM, sort])
    CONTIG_BAM = sort + '.bam'
    log('indexing the sorted bam:', CONTIG_BAM)
    subprocess.call(['samtools', 'index', CONTIG_BAM])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted', EVIDENCE_BAM)
    log('sorting the bam file:', EVIDENCE_BAM)
    subprocess.call(['samtools', 'sort', EVIDENCE_BAM, sort])
    EVIDENCE_BAM = sort + '.bam'
    log('indexing the sorted bam:', EVIDENCE_BAM)
    subprocess.call(['samtools', 'index', EVIDENCE_BAM])

    # write the igv batch file
    with open(IGV_BATCH_FILE, 'w') as fh:
        log('writing:', IGV_BATCH_FILE)
        fh.write('new\ngenome {}\n'.format(args.igv_genome))

        fh.write('load {} name="{}"\n'.format(PASSED_BED_FILE, 'passed events'))
        fh.write('load {} name="{}"\n'.format(CONTIG_BAM, 'aligned contigs'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BED, 'evidence windows'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BAM, 'raw evidence'))
        fh.write('load {} name="{} {} input"\n'.format(args.bamfile, args.library, args.protocol))

    INPUT_BAM_CACHE.close()

if __name__ == '__main__':
    main()
