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
    |       |-- clusterset-#.validation.failed
    |       |-- clusterset-#.validation.passed
    |       |-- clusterset-#.contigs.bam
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
from tools import profile_bam
import math
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
EVIDENCE_SETTINGS = EvidenceSettings(median_insert_size=385, stdev_isize=95)


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
            COLUMNS.cluster_id.name,
            COLUMNS.cluster_size.name,
            COLUMNS.break1_chromosome.name,
            COLUMNS.break1_position_start.name,
            COLUMNS.break1_position_end.name,
            COLUMNS.break1_orientation.name,
            COLUMNS.break1_strand.name,
            COLUMNS.break2_chromosome.name,
            COLUMNS.break2_position_start.name,
            COLUMNS.break2_position_end.name,
            COLUMNS.break2_orientation.name,
            COLUMNS.break2_strand.name,
            COLUMNS.opposing_strands.name
        ],
        cast={
            COLUMNS.cluster_size.name: int,
            COLUMNS.break1_position_start.name: int,
            COLUMNS.break1_position_end.name: int,
            COLUMNS.break2_position_start.name: int,
            COLUMNS.break2_position_end.name: int,
            COLUMNS.opposing_strands.name: TSV.bool
        }
    )
    evidence = []
    for row in rows:
        strands = [(row[COLUMNS.break1_strand.name], row[COLUMNS.break2_strand.name])]
        if is_stranded:
            strands = itertools.product(
                STRAND.expand(row[COLUMNS.break1_strand.name]),
                STRAND.expand(row[COLUMNS.break2_strand.name])
            )
            strands = [(s1, s2) for s1, s2 in strands if row[COLUMNS.opposing_strands.name] == (s1 != s2)]
        if len(strands) == 0:
            raise UserWarning('error in reading input file. could not resolve strands', row)

        for s1, s2 in strands:
            bpp = BreakpointPair(
                Breakpoint(
                    row[COLUMNS.break1_chromosome.name],
                    row[COLUMNS.break1_position_start.name],
                    row[COLUMNS.break1_position_end.name],
                    strand=s1,
                    orient=row[COLUMNS.break1_orientation.name]
                ),
                Breakpoint(
                    row[COLUMNS.break2_chromosome.name],
                    row[COLUMNS.break2_position_start.name],
                    row[COLUMNS.break2_position_end.name],
                    strand=s2,
                    orient=row[COLUMNS.break2_orientation.name]
                ),
                opposing_strands=row[COLUMNS.opposing_strands.name]
            )
            try:
                e = Evidence(
                    bpp,
                    INPUT_BAM_CACHE,
                    HUMAN_REFERENCE_GENOME,
                    annotations=REFERENCE_ANNOTATIONS,
                    protocol=row[COLUMNS.protocol.name],
                    data=row
                )
                evidence.append(e)
            except UserWarning as e:
                warnings.warn('failed to read cluster {}'.format(repr(e)))
    return evidence


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
    parser.add_argument(
        '-m', '--masking_file',
        default='/home/creisle/svn/svmerge/trunk/hg19_masked_regions.tsv',
        help='path to the masking regions file'
    )
    parser.add_argument(
        '-a', '--annotations',
        default='/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv',
        help='path to the reference annotations of genes, transcript, exons, domains, etc.'
    )
    parser.add_argument(
        '-r', '--reference_genome',
        default='/home/pubseq/genomes/Homo_sapiens/TCGA_Special/GRCh37-lite.fa',
        help='path to the human reference genome in fa format'
    )
    parser.add_argument(
        '-p', '--protocol',
        default=PROTOCOL.GENOME,
        choices=[PROTOCOL.GENOME, PROTOCOL.TRANS]
    )
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
            ihist = {}
            for read in itertools.chain.from_iterable(e.flanking_reads):
                isize = abs(read.template_length)
                ihist[isize] = ihist.get(isize, 0) + 1
            try:
                median = profile_bam.histogram_median(ihist)
                stdev = math.sqrt(profile_bam.histogram_stderr(ihist, median))
                log('insert size: {:.0f} +/- {:.2f}'.format(median, stdev), time_stamp=False)
            except:
                pass
        else:
            raise NotImplementedError('currently only genome protocols are supported')
    return evidence


def main():
    global INPUT_BAM_CACHE, REFERENCE_ANNOTATIONS, MASKED_REGIONS, HUMAN_REFERENCE_GENOME
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
    PASSED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validated.passed')
    FAILED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validated.failed')
    MIN_EXTEND_OVERLAP = 6  # on each end
    MIN_CONTIG_READ_REMAP = 3
    MIN_BREAKPOINT_RESOLUTION = 3
    INPUT_BAM_CACHE = BamCache(args.bamfile)
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
            row[COLUMNS.raw_flanking_reads.name] = len(fl)
            row[COLUMNS.raw_break1_split_reads.name] = len(cluster.split_reads[0])
            row[COLUMNS.raw_break2_split_reads.name] = len(cluster.split_reads[1])
            row['failure_comment'] = 'dropped b/c overlapped a masked region {}:{}-{}'.format(
                mask.reference_object, mask.start, mask.end
            )
            failed_cluster_rows.append(row)

    evidence = gather_evidence_from_bam(filtered_clusters)

    blat_sequences = set()
    for e in evidence:
        for c in e.contigs:
            blat_sequences.add(c.seq)
    print()
    log('aligning {} contig sequences'.format(len(blat_sequences)))
    blat_contig_alignments = blat_contigs(
        evidence,
        INPUT_BAM_CACHE,
        reference_genome=HUMAN_REFERENCE_GENOME
    )
    log('alignment complete')
    event_calls = []
    for e in evidence:
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
            row[COLUMNS.raw_flanking_reads.name] = len(fl)
            row[COLUMNS.raw_break1_split_reads.name] = len(e.split_reads[0])
            row[COLUMNS.raw_break2_split_reads.name] = len(e.split_reads[1])
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
                COLUMNS.cluster_id.name: ec.data[COLUMNS.cluster_id.name],
                COLUMNS.validation_id.name: 'validation_{}-{}'.format(id_prefix, id),
                COLUMNS.break1_chromosome.name: ec.break1.chr,
                COLUMNS.break1_position_start.name: ec.break1.start,
                COLUMNS.break1_position_end.name: ec.break1.end,
                COLUMNS.break1_strand.name: STRAND.NS,
                COLUMNS.break1_orientation.name: ec.break1.orient,
                COLUMNS.break1_sequence.name: ec.break1.seq,
                COLUMNS.break2_chromosome.name: ec.break2.chr,
                COLUMNS.break2_position_start.name: ec.break2.start,
                COLUMNS.break2_position_end.name: ec.break2.end,
                COLUMNS.break2_strand.name: STRAND.NS,
                COLUMNS.break2_orientation.name: ec.break2.orient,
                COLUMNS.break2_sequence.name: ec.break2.seq,
                COLUMNS.event_type.name: ec.classification,
                COLUMNS.opposing_strands.name: ec.opposing_strands,
                COLUMNS.stranded.name: ec.stranded,
                COLUMNS.protocol.name: ec.evidence.protocol,
                COLUMNS.tools.name: ec.data[COLUMNS.tools.name],
                COLUMNS.contigs_assembled.name: len(ec.evidence.contigs),
                COLUMNS.contigs_aligned.name: sum([len(c.alignments) for c in ec.evidence.contigs]),
                COLUMNS.contig_sequence.name: None,
                COLUMNS.contig_remap_score.name: None,
                COLUMNS.contig_alignment_score.name: None,
                COLUMNS.break1_call_method.name: ec.call_method[0],
                COLUMNS.break2_call_method.name: ec.call_method[1],
                COLUMNS.flanking_reads.name: flank_count,
                COLUMNS.median_insert_size.name: round(flank_median, 0) if flank_median is not None else None,
                COLUMNS.stdev_insert_size.name: round(flank_stdev, 0) if flank_stdev is not None else None,
                COLUMNS.break1_split_reads.name: b1_count,
                COLUMNS.break1_split_reads_forced.name: b1_custom,
                COLUMNS.break2_split_reads.name: b2_count,
                COLUMNS.break2_split_reads_forced.name: b2_custom,
                COLUMNS.linking_split_reads.name: link_count,
                COLUMNS.untemplated_sequence.name: None,
                COLUMNS.break1_homologous_sequence.name: b1_homseq,
                COLUMNS.break2_homologous_sequence.name: b2_homseq,
                COLUMNS.break1_ewindow.name: '{}-{}'.format(*ec.evidence.window1),
                COLUMNS.break2_ewindow.name: '{}-{}'.format(*ec.evidence.window2),
                COLUMNS.break1_ewindow_count.name: ec.evidence.counts[0],
                COLUMNS.break2_ewindow_count.name: ec.evidence.counts[1]
            }
            if ec.contig:
                row[COLUMNS.contig_sequence.name] = ec.contig.seq
                row[COLUMNS.contig_remap_score.name] = ec.contig.remap_score()
                if ec.break1.strand == STRAND.NEG and not ec.stranded:
                    row[COLUMNS.contig_sequence.name] = reverse_complement(row[COLUMNS.contig_sequence.name])
            if ec.alignment:
                r1, r2 = ec.alignment
                if r2 is None:
                    row[COLUMNS.contig_alignment_score.name] = r1.get_tag('br')
                else:
                    row[COLUMNS.contig_alignment_score.name] = int(round((r1.get_tag('br') + r2.get_tag('br')) / 2, 0))
            if ec.untemplated_sequence is not None:
                row[COLUMNS.untemplated_sequence.name] = ec.untemplated_sequence
            if ec.stranded:
                row[COLUMNS.break1_strand.name] = ec.break1.strand
                row[COLUMNS.break2_strand.name] = ec.break2.strand
            rows.append(row)
            header.update(row.keys())
        header = sort_columns(header)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join([str(row[col]) for col in header]) + '\n')
            id += 1

    with open(FAILED_OUTPUT_FILE, 'w') as fh:
        log('writing:', FAILED_OUTPUT_FILE)
        rows = []
        header = set()
        for row in failed_cluster_rows:
            header.update(row.keys())

        header = sort_columns(header)
        fh.write('#' + '\t'.join(header) + '\n')
        for row in failed_cluster_rows:
            fh.write('\t'.join([str(row.get(col, None)) for col in header]) + '\n')
            id += 1

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
    log(' indexing the sorted bam:', CONTIG_BAM)
    subprocess.call(['samtools', 'index', CONTIG_BAM])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted', EVIDENCE_BAM)
    log('sorting the bam file:', EVIDENCE_BAM)
    subprocess.call(['samtools', 'sort', EVIDENCE_BAM, sort])
    EVIDENCE_BAM = sort + '.bam'
    log('indexing the sorted bam:', EVIDENCE_BAM)
    subprocess.call(['samtools', 'index', EVIDENCE_BAM])

    INPUT_BAM_CACHE.close()

if __name__ == '__main__':
    main()
