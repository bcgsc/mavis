from argparse import Namespace
import os
import pysam
import re
import subprocess
import sys
import itertools
import json
import networkx as nx
from Bio import SeqIO


# local modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from mavis.annotate import load_reference_genes, load_reference_genome, load_masking_regions, load_templates
from mavis.annotate.variant import annotate_events, determine_prime
from mavis.bam.cache import BamCache
from mavis.blat import blat_contigs
from mavis.breakpoint import BreakpointPair
from mavis.cluster import cluster_breakpoint_pairs
from mavis.constants import PROTOCOL, COLUMNS, PRIME, CALL_METHOD, SVTYPE, SPLICE_TYPE
from mavis.error import DrawingFitError, NotSpecifiedError
from mavis.illustrate.constants import DiagramSettings
from mavis.illustrate.diagram import draw_sv_summary_diagram
from mavis.interval import Interval
from mavis.validate.call import call_events
from mavis.validate.constants import VALIDATION_DEFAULTS
from mavis.validate.evidence import GenomeEvidence, TranscriptomeEvidence
import mavis.bam.cigar as cigar_tools
import mavis.pipeline.config as pconf
from mavis.pipeline.util import *
from mavis.pairing import equivalent_events

VALIDATION_PASS_SUFFIX = '.validation-passed.tab'

QSUB_HEADER = """#!/bin/bash
#$ -V
#$ -N {name}
#$ -q {queue}
#$ -o {output}
#$ -l mem_free={memory}G,mem_token={memory}G,h_vmem={memory}G
#$ -j y"""


def main_validate(args):
    """
    - read the evidence
    - assemble contigs from the split reads
    - blat the contigs
    - pair the blatted contigs (where appropriate)
    - TODO: call the breakpoints and summarize the evidence
    """
    FILENAME_PREFIX = re.sub('\.(txt|tsv|tab)$', '', os.path.basename(args.input))
    RAW_EVIDENCE_BAM = os.path.join(args.output, FILENAME_PREFIX + '.raw_evidence.bam')
    CONTIG_BAM = os.path.join(args.output, FILENAME_PREFIX + '.contigs.bam')
    EVIDENCE_BED = os.path.join(args.output, FILENAME_PREFIX + '.evidence.bed')

    PASSED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + VALIDATION_PASS_SUFFIX)
    PASSED_BED_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-passed.bed')
    FAILED_OUTPUT_FILE = os.path.join(args.output, FILENAME_PREFIX + '.validation-failed.tab')
    CONTIG_BLAT_FA = os.path.join(args.output, FILENAME_PREFIX + '.contigs.fa')
    CONTIG_BLAT_OUTPUT = os.path.join(args.output, FILENAME_PREFIX + '.contigs.blat_out.pslx')
    IGV_BATCH_FILE = os.path.join(args.output, FILENAME_PREFIX + '.igv.batch')
    INPUT_BAM_CACHE = BamCache(args.bam_file, args.stranded_bam)

    evidence_reads = set()  # keep track of collected reads to use for ouput

    split_read_contigs = set()
    chr_to_index = {}
    bpps = read_inputs(
        [args.input], add={COLUMNS.protocol: args.protocol, COLUMNS.library: args.library})
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            e = GenomeEvidence(
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                args.reference_genome[1],
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=args.stdev_fragment_size,
                read_length=args.read_length,
                median_fragment_size=args.median_fragment_size
            )
            evidence_clusters.append(e)
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            e = TranscriptomeEvidence(
                args.annotations[1],
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                args.reference_genome[1],
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=args.stdev_fragment_size,
                read_length=args.read_length,
                median_fragment_size=args.median_fragment_size
            )
            evidence_clusters.append(e)
        else:
            raise ValueError('protocol not recognized', bpp.data[COLUMNS.protocol])

    for chr, masks in args.masking[1].items():  # extend masking by read length
        for mask in masks:
            mask.position.start -= args.read_length
            mask.position.end += args.read_length

    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(evidence_clusters, args.masking[1])

    for i, e in enumerate(evidence_clusters):
        print()
        log(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathering evidence for:',
            e.data['cluster_id']
        )
        log(e, time_stamp=False)
        log('possible event type(s):', BreakpointPair.classify(e), time_stamp=False)
        log('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.outer_window1[0], e.outer_window1[1],
            e.break2.chr, e.outer_window2[0], e.outer_window2[1]), time_stamp=False)
        log('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.inner_window1[0], e.inner_window1[1],
            e.break2.chr, e.inner_window2[0], e.inner_window2[1]), time_stamp=False)
        e.load_evidence(log=log)
        log(
            'flanking pairs: {};'.format(len(e.flanking_pairs)),
            'split reads: {}, {};'.format(*[len(a) for a in e.split_reads]),
            'half-mapped reads: {}, {};'.format(*[len(a) for a in e.half_mapped]),
            'spanning-reads: {};'.format(len(e.spanning_reads)),
            'compatible flanking pairs:', len(e.compatible_flanking_pairs),
            time_stamp=False
        )
        e.assemble_contig(log=log)
        log('assembled {} contigs'.format(len(e.contigs)), time_stamp=False)

    log('will output:', CONTIG_BLAT_FA, CONTIG_BLAT_OUTPUT)
    blat_contigs(
        evidence_clusters,
        INPUT_BAM_CACHE,
        reference_genome=args.reference_genome[1],
        blat_2bit_reference=args.blat_2bit_reference,
        blat_fa_input_file=CONTIG_BLAT_FA,
        blat_pslx_output_file=CONTIG_BLAT_OUTPUT,
        clean_files=False,
        blat_min_percent_of_max_score=args.blat_min_percent_of_max_score,
        blat_min_identity=args.blat_min_identity,
        blat_min_query_consumption=args.blat_min_query_consumption
    )
    log('alignment complete')
    event_calls = []
    passes = 0
    write_bed_file(EVIDENCE_BED, itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]))
    for index, e in enumerate(evidence_clusters):
        print()
        log('({} of {}) calling events for:'.format
            (index + 1, len(evidence_clusters)), e.data[COLUMNS.cluster_id], e.putative_event_types())
        log('source:', e, time_stamp=False)
        calls = []
        failure_comment = None
        try:
            calls.extend(call_events(e))
            event_calls.extend(calls)
        except UserWarning as err:
            log('warning: error in calling events', repr(err), time_stamp=False)
            failure_comment = str(err)
        if len(calls) == 0:
            failure_comment = ['zero events were called'] if failure_comment is None else failure_comment
            e.data[COLUMNS.filter_comment] = failure_comment
            filtered_evidence_clusters.append(e)
        else:
            passes += 1

        log('called {} event(s)'.format(len(calls)))
        for ev in calls:
            log(ev, time_stamp=False)
            log(ev.event_type, ev.call_method, time_stamp=False)
            log('remapped reads: {}; split reads: [{}, {}], flanking pairs: {}'.format(
                0 if not ev.contig else len(ev.contig.input_reads),
                len(ev.break1_split_reads), len(ev.break2_split_reads),
                len(ev.flanking_pairs)), time_stamp=False)

    if len(filtered_evidence_clusters) + passes != len(evidence_clusters):
        raise AssertionError(
            'totals do not match pass + fails == total', passes, len(filtered_evidence_clusters), len(evidence_clusters))
    # write the output validated clusters (split by type and contig)
    validation_batch_id = build_batch_id(prefix='validation-')
    for i, ec in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = ec.breakpoint_sequence_homology(args.reference_genome[1])
        except AttributeError:
            pass
        ec.data.update({
            COLUMNS.validation_id: '{}-{}'.format(validation_batch_id, i + 1),
            COLUMNS.break1_homologous_seq: b1_homseq,
            COLUMNS.break2_homologous_seq: b2_homseq,
        })

    output_tabbed_file(event_calls, PASSED_OUTPUT_FILE)
    output_tabbed_file(filtered_evidence_clusters, FAILED_OUTPUT_FILE)

    with pysam.AlignmentFile(CONTIG_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', CONTIG_BAM)
        for ev in evidence_clusters:
            for c in ev.contigs:
                for read1, read2 in c.alignments:
                    read1.cigar = cigar_tools.convert_for_igv(read1.cigar)
                    fh.write(read1)
                    if read2:
                        read2.cigar = cigar_tools.convert_for_igv(read2.cigar)
                        fh.write(read2)

    # write the evidence
    with pysam.AlignmentFile(RAW_EVIDENCE_BAM, 'wb', template=INPUT_BAM_CACHE.fh) as fh:
        log('writing:', RAW_EVIDENCE_BAM)
        reads = set()
        for ev in evidence_clusters:
            temp = ev.supporting_reads()
            reads.update(temp)
        for read in reads:
            read.cigar = cigar_tools.convert_for_igv(read.cigar)
            fh.write(read)
    # now sort the contig bam
    sort = re.sub('.bam$', '.sorted.bam', CONTIG_BAM)
    log('sorting the bam file:', CONTIG_BAM)
    if args.samtools_version[0] < 1:
        subprocess.call(samtools_v0_sort(CONTIG_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(CONTIG_BAM, sort), shell=True)
    CONTIG_BAM = sort
    log('indexing the sorted bam:', CONTIG_BAM)
    subprocess.call(['samtools', 'index', CONTIG_BAM])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted.bam', RAW_EVIDENCE_BAM)
    log('sorting the bam file:', RAW_EVIDENCE_BAM)
    if args.samtools_version[0] < 1:
        subprocess.call(samtools_v0_sort(RAW_EVIDENCE_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(RAW_EVIDENCE_BAM, sort), shell=True)
    RAW_EVIDENCE_BAM = sort
    log('indexing the sorted bam:', RAW_EVIDENCE_BAM)
    subprocess.call(['samtools', 'index', RAW_EVIDENCE_BAM])

    # write the igv batch file
    with open(IGV_BATCH_FILE, 'w') as fh:
        log('writing:', IGV_BATCH_FILE)
        fh.write('new\ngenome {}\n'.format(args.reference_genome))

        fh.write('load {} name="{}"\n'.format(PASSED_BED_FILE, 'passed events'))
        fh.write('load {} name="{}"\n'.format(CONTIG_BAM, 'aligned contigs'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BED, 'evidence windows'))
        fh.write('load {} name="{}"\n'.format(RAW_EVIDENCE_BAM, 'raw evidence'))
        fh.write('load {} name="{} {} input"\n'.format(args.bam_file, args.library, args.protocol))


def main_pipeline(args, configs):
    # read the config
    # set up the directory structure and run svmerge
    annotation_files = []
    annotation_jobs = []
    for sec in configs:
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        log('setting up the directory structure for', sec.library, 'as', base)
        base = os.path.join(args.output, '{}_{}'.format(sec.library, sec.protocol))
        cluster_output = mkdirp(os.path.join(base, 'clustering'))
        validation_output = mkdirp(os.path.join(base, 'validation'))
        annotation_output = mkdirp(os.path.join(base, 'annotation'))

        # run the merge
        log('clustering')
        merge_args = {}
        merge_args.update(sec.__dict__)
        merge_args.update(args.__dict__)
        merge_args['output'] = os.path.join(base, 'clustering')
        output_files = main_cluster(Namespace(**merge_args))
        merge_file_prefix = None
        for f in output_files:
            m = re.match('^(?P<prefix>.*\D)\d+.tab$', f)
            if not m:
                raise UserWarning('cluster file did not match expected format', f)
            if merge_file_prefix is None:
                merge_file_prefix = m.group('prefix')
            elif merge_file_prefix != m.group('prefix'):
                raise UserWarning('merge file prefixes are not consistent', output_files)

        # now set up the qsub script for the validation and the held job for the annotation
        validation_args = {
            'output': validation_output,
            'masking': args.masking[0],
            'reference_genome': args.reference_genome[0],
            'blat_2bit_reference': args.blat_2bit_reference,
            'annotations': args.annotations[0],
            'library': sec.library,
            'bam_file': sec.bam_file,
            'protocol': sec.protocol,
            'read_length': sec.read_length,
            'stdev_fragment_size': sec.stdev_fragment_size,
            'median_fragment_size': sec.median_fragment_size,
            'stranded_bam': sec.stranded_bam,
            'protocol': sec.protocol
        }
        for attr in sorted(VALIDATION_DEFAULTS.__dict__.keys()):
            validation_args[attr] = getattr(sec, attr)

        qsub = os.path.join(validation_output, 'qsub.sh')
        validation_jobname = 'validation_{}_{}'.format(sec.library, sec.protocol)
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=args.queue, memory=args.validate_memory, name=validation_jobname, output=validation_output
                ) + '\n')
            fh.write('#$ -t {}-{}\n'.format(1, len(output_files)))
            temp = ['--{} {}'.format(k, v) for k, v in validation_args.items() if not isinstance(v, str) and v is not None]
            temp.extend(['--{} "{}"'.format(k, v) for k, v in validation_args.items() if isinstance(v, str) and v is not None])
            validation_args = temp
            validation_args.append('-n {}$SGE_TASK_ID.tab'.format(merge_file_prefix))
            fh.write('python {} validate {}\n'.format(__file__, ' \\\n\t'.join(validation_args)))

        # set up the annotations job
        # for all files with the right suffix
        annotation_args = {
            'output': annotation_output,
            'reference_genome': args.reference_genome[0],
            'annotations': args.annotations[0],
            'template_metadata': args.template_metadata[0],
            'min_orf_size': sec.min_orf_size,
            'max_orf_cap': sec.max_orf_cap,
            'min_domain_mapping_match': sec.min_domain_mapping_match,
            'domain_regex_filter': sec.domain_regex_filter,
            'max_proximity': sec.max_proximity
        }
        temp = ['--{} {}'.format(k, v) for k, v in annotation_args.items() if not isinstance(v, str) and v is not None]
        temp.extend(['--{} "{}"'.format(k, v) for k, v in annotation_args.items() if isinstance(v, str) and v is not None])
        annotation_args = temp
        annotation_args.append('--inputs {}/*{}*{}'.format(
            validation_output, os.path.basename(merge_file_prefix), VALIDATION_PASS_SUFFIX))
        qsub = os.path.join(annotation_output, 'qsub.sh')
        annotation_jobname = 'annotation_{}_{}'.format(sec.library, sec.protocol)
        annotation_jobs.append(annotation_jobname)
        annotation_files.append(os.path.join(annotation_output, 'annotations.tab'))
        with open(qsub, 'w') as fh:
            log('writing:', qsub)
            fh.write(
                QSUB_HEADER.format(
                    queue=args.queue, memory=args.default_memory, name=annotation_jobname, output=annotation_output
                ) + '\n')
            fh.write('#$ -hold_jid {}\n'.format(validation_jobname))
            fh.write('python {} annotate {}\n'.format(__file__, ' \\\n\t'.join(annotation_args)))

    # set up scripts for the pairing held on all of the annotation jobs
    pairing_output = mkdirp(os.path.join(args.output, 'pairing'))
    pairing_args = dict(
        output=pairing_output,
        split_call_distance=args.split_call_distance,
        contig_call_distance=args.contig_call_distance,
        flanking_call_distance=args.flanking_call_distance,
        max_proximity=args.max_proximity,
        annotations=args.annotations[0],
        low_memory=args.low_memory
    )
    temp = ['--{} {}'.format(k, v) for k, v in pairing_args.items() if not isinstance(v, str) and v is not None]
    temp.extend(['--{} "{}"'.format(k, v) for k, v in pairing_args.items() if isinstance(v, str) and v is not None])
    temp.append('--inputs {}'.format(' '.join(annotation_files)))
    pairing_args = temp
    qsub = os.path.join(pairing_output, 'qsub.sh')
    with open(qsub, 'w') as fh:
        log('writing:', qsub)
        fh.write(
            QSUB_HEADER.format(
                queue=args.queue, memory=args.default_memory, name='mavis_pairing', output=pairing_output
            ) + '\n')
        fh.write('#$ -hold_jid {}\n'.format(' '.join(annotation_jobs)))
        fh.write('python {} pairing {}\n'.format(__file__, ' \\\n\t'.join(pairing_args)))


def main_cluster(args):
    # output files
    cluster_batch_id = build_batch_id(prefix='cluster-')
    UNINFORM_OUTPUT = os.path.join(args.output, 'uninformative_clusters.txt')
    CLUSTER_ASSIGN_OUTPUT = os.path.join(args.output, 'cluster_assignment.tab')
    CLUSTER_BED_OUTPUT = os.path.join(args.output, 'clusters.bed')
    split_file_name_func = lambda x: os.path.join(args.output, '{}-{}.tab'.format(cluster_batch_id, x))
    # load the input files
    temp = read_inputs(
        args.inputs, args.stranded_bam,
        cast={COLUMNS.tools: lambda x: set(x.split(';')) if x else set()},
        add={COLUMNS.library: args.library, COLUMNS.protocol: args.protocol}
    )
    breakpoint_pairs = []
    for bpp in temp:
        if bpp.data[COLUMNS.library] == args.library and bpp.data[COLUMNS.protocol] == args.protocol:
            breakpoint_pairs.append(bpp)
    # filter by masking file
    breakpoint_pairs, filtered_bpp = filter_on_overlap(breakpoint_pairs, args.masking[1])

    log('computing clusters')
    clusters = cluster_breakpoint_pairs(breakpoint_pairs, r=args.cluster_radius, k=args.cluster_clique_size)

    hist = {}
    length_hist = {}
    for index, cluster in enumerate(clusters):
        input_pairs = clusters[cluster]
        hist[len(input_pairs)] = hist.get(len(input_pairs), 0) + 1
        c1 = round(len(cluster[0]), -2)
        c2 = round(len(cluster[1]), -2)
        length_hist[c1] = length_hist.get(c1, 0) + 1
        length_hist[c2] = length_hist.get(c2, 0) + 1
        cluster.data[COLUMNS.cluster_id] = '{}-{}'.format(cluster_batch_id, index + 1)
        cluster.data[COLUMNS.cluster_size] = len(input_pairs)
        temp = set()
        for p in input_pairs:
            temp.update(p.data[COLUMNS.tools])
        cluster.data[COLUMNS.tools] = ';'.join(sorted(list(temp)))
    log('computed', len(clusters), 'clusters', time_stamp=False)
    log('cluster input pairs distribution', sorted(hist.items()), time_stamp=False)
    log('cluster intervals lengths', sorted(length_hist.items()), time_stamp=False)
    # map input pairs to cluster ids
    # now create the mapping from the original input files to the cluster(s)
    mkdirp(args.output)

    with open(CLUSTER_ASSIGN_OUTPUT, 'w') as fh:
        header = set()
        rows = {}

        for cluster, input_pairs in clusters.items():
            for p in input_pairs:
                if p in rows:
                    rows[p][COLUMNS.tools].update(p.data[COLUMNS.tools])
                else:
                    rows[p] = BreakpointPair.flatten(p)
                rows[p].setdefault('clusters', set()).add(cluster.data[COLUMNS.cluster_id])
        for row in rows.values():
            row['clusters'] = ';'.join([str(c) for c in sorted(list(row['clusters']))])
            row[COLUMNS.tools] = ';'.join(sorted(list(row[COLUMNS.tools])))
            row[COLUMNS.library] = args.library
            row[COLUMNS.protocol] = args.protocol
        output_tabbed_file(rows, CLUSTER_ASSIGN_OUTPUT)

    output_files = []
    # filter clusters based on annotations
    # decide on the number of clusters to validate per job
    pass_clusters = list(clusters)
    fail_clusters = []

    if args.uninformative_filter:
        pass_clusters = []
        for cluster in clusters:
            # loop over the annotations
            overlaps_gene = False
            w1 = Interval(cluster.break1.start - args.max_proximity, cluster.break1.end + args.max_proximity)
            w2 = Interval(cluster.break2.start - args.max_proximity, cluster.break2.end + args.max_proximity)
            for gene in args.annotations[1].get(cluster.break1.chr, []):
                if Interval.overlaps(gene, w1):
                    overlaps_gene = True
                    break
            for gene in args.annotations[1].get(cluster.break2.chr, []):
                if Interval.overlaps(gene, w2):
                    overlaps_gene = True
                    break
            if overlaps_gene:
                pass_clusters.append(cluster)
            else:
                fail_clusters.append(cluster)
    if len(fail_clusters) + len(pass_clusters) != len(clusters):
        raise AssertionError('totals do not add up', len(fail_clusters), len(pass_clusters), 'does not total to', len(clusters))
    log('filtered', len(fail_clusters), 'clusters as not informative')
    output_tabbed_file(fail_clusters, UNINFORM_OUTPUT)

    JOB_SIZE = args.min_clusters_per_file
    if len(pass_clusters) // args.min_clusters_per_file > args.max_files - 1:
        JOB_SIZE = len(pass_clusters) // args.max_files
        assert(len(pass_clusters) // JOB_SIZE == args.max_files)

    bedfile = os.path.join(args.output, 'clusters.bed')
    write_bed_file(bedfile, itertools.chain.from_iterable([b.get_bed_repesentation() for b in pass_clusters]))

    job_ranges = list(range(0, len(pass_clusters), JOB_SIZE))
    if job_ranges[-1] != len(pass_clusters):
        job_ranges.append(len(pass_clusters))
    job_ranges = zip(job_ranges, job_ranges[1::])

    for i, jrange in enumerate(job_ranges):
        # generate an output file
        filename = split_file_name_func(i + 1)
        output_files.append(filename)
        output_tabbed_file(pass_clusters[jrange[0]:jrange[1]], filename)

    return output_files

def main_annotate(args):
    DRAWINGS_DIRECTORY = os.path.join(args.output, 'drawings')
    TABBED_OUTPUT_FILE = os.path.join(args.output, 'annotations.tab')
    FA_OUTPUT_FILE = os.path.join(args.output, 'annotations.fusion-cdna.fa')

    mkdirp(DRAWINGS_DIRECTORY)
    # test that the sequence makes sense for a random transcript
    bpps = read_inputs(args.inputs, in_={COLUMNS.protocol: PROTOCOL})
    log('read {} breakpoint pairs'.format(len(bpps)))

    annotations = annotate_events(
        bpps,
        reference_genome=args.reference_genome[1], annotations=args.annotations[1],
        min_orf_size=args.min_orf_size, min_domain_mapping_match=args.min_domain_mapping_match,
        max_orf_cap=args.max_orf_cap,
        log=log
    )

    id_prefix = build_batch_id(prefix='annotation-', suffix='-')
    rows = []  # hold the row information for the final tsv file
    fa_sequences = {}
    # now try generating the svg
    DS = DiagramSettings()
    DS.DOMAIN_NAME_REGEX_FILTER = args.domain_regex_filter

    for i, ann in enumerate(annotations):
        ann.data[COLUMNS.annotation_id] = id_prefix + str(i + 1)
        row = ann.flatten()
        row[COLUMNS.break1_strand] = ann.transcript1.get_strand()
        row[COLUMNS.break2_strand] = ann.transcript2.get_strand()
        row[COLUMNS.fusion_sequence_fasta_file] = FA_OUTPUT_FILE

        log('current annotation', ann.data[COLUMNS.annotation_id], ann.transcript1, ann.transcript2, ann.event_type)

        # try building the fusion product
        ann_rows = []
        # add fusion information to the current row
        transcripts = [] if not ann.fusion else ann.fusion.transcripts
        for t in transcripts:
            fusion_fa_id = '{}_{}'.format(ann.data[COLUMNS.annotation_id], t.splicing_pattern.splice_type)
            if fusion_fa_id in fa_sequences:
                raise AssertionError('should not be duplicate fa sequence ids', fusion_fa_id)
            fa_sequences[fusion_fa_id] = ann.fusion.get_cdna_seq(t.splicing_pattern)

            # duplicate the row for each translation
            for tl in t.translations:
                nrow = dict()
                nrow.update(row)
                nrow[COLUMNS.fusion_splicing_pattern] = tl.transcript.splicing_pattern.splice_type
                nrow[COLUMNS.fusion_cdna_coding_start] = tl.start
                nrow[COLUMNS.fusion_cdna_coding_end] = tl.end
                nrow[COLUMNS.fusion_sequence_fasta_id] = fusion_fa_id

                domains = []
                for dom in tl.domains:
                    m, t = dom.score_region_mapping()
                    temp = {
                        "name": dom.name,
                        "sequences": dom.get_seqs(),
                        "regions": [
                            {"start": dr.start, "end": dr.end} for dr in sorted(dom.regions, key=lambda x: x.start)
                        ],
                        "mapping_quality": round(m * 100 / t, 0),
                        "matches": m
                    }
                    domains.append(temp)
                nrow[COLUMNS.fusion_mapped_domains] = json.dumps(domains)
                ann_rows.append(nrow)

        drawing = None
        retry_count = 0
        while drawing is None:  # continue if drawing error and increase width
            try:
                canvas, legend = draw_sv_summary_diagram(
                    DS, ann, reference_genome=args.reference_genome[1], templates=args.template_metadata[1])

                gene_aliases1 = 'NA'
                gene_aliases2 = 'NA'
                try:
                    if len(ann.transcript1.gene.aliases) > 0:
                        gene_aliases1 = '-'.join(ann.transcript1.gene.aliases)
                    if ann.transcript1.is_best_transcript:
                        gene_aliases1 = 'b-' + gene_aliases1
                except AttributeError:
                    pass
                try:
                    if len(ann.transcript2.gene.aliases) > 0:
                        gene_aliases2 = '-'.join(ann.transcript2.gene.aliases)
                    if ann.transcript2.is_best_transcript:
                        gene_aliases2 = 'b-' + gene_aliases2
                except AttributeError:
                    pass
                try:
                    if determine_prime(ann.transcript1, ann.break1) == PRIME.THREE:
                        gene_aliases1, gene_aliases2 = gene_aliases2, gene_aliases1
                except NotSpecifiedError:
                    pass

                name = 'mavis_{}_{}.{}_{}.{}'.format(
                    ann.break1.chr, ann.break2.chr, gene_aliases1, gene_aliases2, ann.data[COLUMNS.annotation_id])

                drawing = os.path.join(DRAWINGS_DIRECTORY, name + '.svg')
                l = os.path.join(DRAWINGS_DIRECTORY, name + '.legend.json')
                for r in ann_rows + [row]:
                    r[COLUMNS.annotation_figure] = drawing
                    r[COLUMNS.annotation_figure_legend] = l
                log('generating svg:', drawing, time_stamp=False)
                canvas.saveas(drawing)

                log('generating legend:', l, time_stamp=False)
                with open(l, 'w') as fh:
                    json.dump(legend, fh)
                break
            except DrawingFitError as err:
                log('extending width:', DS.WIDTH, DS.WIDTH + 500, time_stamp=False)
                DS.WIDTH += 500
                retry_count += 1
                if retry_count > 10:
                    raise err
        if len(ann_rows) == 0:
            rows.append(row)
        else:
            rows.extend(ann_rows)

    output_tabbed_file(rows, TABBED_OUTPUT_FILE)

    with open(FA_OUTPUT_FILE, 'w') as fh:
        log('writing:', FA_OUTPUT_FILE)
        for name, seq in sorted(fa_sequences.items()):
            fh.write('> {}\n'.format(name))
            fh.write('{}\n'.format(seq))


def main_pairing(args):
    # load the file
    DISTANCES = {
        CALL_METHOD.FLANK: args.flanking_call_distance,
        CALL_METHOD.SPLIT: args.split_call_distance,
        CALL_METHOD.CONTIG: args.contig_call_distance
    }

    bpps = []
    bpps.extend(read_inputs(
        args.inputs,
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
        in_={
            COLUMNS.protocol: PROTOCOL,
            COLUMNS.event_type: SVTYPE,
            COLUMNS.break1_call_method: CALL_METHOD,
            COLUMNS.break2_call_method: CALL_METHOD,
            COLUMNS.fusion_splicing_pattern: SPLICE_TYPE.values() + [None, 'None']
        },
        add={
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.fusion_sequence_fasta_file: None,
            COLUMNS.fusion_splicing_pattern: None
        }
    ))
    log('read {} breakpoint pairs'.format(len(bpps)))
    libraries = set()

    SEQUENCES = dict()
    sequence_files = set()
    for bpp in bpps:
        if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
            sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
        libraries.add(bpp.data[COLUMNS.library])
    log('pairing between', len(libraries), 'libraries')
    for f in sorted(list(sequence_files)):
        log('loading:', f)
        with open(f, 'rU') as fh:
            temp = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
            for k in temp:
                if k in SEQUENCES:
                    raise AssertionError('sequence identifiers are not unique', k)
            SEQUENCES.update(temp)

    TRANSCRIPTS = dict()

    for chr, genes in args.annotations[1].items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in TRANSCRIPTS:
                    raise AssertionError('transcript is not unique', gene, t)
                TRANSCRIPTS[t.name] = t

    calls_by_lib = dict()
    pairkey_bpp_mapping = dict()
    pairing = dict()

    for bpp in bpps:
        lib_key = (bpp.data[COLUMNS.library], bpp.data[COLUMNS.protocol])
        pair_key = [
            bpp.data[COLUMNS.library],
            bpp.data[COLUMNS.protocol],
            bpp.data[COLUMNS.annotation_id],
            bpp.data[COLUMNS.fusion_splicing_pattern],
            bpp.data[COLUMNS.fusion_cdna_coding_start],
            bpp.data[COLUMNS.fusion_cdna_coding_end]
        ]
        pair_key = '_'.join([str(k) for k in pair_key if k is not None])
        bpp.data[COLUMNS.product_id] = pair_key
        calls_by_lib.setdefault(lib_key, set())
        calls_by_lib[lib_key].add(pair_key)

        if pair_key in calls_by_lib:
            raise KeyError('duplicate bpp is not unique within lib', pair_key, bpp, bpp.data)
        pairkey_bpp_mapping[pair_key] = bpp
        pairing[pair_key] = set()

    # pairwise comparison of breakpoints between all libraries
    for lib1, lib2 in itertools.combinations(calls_by_lib.keys(), 2):
        # for each two libraries pair all calls
        log(len(calls_by_lib[lib1]) * len(calls_by_lib[lib2]), 'comparison(s) between', lib1, 'and', lib2)
        for pkey1, pkey2 in itertools.product(calls_by_lib[lib1], calls_by_lib[lib2]):
            if equivalent_events(
                pairkey_bpp_mapping[pkey1],
                pairkey_bpp_mapping[pkey2],
                DISTANCES=DISTANCES,
                TRANSCRIPTS=TRANSCRIPTS,
                SEQUENCES=SEQUENCES
            ):
                pairing[pkey1].add(pkey2)
                pairing[pkey2].add(pkey1)
    
    for pkey, pairs in pairing.items():
        bpp = pairkey_bpp_mapping[pkey]
        # filter any matches where genes match but transcripts do not
        filtered = []
        for pkey2 in pairs:
            pair_bpp = pairkey_bpp_mapping[pkey2]
            if bpp.data[COLUMNS.gene1] and bpp.data[COLUMNS.gene1] == pair_bpp.data[COLUMNS.gene1]:
                if bpp.data[COLUMNS.transcript1] != pair_bpp.data[COLUMNS.transcript1]:
                    continue
            if bpp.data[COLUMNS.gene2] and bpp.data[COLUMNS.gene2] == pair_bpp.data[COLUMNS.gene2]:
                if bpp.data[COLUMNS.transcript2] != pair_bpp.data[COLUMNS.transcript2]:
                    continue
            filtered.append(pkey2)
        bpp.data[COLUMNS.pairing] = ';'.join(sorted(filtered))

    fname = os.path.join(
        args.output,
        'mavis_paired_{}.tab'.format('_'.join(sorted([l for l, p in calls_by_lib])))
    )
    output_tabbed_file(bpps, fname)


def main_summary(args):
    pass


def main():
    def usage(err):
        name = os.path.basename(__file__)
        print('usage: {} {{cluster,validate,annotate,pairing,summary,pipeline}}'.format(name))
        print('{}: error:'.format(name), err)
        exit(1)

    if len(sys.argv) < 2:
        usage('the <pipeline step> argument is required')

    pstep = sys.argv.pop(1)
    sys.argv[0] = '{} {}'.format(sys.argv[0], pstep)

    if pstep == PIPELINE_STEP.SUMMARY:
        raise NotImplementedError('summary script has not been written')
    args = pconf.parse_arguments(pstep)
    config = []
    log('input arguments')
    for arg, val in sorted(args.__dict__.items()):
        log(arg, '=', repr(val), time_stamp=False)
    if pstep == PIPELINE_STEP.PIPELINE:
        if args.write:
            log('writing:', args.config)
            pconf.write_config(args.config, include_defaults=True)
            exit()
        else:
            temp, config = pconf.read_config(args.config)
            args.__dict__.update(temp.__dict__)
            for sec in config:
                sec.output = args.output
    # load the reference files if they have been given and reset the arguments to hold the original file name and the
    # loaded data
    if any([
        pstep not in [PIPELINE_STEP.PIPELINE, PIPELINE_STEP.VALIDATE],
        hasattr(args, 'uninformative_filter') and args.uninformative_filter,
        pstep in [PIPELINE_STEP.VALIDATE, PIPELINE_STEP.CLUSTER] and args.protocol == PROTOCOL.TRANS,
        pstep == PIPELINE_STEP.PIPELINE and any([sec.protocol == PROTOCOL.TRANS for sec in config])
    ]):
        log('loading:', args.annotations)
        args.__dict__['annotations'] = args.annotations, load_reference_genes(args.annotations)
    else:
        args.__dict__['annotations'] = args.annotations, None

    try:
        if pstep == PIPELINE_STEP.PIPELINE:
            raise AttributeError()
        log('loading:' if not args.low_memory else 'indexing:', args.reference_genome)
        args.__dict__['reference_genome'] = args.reference_genome, load_reference_genome(
            args.reference_genome, args.low_memory)
    except AttributeError:
        args.__dict__['reference_genome'] = args.__dict__.get('reference_genome', None), None
    try:
        log('loading:', args.masking)
        args.__dict__['masking'] = args.masking, load_masking_regions(args.masking)
    except AttributeError as err:
        args.__dict__['masking'] = args.__dict__.get('masking', None), None
    try:
        if pstep == PIPELINE_STEP.PIPELINE:
            raise AttributeError()
        log('loading:', args.template_metadata)
        args.__dict__['template_metadata'] = args.template_metadata, load_templates(args.template_metadata)
    except AttributeError:
        args.__dict__['template_metadata'] = args.__dict__.get('template_metadata', None), None

    
    # decide which main function to execute
    if pstep == PIPELINE_STEP.CLUSTER:
        main_cluster(args)
    elif pstep == PIPELINE_STEP.VALIDATE:
        main_validate(args)
    elif pstep == PIPELINE_STEP.ANNOTATE:
        main_annotate(args)
    elif pstep == PIPELINE_STEP.PAIR:
        main_pairing(args)
    elif pstep == PIPELINE_STEP.SUMMARY:
        main_summary(args)
    else:  # PIPELINE
        main_pipeline(args, config)

if __name__ == '__main__':
    main()
