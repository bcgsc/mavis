import os
import pysam
import re
import subprocess
import uuid
import itertools
from ..bam.cache import BamCache
from ..blat import blat_contigs
from ..breakpoint import BreakpointPair
from ..constants import PROTOCOL, COLUMNS
from .call import call_events
from .evidence import GenomeEvidence, TranscriptomeEvidence
from .base import Evidence
from .constants import DEFAULTS
from ..annotate.base import BioInterval
from ..bam.read import get_samtools_version, samtools_v0_sort, samtools_v1_sort
from ..bam import cigar as cigar_tools
from ..util import read_inputs, log, output_tabbed_file, filter_on_overlap, write_bed_file
from ..util import generate_complete_stamp

VALIDATION_PASS_SUFFIX = '.validation-passed.tab'


def main(
    input, output,
    bam_file, stranded_bam,
    library, protocol, median_fragment_size, stdev_fragment_size, read_length,
    reference_genome, reference_genome_filename, annotations, masking, blat_2bit_reference,
    samtools_version, **kwargs
):
    """
    Args:
        input (str): the input file containing the breakpoint pairs
        output (str): path to the output directory
        bam_file (str): path the bam file
        stranded_bam (bool): flag to indicate the input bam is using a strand specific protocol
        median_fragment_size (int): the median fragment size
        stdev_fragment_size (int): the standard deviation in fragment size
        read_length (int): read length
        reference_genome (Object): see :func:`~mavis.annotate.file_io.load_reference_genome`
        annotations (object): see :func:`~mavis.annotate.file_io.load_reference_genes`
        masking (object): see :func:`~mavis.annotate.file_io.load_masking_regions`
        blat_2bit_reference (str): path to the 2bit reference file
    """
    FILENAME_PREFIX = re.sub('\.(txt|tsv|tab)$', '', os.path.basename(input))
    RAW_EVIDENCE_BAM = os.path.join(output, FILENAME_PREFIX + '.raw_evidence.bam')
    CONTIG_BAM = os.path.join(output, FILENAME_PREFIX + '.contigs.bam')
    EVIDENCE_BED = os.path.join(output, FILENAME_PREFIX + '.evidence.bed')

    PASSED_OUTPUT_FILE = os.path.join(output, FILENAME_PREFIX + VALIDATION_PASS_SUFFIX)
    PASSED_BED_FILE = os.path.join(output, FILENAME_PREFIX + '.validation-passed.bed')
    FAILED_OUTPUT_FILE = os.path.join(output, FILENAME_PREFIX + '.validation-failed.tab')
    CONTIG_BLAT_FA = os.path.join(output, FILENAME_PREFIX + '.contigs.fa')
    CONTIG_BLAT_OUTPUT = os.path.join(output, FILENAME_PREFIX + '.contigs.blat_out.pslx')
    IGV_BATCH_FILE = os.path.join(output, FILENAME_PREFIX + '.igv.batch')
    INPUT_BAM_CACHE = BamCache(bam_file, stranded_bam)

    if samtools_version is None:
        samtools_version = get_samtools_version()

    validation_settings = {}
    validation_settings.update(DEFAULTS.__dict__)
    validation_settings.update({k: v for k, v in kwargs.items() if k in DEFAULTS.__dict__})
    evidence_reads = set()  # keep track of collected reads to use for ouput

    split_read_contigs = set()
    chr_to_index = {}
    bpps = read_inputs(
        [input], add={COLUMNS.protocol: protocol, COLUMNS.library: library, COLUMNS.cluster_id: None},
        expand_ns=False, explicit_strand=False,
        cast={COLUMNS.cluster_id: lambda x: str(uuid.uuid4()) if not x else x}
    )
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            e = GenomeEvidence(
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                reference_genome,
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=stdev_fragment_size,
                read_length=read_length,
                median_fragment_size=median_fragment_size,
                **validation_settings
            )
            evidence_clusters.append(e)
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            e = TranscriptomeEvidence(
                annotations,
                bpp.break1, bpp.break2,
                INPUT_BAM_CACHE,
                reference_genome,
                opposing_strands=bpp.opposing_strands,
                stranded=bpp.stranded,
                untemplated_seq=bpp.untemplated_seq,
                data=bpp.data,
                stdev_fragment_size=stdev_fragment_size,
                read_length=read_length,
                median_fragment_size=median_fragment_size,
                **validation_settings
            )
            evidence_clusters.append(e)
        else:
            raise ValueError('protocol not recognized', bpp.data[COLUMNS.protocol])

    extended_masks = {}
    for chr, masks in masking.items():  # extend masking by read length
        extended_masks[chr] = []
        for mask in masks:
            extended_masks[chr].append(BioInterval(
                chr, mask.start - read_length, mask.end + read_length,
                name=mask.name
            ))

    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(evidence_clusters, extended_masks)
    if not validation_settings['fetch_method_individual']:
        Evidence.load_multiple(evidence_clusters, log)
    for i, e in enumerate(evidence_clusters):
        print()
        log(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathered evidence for:', e.cluster_id
        )
        log(e, time_stamp=False)
        log('possible event type(s):', BreakpointPair.classify(e), time_stamp=False)
        log('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.outer_window1[0], e.outer_window1[1],
            e.break2.chr, e.outer_window2[0], e.outer_window2[1]), time_stamp=False)
        log('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            e.break1.chr, e.inner_window1[0], e.inner_window1[1],
            e.break2.chr, e.inner_window2[0], e.inner_window2[1]), time_stamp=False)
        if validation_settings['fetch_method_individual']:
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
        for contig in e.contigs:
            log('>', contig.seq, time_stamp=False)

    log('will output:', CONTIG_BLAT_FA, CONTIG_BLAT_OUTPUT)
    blat_contigs(
        evidence_clusters,
        INPUT_BAM_CACHE,
        reference_genome=reference_genome,
        blat_2bit_reference=blat_2bit_reference,
        blat_fa_input_file=CONTIG_BLAT_FA,
        blat_pslx_output_file=CONTIG_BLAT_OUTPUT,
        clean_files=False,
        blat_min_percent_of_max_score=kwargs.get(
            'blat_min_percent_of_max_score', DEFAULTS.blat_min_percent_of_max_score),
        blat_min_identity=kwargs.get(
            'blat_min_identity', DEFAULTS.blat_min_identity),
        contig_aln_min_query_consumption=kwargs.get(
            'contig_aln_min_query_consumption', DEFAULTS.contig_aln_min_query_consumption),
        contig_aln_max_event_size=kwargs.get(
            'contig_aln_max_event_size', DEFAULTS.contig_aln_max_event_size),
        contig_aln_min_anchor_size=kwargs.get(
            'contig_aln_min_anchor_size', DEFAULTS.contig_aln_min_anchor_size),
        contig_aln_merge_inner_anchor=kwargs.get(
            'contig_aln_merge_inner_anchor', DEFAULTS.contig_aln_merge_inner_anchor),
        contig_aln_merge_outer_anchor=kwargs.get(
            'contig_aln_merge_outer_anchor', DEFAULTS.contig_aln_merge_outer_anchor)
    )
    log('alignment complete')
    event_calls = []
    write_bed_file(EVIDENCE_BED, itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]))
    for index, e in enumerate(evidence_clusters):
        print()
        log('({} of {}) calling events for:'.format
            (index + 1, len(evidence_clusters)), e.cluster_id, e.putative_event_types())
        log('source:', e, time_stamp=False)
        calls = []
        failure_comment = None
        try:
            calls = call_events(e)
            event_calls.extend(calls)
        except UserWarning as err:
            log('warning: error in calling events', repr(err), time_stamp=False)
            failure_comment = str(err)

        if len(calls) == 0:
            failure_comment = ['zero events were called'] if failure_comment is None else failure_comment
            e.data[COLUMNS.filter_comment] = failure_comment
            filtered_evidence_clusters.append(e)

        log('called {} event(s)'.format(len(calls)))
        for i, ev in enumerate(calls):
            log(ev, time_stamp=False)
            log(ev.event_type, ev.call_method, time_stamp=False)
            ev.data[COLUMNS.validation_id] = '{}-v{}'.format(ev.cluster_id, i + 1)
            log('remapped reads: {}; spanning reads: {}; split reads: [{} ({}), {} ({}), {}], flanking pairs: {}'.format(
                0 if not ev.contig else len(ev.contig.input_reads),
                len(ev.spanning_reads),
                len(ev.break1_split_reads), len(ev.break1_tgt_align_split_read_names()),
                len(ev.break2_split_reads), len(ev.break2_tgt_align_split_read_names()),
                len(ev.linking_split_read_names()),
                len(ev.flanking_pairs)), time_stamp=False)

    # write the output validated clusters (split by type and contig)
    for i, ec in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = ec.breakpoint_sequence_homology(reference_genome)
        except AttributeError:
            pass
        ec.data.update({
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
    if samtools_version <= (1, 2, 0):
        subprocess.call(samtools_v0_sort(CONTIG_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(CONTIG_BAM, sort), shell=True)
    CONTIG_BAM = sort
    log('indexing the sorted bam:', CONTIG_BAM)
    subprocess.call(['samtools', 'index', CONTIG_BAM])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted.bam', RAW_EVIDENCE_BAM)
    log('sorting the bam file:', RAW_EVIDENCE_BAM)
    if samtools_version <= (1, 2, 0):
        subprocess.call(samtools_v0_sort(RAW_EVIDENCE_BAM, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(RAW_EVIDENCE_BAM, sort), shell=True)
    RAW_EVIDENCE_BAM = sort
    log('indexing the sorted bam:', RAW_EVIDENCE_BAM)
    subprocess.call(['samtools', 'index', RAW_EVIDENCE_BAM])

    # write the igv batch file
    with open(IGV_BATCH_FILE, 'w') as fh:
        log('writing:', IGV_BATCH_FILE)
        fh.write('new\ngenome {}\n'.format(reference_genome_filename))

        fh.write('load {} name="{}"\n'.format(PASSED_BED_FILE, 'passed events'))
        fh.write('load {} name="{}"\n'.format(CONTIG_BAM, 'aligned contigs'))
        fh.write('load {} name="{}"\n'.format(EVIDENCE_BED, 'evidence windows'))
        fh.write('load {} name="{}"\n'.format(RAW_EVIDENCE_BAM, 'raw evidence'))
        fh.write('load {} name="{} {} input"\n'.format(bam_file, library, protocol))

    generate_complete_stamp(output, log, prefix=FILENAME_PREFIX + '.')
