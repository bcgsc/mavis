import itertools
import os
import re
import subprocess
import uuid

import pysam

from .base import Evidence
from .call import call_events
from .constants import DEFAULTS
from .evidence import GenomeEvidence, TranscriptomeEvidence
from ..align import align_contigs
from ..annotate.base import BioInterval
from ..bam import cigar as cigar_tools
from ..bam.cache import BamCache
from ..bam.read import get_samtools_version, samtools_v0_sort, samtools_v1_sort
from ..breakpoint import BreakpointPair
from ..constants import COLUMNS, PROTOCOL
from ..util import filter_on_overlap, generate_complete_stamp, log, mkdirp, output_tabbed_file, read_inputs, write_bed_file

VALIDATION_PASS_SUFFIX = '.validation-passed.tab'


def main(
    input, output,
    bam_file, stranded_bam,
    library, protocol, median_fragment_size, stdev_fragment_size, read_length,
    reference_genome, reference_genome_filename, annotations, masking, aligner_reference,
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
        aligner_reference (str): path to the aligner reference file (e.g 2bit file for blat)
    """
    mkdirp(output)
    filename_prefix = re.sub(r'\.(txt|tsv|tab)$', '', os.path.basename(input))
    raw_evidence_bam = os.path.join(output, filename_prefix + '.raw_evidence.bam')
    contig_bam = os.path.join(output, filename_prefix + '.contigs.bam')
    evidence_bed = os.path.join(output, filename_prefix + '.evidence.bed')

    passed_output_file = os.path.join(output, filename_prefix + VALIDATION_PASS_SUFFIX)
    passed_bed_file = os.path.join(output, filename_prefix + '.validation-passed.bed')
    failed_output_file = os.path.join(output, filename_prefix + '.validation-failed.tab')
    contig_aligner_fa = os.path.join(output, filename_prefix + '.contigs.fa')
    contig_aligner_output = os.path.join(output, filename_prefix + '.contigs.blat_out.pslx')
    igv_batch_file = os.path.join(output, filename_prefix + '.igv.batch')
    input_bam_cache = BamCache(bam_file, stranded_bam)

    if samtools_version is None:
        samtools_version = get_samtools_version()

    validation_settings = {}
    validation_settings.update(DEFAULTS.__dict__)
    validation_settings.update({k: v for k, v in kwargs.items() if k in DEFAULTS.__dict__})

    bpps = read_inputs(
        [input],
        add_default={
            COLUMNS.protocol: protocol,
            COLUMNS.library: library,
            COLUMNS.cluster_id: None,
            COLUMNS.stranded: stranded_bam
        },
        expand_ns=False, explicit_strand=False,
        cast={COLUMNS.cluster_id: lambda x: str(uuid.uuid4()) if not x else x}
    )
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            evidence = GenomeEvidence(
                bpp.break1, bpp.break2,
                input_bam_cache,
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
            evidence_clusters.append(evidence)
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            evidence = TranscriptomeEvidence(
                annotations,
                bpp.break1, bpp.break2,
                input_bam_cache,
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
            evidence_clusters.append(evidence)
        else:
            raise ValueError('protocol not recognized', bpp.data[COLUMNS.protocol])

    extended_masks = {}
    for chrom, masks in masking.items():  # extend masking by read length
        extended_masks[chrom] = []
        for mask in masks:
            extended_masks[chrom].append(BioInterval(
                chrom, mask.start - read_length, mask.end + read_length,
                name=mask.name
            ))

    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(evidence_clusters, extended_masks)
    if not validation_settings['fetch_method_individual']:
        Evidence.load_multiple(evidence_clusters, log)
    for i, evidence in enumerate(evidence_clusters):
        print()
        log(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathered evidence for:', evidence.cluster_id,
            '' if 'input_id' not in evidence.data else '(input_id: {})'.format(evidence.data['input_id'])
        )
        log(evidence, time_stamp=False)
        log('possible event type(s):', BreakpointPair.classify(evidence), time_stamp=False)
        log('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.outer_window1[0], evidence.outer_window1[1],
            evidence.break2.chr, evidence.outer_window2[0], evidence.outer_window2[1]), time_stamp=False)
        log('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.inner_window1[0], evidence.inner_window1[1],
            evidence.break2.chr, evidence.inner_window2[0], evidence.inner_window2[1]), time_stamp=False)
        if validation_settings['fetch_method_individual']:
            evidence.load_evidence(log=log)
        log(
            'flanking pairs: {};'.format(len(evidence.flanking_pairs)),
            'split reads: {}, {};'.format(*[len(a) for a in evidence.split_reads]),
            'half-mapped reads: {}, {};'.format(*[len(a) for a in evidence.half_mapped]),
            'spanning-reads: {};'.format(len(evidence.spanning_reads)),
            'compatible flanking pairs:', len(evidence.compatible_flanking_pairs),
            time_stamp=False
        )
        evidence.assemble_contig(log=log)
        log('assembled {} contigs'.format(len(evidence.contigs)), time_stamp=False)
        for contig in evidence.contigs:
            log('>', contig.seq, time_stamp=False)

    log('will output:', contig_aligner_fa, contig_aligner_output)
    align_contigs(
        evidence_clusters,
        input_bam_cache,
        reference_genome=reference_genome,
        aligner_fa_input_file=contig_aligner_fa,
        aligner_output_file=contig_aligner_output,
        clean_files=False,
        aligner=kwargs.get(
            'aligner', DEFAULTS.aligner),
        aligner_reference=aligner_reference,
        blat_min_identity=kwargs.get(
            'blat_min_identity', DEFAULTS.blat_min_identity),
        blat_limit_top_aln=kwargs.get(
            'blat_limit_top_aln', DEFAULTS.blat_limit_top_aln),
        contig_aln_min_query_consumption=kwargs.get(
            'contig_aln_min_query_consumption', DEFAULTS.contig_aln_min_query_consumption),
        contig_aln_max_event_size=kwargs.get(
            'contig_aln_max_event_size', DEFAULTS.contig_aln_max_event_size),
        contig_aln_min_anchor_size=kwargs.get(
            'contig_aln_min_anchor_size', DEFAULTS.contig_aln_min_anchor_size),
        contig_aln_merge_inner_anchor=kwargs.get(
            'contig_aln_merge_inner_anchor', DEFAULTS.contig_aln_merge_inner_anchor),
        contig_aln_merge_outer_anchor=kwargs.get(
            'contig_aln_merge_outer_anchor', DEFAULTS.contig_aln_merge_outer_anchor),
    )
    log('alignment complete')
    event_calls = []
    total_pass = 0
    write_bed_file(evidence_bed, itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]))
    for index, evidence in enumerate(evidence_clusters):
        print()
        log('({} of {}) calling events for:'.format
            (index + 1, len(evidence_clusters)), evidence.cluster_id, evidence.putative_event_types())
        log('source:', evidence, time_stamp=False)
        calls = []
        failure_comment = None
        try:
            calls = call_events(evidence)
            event_calls.extend(calls)
        except UserWarning as err:
            log('warning: error in calling events', repr(err), time_stamp=False)
            failure_comment = str(err)

        if not calls:
            failure_comment = ['zero events were called'] if failure_comment is None else failure_comment
            evidence.data[COLUMNS.filter_comment] = failure_comment
            filtered_evidence_clusters.append(evidence)
        else:
            total_pass += 1

        log('called {} event(s)'.format(len(calls)))
        for i, call in enumerate(calls):
            log(call, time_stamp=False)
            log(call.event_type, call.call_method, time_stamp=False)
            call.data[COLUMNS.validation_id] = '{}-v{}'.format(call.cluster_id, i + 1)
            log(
                'remapped reads: {}; spanning reads: {}; split reads: [{} ({}), {} ({}), {}]'
                ', flanking pairs: {}{}'.format(
                    0 if not call.contig else len(call.contig.input_reads),
                    len(call.spanning_reads),
                    len(call.break1_split_reads), len(call.break1_tgt_align_split_read_names()),
                    len(call.break2_split_reads), len(call.break2_tgt_align_split_read_names()),
                    len(call.linking_split_read_names()),
                    len(call.flanking_pairs),
                    '' if not call.has_compatible else '(' + str(len(call.compatible_flanking_pairs)) + ')'
                ), time_stamp=False)

    # write the output validated clusters (split by type and contig)
    for i, call in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = call.breakpoint_sequence_homology(reference_genome)
        except AttributeError:
            pass
        call.data.update({
            COLUMNS.break1_homologous_seq: b1_homseq,
            COLUMNS.break2_homologous_seq: b2_homseq,
        })
    log('{} putative calls resulted in {} events with 1 or more event call'.format(len(evidence_clusters), total_pass))
    output_tabbed_file(event_calls, passed_output_file)
    output_tabbed_file(filtered_evidence_clusters, failed_output_file)
    write_bed_file(passed_bed_file, itertools.chain.from_iterable([e.get_bed_repesentation() for e in event_calls]))

    with pysam.AlignmentFile(contig_bam, 'wb', template=input_bam_cache.fh) as fh:
        log('writing:', contig_bam)
        for evidence in evidence_clusters:
            for contig in evidence.contigs:
                for read1, read2 in contig.alignments:
                    read1.cigar = cigar_tools.convert_for_igv(read1.cigar)
                    fh.write(read1)
                    if read2:
                        read2.cigar = cigar_tools.convert_for_igv(read2.cigar)
                        fh.write(read2)

    # write the evidence
    with pysam.AlignmentFile(raw_evidence_bam, 'wb', template=input_bam_cache.fh) as fh:
        log('writing:', raw_evidence_bam)
        reads = set()
        for evidence in evidence_clusters:
            reads.update(evidence.supporting_reads())
        for read in reads:
            read.cigar = cigar_tools.convert_for_igv(read.cigar)
            fh.write(read)
    # now sort the contig bam
    sort = re.sub('.bam$', '.sorted.bam', contig_bam)
    log('sorting the bam file:', contig_bam)
    if samtools_version <= (1, 2, 0):
        subprocess.call(samtools_v0_sort(contig_bam, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(contig_bam, sort), shell=True)
    contig_bam = sort
    log('indexing the sorted bam:', contig_bam)
    subprocess.call(['samtools', 'index', contig_bam])

    # then sort the evidence bam file
    sort = re.sub('.bam$', '.sorted.bam', raw_evidence_bam)
    log('sorting the bam file:', raw_evidence_bam)
    if samtools_version <= (1, 2, 0):
        subprocess.call(samtools_v0_sort(raw_evidence_bam, sort), shell=True)
    else:
        subprocess.call(samtools_v1_sort(raw_evidence_bam, sort), shell=True)
    raw_evidence_bam = sort
    log('indexing the sorted bam:', raw_evidence_bam)
    subprocess.call(['samtools', 'index', raw_evidence_bam])

    # write the igv batch file
    with open(igv_batch_file, 'w') as fh:
        log('writing:', igv_batch_file)
        fh.write('new\ngenome {}\n'.format(reference_genome_filename))

        fh.write('load {} name="{}"\n'.format(passed_bed_file, 'passed events'))
        fh.write('load {} name="{}"\n'.format(contig_bam, 'aligned contigs'))
        fh.write('load {} name="{}"\n'.format(evidence_bed, 'evidence windows'))
        fh.write('load {} name="{}"\n'.format(raw_evidence_bam, 'raw evidence'))
        fh.write('load {} name="{} {} input"\n'.format(bam_file, library, protocol))

    generate_complete_stamp(output, log, prefix=filename_prefix + '.')
