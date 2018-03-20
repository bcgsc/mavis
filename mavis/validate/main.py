import hashlib
import itertools
import os
import re
import subprocess
import time
import warnings

import pysam
from shortuuid import uuid

from .base import Evidence
from .call import call_events
from .constants import DEFAULTS
from .evidence import GenomeEvidence, TranscriptomeEvidence
from ..align import align_sequences, select_contig_alignments, SUPPORTED_ALIGNER
from ..annotate.base import BioInterval
from ..bam import cigar as _cigar
from ..bam.cache import BamCache
from ..breakpoint import BreakpointPair
from ..constants import COLUMNS, MavisNamespace, PROTOCOL
from ..util import filter_on_overlap, generate_complete_stamp, log, mkdirp, output_tabbed_file, read_inputs, write_bed_file

VALIDATION_PASS_SUFFIX = '.validation-passed.tab'


def main(
    inputs, output,
    bam_file, strand_specific,
    library, protocol, median_fragment_size, stdev_fragment_size, read_length,
    reference_genome, reference_genome_filename, annotations, masking, aligner_reference,
    start_time=int(time.time()), filename_prefix='validate', **kwargs
):
    """
    Args:
        inputs (list): list of input files containing the breakpoint pairs
        output (str): path to the output directory
        bam_file (str): path the bam file
        strand_specific (bool): flag to indicate the input bam is using a strand specific protocol
        median_fragment_size (int): the median fragment size
        stdev_fragment_size (int): the standard deviation in fragment size
        read_length (int): read length
        reference_genome (Object): see :func:`~mavis.annotate.file_io.load_reference_genome`
        annotations (object): see :func:`~mavis.annotate.file_io.load_reference_genes`
        masking (object): see :func:`~mavis.annotate.file_io.load_masking_regions`
        aligner_reference (str): path to the aligner reference file (e.g 2bit file for blat)
    """
    mkdirp(output)
    validation_settings = {}
    validation_settings.update(DEFAULTS.flatten())
    validation_settings.update({k: v for k, v in kwargs.items() if k in DEFAULTS})
    validation_settings = MavisNamespace(**validation_settings)

    raw_evidence_bam = os.path.join(output, filename_prefix + '.raw_evidence.bam')
    contig_bam = os.path.join(output, filename_prefix + '.contigs.bam')
    evidence_bed = os.path.join(output, filename_prefix + '.evidence.bed')

    passed_output_file = os.path.join(output, filename_prefix + VALIDATION_PASS_SUFFIX)
    passed_bed_file = os.path.join(output, filename_prefix + '.validation-passed.bed')
    failed_output_file = os.path.join(output, filename_prefix + '.validation-failed.tab')
    contig_aligner_fa = os.path.join(output, filename_prefix + '.contigs.fa')
    if validation_settings.aligner == SUPPORTED_ALIGNER.BLAT:
        contig_aligner_output = os.path.join(output, filename_prefix + '.contigs.blat_out.pslx')
        contig_aligner_log = os.path.join(output, filename_prefix + '.contigs.blat.log')
    elif validation_settings.aligner == SUPPORTED_ALIGNER.BWA_MEM:
        contig_aligner_output = os.path.join(output, filename_prefix + '.contigs.bwa_mem.sam')
        contig_aligner_log = os.path.join(output, filename_prefix + '.contigs.bwa_mem.log')
    else:
        raise NotImplementedError('unsupported aligner', validation_settings.aligner)
    igv_batch_file = os.path.join(output, filename_prefix + '.igv.batch')
    input_bam_cache = BamCache(bam_file, strand_specific)

    bpps = read_inputs(
        inputs,
        add_default={
            COLUMNS.protocol: protocol,
            COLUMNS.library: library,
            COLUMNS.cluster_id: None,
            COLUMNS.stranded: strand_specific
        },
        expand_strand=False, expand_orient=True,
        cast={COLUMNS.cluster_id: lambda x: str(uuid()) if not x else x}
    )
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            try:
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
                    **validation_settings.flatten()
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                warnings.warn('Dropping breakpoint pair ({}) as bad input {}'.format(str(bpp), str(err)))
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            try:
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
                    **validation_settings.flatten()
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                warnings.warn('Dropping breakpoint pair ({}) as bad input {}'.format(str(bpp), str(err)))
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
    contig_sequences = {}
    for i, evidence in enumerate(evidence_clusters):
        print()
        log(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathered evidence for:', evidence.cluster_id,
            '' if COLUMNS.tracking_id not in evidence.data else '(tracking_id: {})'.format(evidence.tracking_id)
        )
        log(evidence, time_stamp=False)
        log('possible event type(s):', BreakpointPair.classify(evidence), time_stamp=False)
        log('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.outer_window1[0], evidence.outer_window1[1],
            evidence.break2.chr, evidence.outer_window2[0], evidence.outer_window2[1]), time_stamp=False)
        log('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.inner_window1[0], evidence.inner_window1[1],
            evidence.break2.chr, evidence.inner_window2[0], evidence.inner_window2[1]), time_stamp=False)
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
            name = 'seq-{}'.format(hashlib.md5(contig.seq.encode('utf-8')).hexdigest())
            log('>', name, '(size={}; reads={:.0f}; coverage={:.2f})'.format(
                len(contig.seq), contig.remap_score(), contig.remap_coverage()), time_stamp=False)
            log(contig.seq[:140], time_stamp=False)
            contig_sequences[name] = contig.seq

    log('will output:', contig_aligner_fa, contig_aligner_output)
    raw_contig_alignments = align_sequences(
        contig_sequences,
        input_bam_cache,
        reference_genome=reference_genome,
        aligner_fa_input_file=contig_aligner_fa,
        aligner_output_file=contig_aligner_output,
        clean_files=validation_settings.clean_aligner_files,
        aligner=kwargs.get('aligner', validation_settings.aligner),
        aligner_reference=aligner_reference,
        aligner_output_log=contig_aligner_log,
        blat_min_identity=kwargs.get('blat_min_identity', validation_settings.blat_min_identity),
        blat_limit_top_aln=kwargs.get('blat_limit_top_aln', validation_settings.blat_limit_top_aln),
        log=log
    )
    for evidence in evidence_clusters:
        select_contig_alignments(evidence, raw_contig_alignments)
    log('alignment complete')
    event_calls = []
    total_pass = 0
    write_bed_file(evidence_bed, itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]))
    for index, evidence in enumerate(evidence_clusters):
        print()
        log('({} of {}) calling events for: {} {} (tracking_id: {})'.format(
            index + 1, len(evidence_clusters), evidence.cluster_id, evidence.putative_event_types(), evidence.tracking_id))
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
            if call.contig_alignment:
                log('{} {} [{}] contig_alignment_score: {}, contig_alignment_mq: {} contig_alignment_rank: {}'.format(
                    call.event_type, call.call_method, call.contig_alignment.query_name,
                    round(call.contig_alignment.score(), 2), tuple(call.contig_alignment.mapping_quality()),
                    tuple(call.contig_alignment.alignment_rank())
                ), time_stamp=False)
                log('alignment: ({}, {})'.format(call.contig_alignment.read1.alignment_id,
                    None if not call.contig_alignment.read2 else call.contig_alignment.read2.alignment_id),
                    time_stamp=False)
            else:
                log(call.event_type, call.call_method, time_stamp=False)
            call.data[COLUMNS.validation_id] = '{}-v{}'.format(call.cluster_id, i + 1)
            log(
                'remapped reads: {}; spanning reads: {}; split reads: [{} ({}), {} ({}), {}]'
                ', flanking pairs: {}{}'.format(
                    0 if not call.contig else len(call.contig.input_reads),
                    len(call.spanning_reads),
                    len(call.break1_split_read_names()), len(call.break1_split_read_names(tgt=True)),
                    len(call.break2_split_read_names()), len(call.break2_split_read_names(tgt=True)),
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

    if validation_settings.write_evidence_files:
        with pysam.AlignmentFile(contig_bam, 'wb', template=input_bam_cache.fh) as fh:
            log('writing:', contig_bam)
            for evidence in evidence_clusters:
                for contig in evidence.contigs:
                    for aln in contig.alignments:
                        aln.read1.cigar = _cigar.convert_for_igv(aln.read1.cigar)
                        fh.write(aln.read1)
                        if aln.read2:
                            aln.read2.cigar = _cigar.convert_for_igv(aln.read2.cigar)
                            fh.write(aln.read2)

        # write the evidence
        with pysam.AlignmentFile(raw_evidence_bam, 'wb', template=input_bam_cache.fh) as fh:
            log('writing:', raw_evidence_bam)
            reads = set()
            for evidence in evidence_clusters:
                reads.update(evidence.supporting_reads())
            for read in reads:
                read.cigar = _cigar.convert_for_igv(read.cigar)
                fh.write(read)
        # now sort the contig bam
        sort = re.sub('.bam$', '.sorted.bam', contig_bam)
        log('sorting the bam file:', contig_bam)
        pysam.sort('-o', sort, contig_bam)
        contig_bam = sort
        log('indexing the sorted bam:', contig_bam)
        pysam.index(contig_bam)

        # then sort the evidence bam file
        sort = re.sub('.bam$', '.sorted.bam', raw_evidence_bam)
        log('sorting the bam file:', raw_evidence_bam)
        pysam.sort('-o', sort, raw_evidence_bam)
        raw_evidence_bam = sort
        log('indexing the sorted bam:', raw_evidence_bam)
        pysam.index(raw_evidence_bam)

        # write the igv batch file
        with open(igv_batch_file, 'w') as fh:
            log('writing:', igv_batch_file)
            fh.write('new\ngenome {}\n'.format(reference_genome_filename))

            fh.write('load {} name="{}"\n'.format(passed_bed_file, 'passed events'))
            fh.write('load {} name="{}"\n'.format(contig_bam, 'aligned contigs'))
            fh.write('load {} name="{}"\n'.format(evidence_bed, 'evidence windows'))
            fh.write('load {} name="{}"\n'.format(raw_evidence_bam, 'raw evidence'))
            fh.write('load {} name="{} {} input"\n'.format(bam_file, library, protocol))

    generate_complete_stamp(output, log, start_time=start_time)
