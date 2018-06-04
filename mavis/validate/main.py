import hashlib
import itertools
import os
import re
import subprocess
import time
import warnings

import pysam
from shortuuid import uuid

from .call import call_events
from .constants import DEFAULTS, PASS_FILENAME
from .evidence import GenomeEvidence, TranscriptomeEvidence
from ..align import align_sequences, select_contig_alignments, SUPPORTED_ALIGNER
from ..annotate.base import BioInterval
from ..annotate import file_io as _file_io
from ..bam import cigar as _cigar
from ..bam.cache import BamCache
from ..breakpoint import BreakpointPair
from ..constants import CALL_METHOD, COLUMNS, MavisNamespace, PROTOCOL
from ..util import filter_on_overlap, LOG, mkdirp, output_tabbed_file, read_inputs, write_bed_file


def main(
    inputs, output,
    bam_file, strand_specific,
    library, protocol, median_fragment_size, stdev_fragment_size, read_length,
    reference_genome, annotations, masking, aligner_reference,
    start_time=int(time.time()), **kwargs
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
        reference_genome (:class:`~mavis.annotate.file_io.ReferenceFile`): see :func:`~mavis.annotate.file_io.load_reference_genome`
        annotations (:class:`~mavis.annotate.file_io.ReferenceFile`): see :func:`~mavis.annotate.file_io.load_reference_genes`
        masking (:class:`~mavis.annotate.file_io.ReferenceFile`): see :func:`~mavis.annotate.file_io.load_masking_regions`
        aligner_reference (:class:`~mavis.annotate.file_io.ReferenceFile`): path to the aligner reference file (e.g 2bit file for blat)
    """
    mkdirp(output)
    # check the files exist early to avoid waiting for errors
    if protocol == PROTOCOL.TRANS:
        annotations.load()
    reference_genome.load()
    masking.load()

    validation_settings = {}
    validation_settings.update(DEFAULTS.items())
    validation_settings.update({k: v for k, v in kwargs.items() if k in DEFAULTS})
    validation_settings = MavisNamespace(**validation_settings)

    raw_evidence_bam = os.path.join(output, 'raw_evidence.bam')
    contig_bam = os.path.join(output, 'contigs.bam')
    evidence_bed = os.path.join(output, 'evidence.bed')

    passed_output_file = os.path.join(output, PASS_FILENAME)
    passed_bed_file = os.path.join(output, 'validation-passed.bed')
    failed_output_file = os.path.join(output, 'validation-failed.tab')
    contig_aligner_fa = os.path.join(output, 'contigs.fa')
    if validation_settings.aligner == SUPPORTED_ALIGNER.BLAT:
        contig_aligner_output = os.path.join(output, 'contigs.blat_out.pslx')
        contig_aligner_log = os.path.join(output, 'contigs.blat.log')
    elif validation_settings.aligner == SUPPORTED_ALIGNER.BWA_MEM:
        contig_aligner_output = os.path.join(output, 'contigs.bwa_mem.sam')
        contig_aligner_log = os.path.join(output, 'contigs.bwa_mem.log')
    else:
        raise NotImplementedError('unsupported aligner', validation_settings.aligner)
    igv_batch_file = os.path.join(output, 'igv.batch')
    input_bam_cache = BamCache(bam_file, strand_specific)

    bpps = read_inputs(
        inputs,
        add_default={
            COLUMNS.cluster_id: None,
            COLUMNS.stranded: False
        },
        add={
            COLUMNS.protocol: protocol,
            COLUMNS.library: library
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
                    reference_genome.content,
                    opposing_strands=bpp.opposing_strands,
                    stranded=bpp.stranded,
                    untemplated_seq=bpp.untemplated_seq,
                    data=bpp.data,
                    stdev_fragment_size=stdev_fragment_size,
                    read_length=read_length,
                    median_fragment_size=median_fragment_size,
                    **dict(validation_settings.items())
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                warnings.warn('Dropping breakpoint pair ({}) as bad input {}'.format(str(bpp), str(err)))
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            try:
                evidence = TranscriptomeEvidence(
                    annotations.content,
                    bpp.break1, bpp.break2,
                    input_bam_cache,
                    reference_genome.content,
                    opposing_strands=bpp.opposing_strands,
                    stranded=bpp.stranded,
                    untemplated_seq=bpp.untemplated_seq,
                    data=bpp.data,
                    stdev_fragment_size=stdev_fragment_size,
                    read_length=read_length,
                    median_fragment_size=median_fragment_size,
                    **dict(validation_settings.items())
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                warnings.warn('Dropping ({}) as bad input {}'.format(str(bpp), str(err)))
        else:
            raise ValueError('protocol error', bpp.data[COLUMNS.protocol])

    extended_masks = {}
    for chrom, masks in masking.content.items():  # extend masking by read length
        extended_masks[chrom] = []
        for mask in masks:
            extended_masks[chrom].append(BioInterval(
                chrom, mask.start - read_length, mask.end + read_length,
                name=mask.name
            ))

    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(evidence_clusters, extended_masks)
    contig_sequences = {}
    for i, evidence in enumerate(evidence_clusters):
        LOG()
        LOG(
            '({} of {})'.format(i + 1, len(evidence_clusters)),
            'gathered evidence for:', evidence.cluster_id,
            '' if COLUMNS.tracking_id not in evidence.data else '(tracking_id: {})'.format(evidence.tracking_id),
            time_stamp=True
        )
        LOG(evidence, time_stamp=False)
        LOG('possible event type(s):', BreakpointPair.classify(evidence), time_stamp=False)
        LOG('outer window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.outer_window1[0], evidence.outer_window1[1],
            evidence.break2.chr, evidence.outer_window2[0], evidence.outer_window2[1]), time_stamp=False)
        LOG('inner window regions:  {}:{}-{}  {}:{}-{}'.format(
            evidence.break1.chr, evidence.inner_window1[0], evidence.inner_window1[1],
            evidence.break2.chr, evidence.inner_window2[0], evidence.inner_window2[1]), time_stamp=False)
        evidence.load_evidence(log=LOG)
        LOG(
            'flanking pairs: {};'.format(len(evidence.flanking_pairs)),
            'split reads: {}, {};'.format(*[len(a) for a in evidence.split_reads]),
            'half-mapped reads: {}, {};'.format(*[len(a) for a in evidence.half_mapped]),
            'spanning-reads: {};'.format(len(evidence.spanning_reads)),
            'compatible flanking pairs:', len(evidence.compatible_flanking_pairs),
            time_stamp=False
        )
        evidence.assemble_contig(log=LOG)
        LOG('assembled {} contigs'.format(len(evidence.contigs)), time_stamp=False)
        for contig in evidence.contigs:
            name = 'seq-{}'.format(hashlib.md5(contig.seq.encode('utf-8')).hexdigest())
            LOG('>', name, '(size={}; reads={:.0f}; coverage={:.2f})'.format(
                len(contig.seq), contig.remap_score(), contig.remap_coverage()), time_stamp=False)
            LOG(contig.seq[:140], time_stamp=False)
            contig_sequences[name] = contig.seq

    LOG('will output:', contig_aligner_fa, contig_aligner_output)
    raw_contig_alignments = align_sequences(
        contig_sequences,
        input_bam_cache,
        reference_genome=reference_genome.content,
        aligner_fa_input_file=contig_aligner_fa,
        aligner_output_file=contig_aligner_output,
        clean_files=validation_settings.clean_aligner_files,
        aligner=kwargs.get('aligner', validation_settings.aligner),
        aligner_reference=aligner_reference.name[0],
        aligner_output_log=contig_aligner_log,
        blat_min_identity=kwargs.get('blat_min_identity', validation_settings.blat_min_identity),
        blat_limit_top_aln=kwargs.get('blat_limit_top_aln', validation_settings.blat_limit_top_aln),
        log=LOG
    )
    for evidence in evidence_clusters:
        select_contig_alignments(evidence, raw_contig_alignments)
    LOG('alignment complete', time_stamp=True)
    event_calls = []
    total_pass = 0
    write_bed_file(evidence_bed, itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]))
    validation_counts = {}
    for index, evidence in enumerate(evidence_clusters):
        LOG()
        LOG('({} of {}) calling events for: {} {} (tracking_id: {})'.format(
            index + 1, len(evidence_clusters), evidence.cluster_id, evidence.putative_event_types(), evidence.tracking_id), time_stamp=True)
        LOG('source:', evidence)
        calls = []
        failure_comment = None
        try:
            calls = call_events(evidence)
            event_calls.extend(calls)
        except UserWarning as err:
            LOG('warning: error in calling events', repr(err))
            failure_comment = str(err)

        if not calls:
            failure_comment = ['zero events were called'] if failure_comment is None else failure_comment
            evidence.data[COLUMNS.filter_comment] = failure_comment
            filtered_evidence_clusters.append(evidence)
        else:
            total_pass += 1

        LOG('called {} event(s)'.format(len(calls)), time_stamp=True)
        for call in calls:
            LOG(call)
            if call.call_method == CALL_METHOD.CONTIG:
                LOG('\t{} {} [{}] contig_alignment_score: {}, contig_alignment_mq: {} contig_alignment_rank: {}'.format(
                    call.event_type, call.call_method, call.contig_alignment.query_name,
                    round(call.contig_alignment.score(), 2), tuple(call.contig_alignment.mapping_quality()),
                    tuple(call.contig_alignment.alignment_rank())
                ))
                LOG('\talignment:', call.contig_alignment.alignment_id())
            elif call.contig_alignment:
                LOG('\t{} {} alignment:'.format(
                    call.event_type, call.call_method), call.contig_alignment.alignment_id())
            else:
                LOG('\t{} {}'.format(call.event_type, call.call_method), time_stamp=False)
            validation_counts[call.cluster_id] = validation_counts.get(call.cluster_id, 0) + 1
            call.data[COLUMNS.validation_id] = '{}-v{}'.format(call.cluster_id, validation_counts[call.cluster_id])
            LOG(
                '\tremapped reads: {}; spanning reads: {}; split reads: [{} ({}), {} ({}), {}]'
                ', flanking pairs: {}{}'.format(
                    0 if not call.contig else len(call.contig.input_reads),
                    len(call.spanning_reads),
                    len(call.break1_split_read_names()), len(call.break1_split_read_names(tgt=True)),
                    len(call.break2_split_read_names()), len(call.break2_split_read_names(tgt=True)),
                    len(call.linking_split_read_names()),
                    len(call.flanking_pairs),
                    '' if not call.has_compatible else '(' + str(len(call.compatible_flanking_pairs)) + ')'
                ))

    # write the output validated clusters (split by type and contig)
    for i, call in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = call.breakpoint_sequence_homology(reference_genome.content)
        except AttributeError:
            pass
        call.data.update({
            COLUMNS.break1_homologous_seq: b1_homseq,
            COLUMNS.break2_homologous_seq: b2_homseq,
        })
    LOG('{} putative calls resulted in {} events with 1 or more event call'.format(len(evidence_clusters), total_pass), time_stamp=True)
    output_tabbed_file(event_calls, passed_output_file)
    output_tabbed_file(filtered_evidence_clusters, failed_output_file)
    write_bed_file(passed_bed_file, itertools.chain.from_iterable([e.get_bed_repesentation() for e in event_calls]))

    if validation_settings.write_evidence_files:
        with pysam.AlignmentFile(contig_bam, 'wb', template=input_bam_cache.fh) as fh:
            LOG('writing:', contig_bam, time_stamp=True)
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
            LOG('writing:', raw_evidence_bam, time_stamp=True)
            reads = set()
            for evidence in evidence_clusters:
                reads.update(evidence.supporting_reads())
            for read in reads:
                read.cigar = _cigar.convert_for_igv(read.cigar)
                fh.write(read)
        # now sort the contig bam
        sort = re.sub(r'.bam$', '.sorted.bam', contig_bam)
        LOG('sorting the bam file:', contig_bam, time_stamp=True)
        pysam.sort('-o', sort, contig_bam)
        contig_bam = sort
        LOG('indexing the sorted bam:', contig_bam)
        pysam.index(contig_bam)

        # then sort the evidence bam file
        sort = re.sub(r'.bam$', '.sorted.bam', raw_evidence_bam)
        LOG('sorting the bam file:', raw_evidence_bam, time_stamp=True)
        pysam.sort('-o', sort, raw_evidence_bam)
        raw_evidence_bam = sort
        LOG('indexing the sorted bam:', raw_evidence_bam)
        pysam.index(raw_evidence_bam)

        # write the igv batch file
        with open(igv_batch_file, 'w') as fh:
            LOG('writing:', igv_batch_file, time_stamp=True)

            fh.write('load {} name="{}"\n'.format(passed_bed_file, 'passed events'))
            fh.write('load {} name="{}"\n'.format(contig_bam, 'aligned contigs'))
            fh.write('load {} name="{}"\n'.format(evidence_bed, 'evidence windows'))
            fh.write('load {} name="{}"\n'.format(raw_evidence_bam, 'raw evidence'))
            fh.write('load {} name="{} {} input"\n'.format(bam_file, library, protocol))
