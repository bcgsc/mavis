import hashlib
import itertools
import os
import re
import time
from typing import Dict, List

import pysam
from shortuuid import uuid

from ..align import SUPPORTED_ALIGNER, align_sequences, select_contig_alignments
from ..annotate.base import BioInterval
from ..annotate.file_io import ReferenceFile
from ..bam import cigar as _cigar
from ..bam.cache import BamCache
from ..breakpoint import BreakpointPair
from ..constants import CALL_METHOD, COLUMNS, PROTOCOL
from ..util import (
    filter_on_overlap,
    generate_complete_stamp,
    logger,
    mkdirp,
    output_tabbed_file,
    read_inputs,
    write_bed_file,
)
from .call import call_events
from .constants import PASS_FILENAME
from .evidence import GenomeEvidence, TranscriptomeEvidence


def main(
    inputs: List[str],
    output: str,
    library: str,
    config: Dict,
    start_time=int(time.time()),
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
        reference_genome (mavis.annotate.file_io.ReferenceFile): see :func:`mavis.annotate.file_io.load_reference_genome`
        annotations (mavis.annotate.file_io.ReferenceFile): see :func:`mavis.annotate.file_io.load_annotations`
        masking (mavis.annotate.file_io.ReferenceFile): see :func:`mavis.annotate.file_io.load_masking_regions`
        aligner_reference (mavis.annotate.file_io.ReferenceFile): path to the aligner reference file (e.g 2bit file for blat)
    """
    mkdirp(output)
    reference_genome = ReferenceFile.load_from_config(config, 'reference_genome', eager_load=True)
    annotations = ReferenceFile.load_from_config(
        config,
        'annotations',
        eager_load=bool(config['libraries'][library]['protocol'] == PROTOCOL.TRANS),
    )
    masking = ReferenceFile.load_from_config(config, 'masking')
    if not masking.is_empty():
        masking.load()

    raw_evidence_bam = os.path.join(output, 'raw_evidence.bam')
    contig_bam = os.path.join(output, 'contigs.bam')
    evidence_bed = os.path.join(output, 'evidence.bed')

    passed_output_file = os.path.join(output, PASS_FILENAME)
    passed_bed_file = os.path.join(output, 'validation-passed.bed')
    failed_output_file = os.path.join(output, 'validation-failed.tab')
    contig_aligner_fa = os.path.join(output, 'contigs.fa')
    if config['validate.aligner'] == SUPPORTED_ALIGNER.BLAT:
        contig_aligner_output = os.path.join(output, 'contigs.blat_out.pslx')
        contig_aligner_log = os.path.join(output, 'contigs.blat.log')
    elif config['validate.aligner'] == SUPPORTED_ALIGNER.BWA_MEM:
        contig_aligner_output = os.path.join(output, 'contigs.bwa_mem.sam')
        contig_aligner_log = os.path.join(output, 'contigs.bwa_mem.log')
    else:
        raise NotImplementedError('unsupported aligner', config['validate.aligner'])
    igv_batch_file = os.path.join(output, 'igv.batch')
    input_bam_cache = BamCache(
        config['libraries'][library]['bam_file'], config['libraries'][library]['strand_specific']
    )

    bpps = read_inputs(
        inputs,
        add_default={COLUMNS.cluster_id: str(uuid()), COLUMNS.stranded: False},
        overwrite={
            COLUMNS.protocol: config['libraries'][library]['protocol'],
            COLUMNS.library: library,
        },
        expand_strand=False,
        expand_orient=True,
    )
    evidence_clusters = []
    for bpp in bpps:
        if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
            try:
                evidence = GenomeEvidence(
                    bpp.break1,
                    bpp.break2,
                    input_bam_cache,
                    reference_genome.content,
                    opposing_strands=bpp.opposing_strands,
                    stranded=bpp.stranded,
                    untemplated_seq=bpp.untemplated_seq,
                    stdev_fragment_size=config['libraries'][library]['stdev_fragment_size'],
                    read_length=config['libraries'][library]['read_length'],
                    median_fragment_size=config['libraries'][library]['median_fragment_size'],
                    config=config,
                    **bpp.data,
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                logger.warning(f'Dropping breakpoint pair ({bpp}) as bad input {err}')
        elif bpp.data[COLUMNS.protocol] == PROTOCOL.TRANS:
            try:
                evidence = TranscriptomeEvidence(
                    annotations.content,
                    bpp.break1,
                    bpp.break2,
                    input_bam_cache,
                    reference_genome.content,
                    opposing_strands=bpp.opposing_strands,
                    stranded=bpp.stranded,
                    untemplated_seq=bpp.untemplated_seq,
                    stdev_fragment_size=config['libraries'][library]['stdev_fragment_size'],
                    read_length=config['libraries'][library]['read_length'],
                    median_fragment_size=config['libraries'][library]['median_fragment_size'],
                    strand_determining_read=config['libraries'][library]['strand_determining_read'],
                    config=config,
                    **bpp.data,
                )
                evidence_clusters.append(evidence)
            except ValueError as err:
                logger.warning(f'Dropping ({bpp}) as bad input {err}')
        else:
            raise ValueError('protocol error', bpp.data[COLUMNS.protocol])

    extended_masks = {}
    for chrom, masks in masking.content.items():  # extend masking by read length
        extended_masks[chrom] = []
        for mask in masks:
            extended_masks[chrom].append(
                BioInterval(
                    chrom,
                    mask.start - config['libraries'][library]['read_length'],
                    mask.end + config['libraries'][library]['read_length'],
                    name=mask.name,
                )
            )

    evidence_clusters, filtered_evidence_clusters = filter_on_overlap(
        evidence_clusters, extended_masks
    )
    contig_sequences = {}
    for i, evidence in enumerate(evidence_clusters):
        logger.info(
            f'({i + 1} of {len(evidence_clusters)}) gathered evidence for: {evidence.cluster_id}'
            + (
                ''
                if COLUMNS.tracking_id not in evidence.data
                else f' (tracking_id: {evidence.tracking_id})'
            ),
        )
        logger.info(str(evidence))
        logger.info(f'possible event type(s): {BreakpointPair.classify(evidence)}')
        logger.info(
            f'outer window regions: {evidence.break1.chr}:{evidence.outer_window1[0]}-{evidence.outer_window1[1]}  {evidence.break2.chr}:{evidence.outer_window2[0]}-{evidence.outer_window2[1]}'
        )
        logger.info(
            f'inner window regions: {evidence.break1.chr}:{evidence.inner_window1[0]}-{evidence.inner_window1[1]}  {evidence.break2.chr}:{evidence.inner_window2[0]}-{evidence.inner_window2[1]}'
        )
        evidence.load_evidence()
        logger.info(
            f'flanking pairs: {len(evidence.flanking_pairs)}'
            + '; split reads: {}, {}'.format(*[len(a) for a in evidence.split_reads])
            + '; half-mapped reads: {}, {}'.format(*[len(a) for a in evidence.half_mapped])
            + f'; spanning-reads: {len(evidence.spanning_reads)}; compatible flanking pairs: {len(evidence.compatible_flanking_pairs)}',
        )
        evidence.assemble_contig()
        logger.info(f'assembled {len(evidence.contigs)} contigs')
        for contig in evidence.contigs:
            name = 'seq-{}'.format(hashlib.md5(contig.seq.encode('utf-8')).hexdigest())
            logger.info(
                f'> {name} (size={len(contig.seq)}; reads={contig.remap_score():.0f}; coverage={contig.remap_coverage():.2f})'
            )
            logger.info(contig.seq[:140])
            contig_sequences[name] = contig.seq

    logger.info(f'will output: {contig_aligner_fa} ${contig_aligner_output}')
    raw_contig_alignments = align_sequences(
        contig_sequences,
        input_bam_cache,
        reference_genome=reference_genome.content,
        aligner_fa_input_file=contig_aligner_fa,
        aligner_output_file=contig_aligner_output,
        clean_files=config['validate.clean_aligner_files'],
        aligner=config['validate.aligner'],
        aligner_reference=config['reference.aligner_reference'][0],
        aligner_output_log=contig_aligner_log,
        blat_min_identity=config['validate.blat_min_identity'],
        blat_limit_top_aln=config['validate.blat_limit_top_aln'],
    )
    for evidence in evidence_clusters:
        select_contig_alignments(evidence, raw_contig_alignments)
    logger.info('alignment complete')
    event_calls = []
    total_pass = 0
    write_bed_file(
        evidence_bed,
        itertools.chain.from_iterable([e.get_bed_repesentation() for e in evidence_clusters]),
    )
    validation_counts = {}
    for index, evidence in enumerate(evidence_clusters):
        logger.info(
            f'({index + 1} of {len(evidence_clusters)}) calling events for: {evidence.cluster_id} {evidence.putative_event_types()} (tracking_id: {evidence.tracking_id})'
        )
        logger.info(f'source: {evidence}')
        calls = []
        failure_comment = None
        try:
            calls = call_events(evidence)
            event_calls.extend(calls)
        except UserWarning as err:
            logger.warning('error in calling events {repr(err)}')
            failure_comment = str(err)

        if not calls:
            failure_comment = (
                ['zero events were called'] if failure_comment is None else failure_comment
            )
            evidence.data[COLUMNS.filter_comment] = failure_comment
            filtered_evidence_clusters.append(evidence)
        else:
            total_pass += 1

        logger.info(f'called {len(calls)} event(s)')
        for call in calls:
            logger.info(call)
            if call.call_method == CALL_METHOD.CONTIG:
                logger.info(
                    f'{call.event_type} {call.call_method} [{call.contig_alignment.query_name}] contig_alignment_score: {round(call.contig_alignment.score(), 2)}, contig_alignment_mq: {tuple(call.contig_alignment.mapping_quality())} contig_alignment_rank: {tuple(call.contig_alignment.alignment_rank())}'
                )
                logger.info(f'alignment: {call.contig_alignment.alignment_id()}')
            elif call.contig_alignment:
                logger.info(
                    f'{call.event_type} {call.call_method} alignment: {call.contig_alignment.alignment_id()}'
                )
            else:
                logger.info('{call.event_type} {call.call_method}')
            validation_counts[call.cluster_id] = validation_counts.get(call.cluster_id, 0) + 1
            call.data[COLUMNS.validation_id] = '{}-v{}'.format(
                call.cluster_id, validation_counts[call.cluster_id]
            )
            logger.info(
                'remapped reads: {}; spanning reads: {}; split reads: [{} ({}), {} ({}), {}]'
                ', flanking pairs: {}{}'.format(
                    0 if not call.contig else len(call.contig.input_reads),
                    len(call.spanning_reads),
                    len(call.break1_split_read_names()),
                    len(call.break1_split_read_names(tgt=True)),
                    len(call.break2_split_read_names()),
                    len(call.break2_split_read_names(tgt=True)),
                    len(call.linking_split_read_names()),
                    len(call.flanking_pairs),
                    ''
                    if not call.has_compatible
                    else '(' + str(len(call.compatible_flanking_pairs)) + ')',
                )
            )

    # write the output validated clusters (split by type and contig)
    for i, call in enumerate(event_calls):
        b1_homseq = None
        b2_homseq = None
        try:
            b1_homseq, b2_homseq = call.breakpoint_sequence_homology(reference_genome.content)
        except AttributeError:
            pass
        call.data.update(
            {COLUMNS.break1_homologous_seq: b1_homseq, COLUMNS.break2_homologous_seq: b2_homseq}
        )
    logger.info(
        f'{len(evidence_clusters)} putative calls resulted in {total_pass} events with 1 or more event call'
    )
    output_tabbed_file(event_calls, passed_output_file)
    output_tabbed_file(filtered_evidence_clusters, failed_output_file)
    write_bed_file(
        passed_bed_file,
        itertools.chain.from_iterable([e.get_bed_repesentation() for e in event_calls]),
    )

    if config['validate.write_evidence_files']:
        with pysam.AlignmentFile(contig_bam, 'wb', template=input_bam_cache.fh) as fh:
            logger.info(f'writing: {contig_bam}')
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
            logger.info(f'writing: {raw_evidence_bam}')
            reads = set()
            for evidence in evidence_clusters:
                reads.update(evidence.supporting_reads())
            for read in reads:
                read.cigar = _cigar.convert_for_igv(read.cigar)
                fh.write(read)
        # now sort the contig bam
        sort = re.sub(r'.bam$', '.sorted.bam', contig_bam)
        logger.info(f'sorting the bam file: {contig_bam}')
        pysam.sort('-o', sort, contig_bam)
        contig_bam = sort
        logger.info(f'indexing the sorted bam: {contig_bam}')
        pysam.index(contig_bam)

        # then sort the evidence bam file
        sort = re.sub(r'.bam$', '.sorted.bam', raw_evidence_bam)
        logger.info(f'sorting the bam file: {raw_evidence_bam}')
        pysam.sort('-o', sort, raw_evidence_bam)
        raw_evidence_bam = sort
        logger.info(f'indexing the sorted bam: {raw_evidence_bam}')
        pysam.index(raw_evidence_bam)

        # write the igv batch file
        with open(igv_batch_file, 'w') as fh:
            logger.info(f'writing: {igv_batch_file}')

            fh.write('load {} name="{}"\n'.format(passed_bed_file, 'passed events'))
            fh.write('load {} name="{}"\n'.format(contig_bam, 'aligned contigs'))
            fh.write('load {} name="{}"\n'.format(evidence_bed, 'evidence windows'))
            fh.write('load {} name="{}"\n'.format(raw_evidence_bam, 'raw evidence'))
            fh.write(
                'load {} name="{} {} input"\n'.format(
                    config['libraries'][library]['bam_file'],
                    library,
                    config['libraries'][library]['protocol'],
                )
            )
        generate_complete_stamp(output, start_time=start_time)
