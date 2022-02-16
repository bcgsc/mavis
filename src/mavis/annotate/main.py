import hashlib
import json
import os
import time
from typing import Dict, List

from mavis_config import get_by_prefix

from ..constants import COLUMNS, PRIME, sort_columns
from ..error import DrawingFitError, NotSpecifiedError
from ..illustrate.constants import DiagramSettings
from ..illustrate.diagram import draw_sv_summary_diagram
from ..types import ReferenceGenome
from ..util import generate_complete_stamp, logger, mkdirp, read_inputs
from .constants import PASS_FILENAME
from .file_io import ReferenceFile
from .fusion import determine_prime
from .genomic import PreTranscript, Template
from .variant import (
    annotate_events,
    call_protein_indel,
    choose_more_annotated,
    choose_transcripts_by_priority,
    flatten_fusion_transcript,
    flatten_fusion_translation,
)

ACCEPTED_FILTERS = {
    'choose_more_annotated': choose_more_annotated,
    'choose_transcripts_by_priority': choose_transcripts_by_priority,
}


def draw(
    drawing_config: DiagramSettings,
    ann,
    reference_genome: ReferenceGenome,
    template_metadata: Dict[str, Template],
    drawings_directory: str,
):
    """
    produces the svg diagram and json legend for a given annotation
    """
    drawing = None
    legend = None
    initial_width = drawing_config.width

    drawing_attempts = []
    for attempt in range(0, drawing_config.max_drawing_retries):
        drawing_attempts.append(
            (initial_width + attempt * drawing_config.drawing_width_iter_increase, {})
        )
        drawing_attempts.append(
            (
                initial_width + attempt * drawing_config.drawing_width_iter_increase,
                {'stack_reference_transcripts': True},
            )
        )
    drawing_attempts.append(
        (initial_width, {'draw_fusion_transcript': False, 'draw_reference_transcripts': False})
    )

    for i, (curr_width, other_settings) in enumerate(drawing_attempts):
        logger.info(
            f'drawing attempt: {i + 1} {curr_width}px {other_settings if other_settings else ""}'
        )
        try:
            drawing_config.width = curr_width
            canvas, legend_json = draw_sv_summary_diagram(
                drawing_config,
                ann,
                reference_genome=reference_genome,
                templates=template_metadata,
                **other_settings,
            )

            gene_aliases1 = 'NA'
            gene_aliases2 = 'NA'
            try:
                if ann.transcript1.gene.aliases:
                    gene_aliases1 = '-'.join(ann.transcript1.gene.aliases)
                if ann.transcript1.is_best_transcript:
                    gene_aliases1 = 'b-' + gene_aliases1
            except AttributeError:
                pass
            try:
                if ann.transcript2.gene.aliases:
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

            name = 'mavis_{}-chr{}_chr{}-{}_{}'.format(
                ann.annotation_id, ann.break1.chr, ann.break2.chr, gene_aliases1, gene_aliases2
            )

            drawing = os.path.join(drawings_directory, name + '.svg')
            legend = os.path.join(drawings_directory, name + '.legend.json')
            logger.info(f'generating svg: {drawing}')
            canvas.saveas(drawing)

            logger.info(f'generating legend: {legend}')
            with open(legend, 'w') as fh:
                json.dump(legend_json, fh)
            break
        except DrawingFitError:
            pass
    drawing_config.width = initial_width  # reset the width
    return drawing, legend


def main(
    inputs: List[str],
    output: str,
    library: str,
    config: Dict,
    start_time=int(time.time()),
    **kwargs,
):
    """
    Args:
        inputs: list of input files to read
        output: path to the output directory
    """
    reference_genome = ReferenceFile.load_from_config(config, 'reference_genome')
    annotations = ReferenceFile.load_from_config(config, 'annotations')
    template_metadata = ReferenceFile.load_from_config(config, 'template_metadata', eager_load=True)

    drawings_directory = os.path.join(output, 'drawings')
    tabbed_output_file = os.path.join(output, PASS_FILENAME)
    fa_output_file = os.path.join(output, 'annotations.fusion-cdna.fa')

    annotation_filters = [ACCEPTED_FILTERS[a] for a in config['annotate.annotation_filters']]

    mkdirp(drawings_directory)
    # test that the sequence makes sense for a random transcript
    bpps = read_inputs(
        inputs,
        add_default={
            COLUMNS.protocol: config['libraries'][library]['protocol'],
            COLUMNS.library: library,
            COLUMNS.stranded: False,
        },
        expand_strand=False,
        expand_orient=True,
        expand_svtype=True,
    )
    logger.info(f'read {len(bpps)} breakpoint pairs')

    annotations.load()
    reference_genome.load()

    annotated_events = annotate_events(
        bpps,
        reference_genome=reference_genome.content,
        annotations=annotations.content,
        min_orf_size=config['annotate.min_orf_size'],
        min_domain_mapping_match=config['annotate.min_domain_mapping_match'],
        max_proximity=config['cluster.max_proximity'],
        max_orf_cap=config['annotate.max_orf_cap'],
        filters=annotation_filters,
    )

    # now try generating the svg
    drawing_config = DiagramSettings(**get_by_prefix(config, 'illustrate.'))

    header_req = {
        COLUMNS.break1_strand,
        COLUMNS.break2_strand,
        COLUMNS.fusion_sequence_fasta_file,
        COLUMNS.fusion_splicing_pattern,
        COLUMNS.fusion_cdna_coding_start,
        COLUMNS.fusion_cdna_coding_end,
        COLUMNS.fusion_sequence_fasta_id,
        COLUMNS.fusion_mapped_domains,
        COLUMNS.fusion_protein_hgvs,
        COLUMNS.exon_first_3prime,
        COLUMNS.exon_last_5prime,
        COLUMNS.annotation_id,
        COLUMNS.annotation_figure,
        COLUMNS.annotation_figure_legend,
        COLUMNS.cdna_synon,
        COLUMNS.protein_synon,
    }
    header = None
    logger.info(f'opening for write: {tabbed_output_file}')
    tabbed_fh = open(tabbed_output_file, 'w')
    logger.info(f'opening for write: {fa_output_file}')
    fasta_fh = open(fa_output_file, 'w')

    try:
        total = len(annotated_events)
        for i, ann in enumerate(annotated_events):
            ann_row = ann.flatten()
            ann_row[COLUMNS.fusion_sequence_fasta_file] = fa_output_file
            if header is None:
                header_req.update(ann_row.keys())
                header = sort_columns(header_req)
                tabbed_fh.write('\t'.join([str(c) for c in header]) + '\n')
            logger.info(
                f'({i + 1} of {total}) current annotation {ann.annotation_id} {ann.transcript1} {ann.transcript2} {ann.event_type}'
            )
            logger.info(str(ann))
            # get the reference sequences for either transcript
            ref_cdna_seq = {}
            ref_protein_seq = {}

            for pre_transcript in [
                x for x in [ann.transcript1, ann.transcript2] if isinstance(x, PreTranscript)
            ]:
                name = pre_transcript.name
                for spl_tx in pre_transcript.spliced_transcripts:
                    ref_cdna_seq.setdefault(spl_tx.get_seq(reference_genome.content), set()).add(
                        name
                    )
                    for translation in spl_tx.translations:
                        ref_protein_seq.setdefault(
                            translation.get_aa_seq(reference_genome.content), set()
                        ).add(name)

            # try building the fusion product
            rows = []
            cdna_synon_all = True
            # add fusion information to the current ann_row
            for spl_fusion_tx in [] if not ann.fusion else ann.fusion.transcripts:
                seq = ann.fusion.get_cdna_seq(spl_fusion_tx.splicing_pattern)
                # make the fasta id a hex of the string to avoid having to load the sequences later
                fusion_fa_id = 'seq-{}'.format(hashlib.md5(seq.encode('utf-8')).hexdigest())
                fasta_fh.write('> {}\n{}\n'.format(fusion_fa_id, seq))
                cdna_synon = ';'.join(sorted(list(ref_cdna_seq.get(seq, set()))))

                temp_row = {}
                temp_row.update(ann_row)
                temp_row.update(flatten_fusion_transcript(spl_fusion_tx))
                temp_row[COLUMNS.fusion_sequence_fasta_id] = fusion_fa_id
                temp_row[COLUMNS.cdna_synon] = cdna_synon if cdna_synon else None
                if not cdna_synon:
                    cdna_synon_all = False
                if spl_fusion_tx.translations:
                    # duplicate the ann_row for each translation
                    for fusion_translation in spl_fusion_tx.translations:
                        nrow = dict()
                        nrow.update(ann_row)
                        nrow.update(temp_row)
                        aa_seq = fusion_translation.get_aa_seq()
                        protein_synon = ';'.join(sorted(list(ref_protein_seq.get(aa_seq, set()))))
                        nrow[COLUMNS.protein_synon] = protein_synon if protein_synon else None
                        # select the exon
                        nrow.update(flatten_fusion_translation(fusion_translation))
                        if ann.single_transcript() and ann.transcript1.translations:
                            nrow[COLUMNS.fusion_protein_hgvs] = call_protein_indel(
                                ann.transcript1.translations[0],
                                fusion_translation,
                                reference_genome.content,
                            )
                        rows.append(nrow)
                else:
                    temp_row.update(ann_row)
                    rows.append(temp_row)
            # draw the annotation and add the path to all applicable rows (one drawing for multiple annotated_events)
            if any(
                [
                    not ann.fusion and not config['annotate.draw_fusions_only'],
                    ann.fusion and not config['annotate.draw_non_synonymous_cdna_only'],
                    ann.fusion
                    and config['annotate.draw_non_synonymous_cdna_only']
                    and not cdna_synon_all,
                ]
            ):
                drawing, legend = draw(
                    drawing_config,
                    ann,
                    reference_genome.content,
                    template_metadata.content,
                    drawings_directory,
                )
                for row in rows + [ann_row]:
                    row[COLUMNS.annotation_figure] = drawing
                    row[COLUMNS.annotation_figure_legend] = legend
            if not rows:
                rows = [ann_row]
            for row in rows:
                tabbed_fh.write('\t'.join([str(row.get(k, None)) for k in header]) + '\n')
        generate_complete_stamp(output, start_time=start_time)
    finally:
        logger.info(f'closing: {tabbed_output_file}')
        tabbed_fh.close()
        logger.info(f'closing: {fa_output_file}')
        fasta_fh.close()
