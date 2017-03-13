import os
import sys
import json


# local modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from .variant import annotate_events, determine_prime
from ..constants import PROTOCOL, COLUMNS, PRIME
from ..error import DrawingFitError, NotSpecifiedError
from ..illustrate.constants import DiagramSettings
from ..illustrate.constants import DEFAULTS as ILLUSTRATION_DEFAULTS
from ..illustrate.diagram import draw_sv_summary_diagram
from ..pipeline.util import output_tabbed_file, log, build_batch_id, mkdirp, read_inputs


def main(
    inputs, output,
    reference_genome, annotations, template_metadata,
    min_domain_mapping_match, min_orf_size, max_orf_cap,
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        reference_genome (object): see :func:`~mavis.annotate.file_io.load_reference_genome`
        annotations(object): see :func:`~mavis.annotate.file_io.load_reference_genes`
        template_metadata (object): see :func:`~mavis.annotate.file_io.load_templates`
        min_domain_mapping_match (float): min mapping match percent (0-1) to count a domain as mapped
        min_orf_size (int): minimum size of an :term:`open reading frame` to keep as a putative translation
        max_orf_cap (int): the maximum number of :term:`open reading frame` s to collect for any given event
    """
    DRAWINGS_DIRECTORY = os.path.join(output, 'drawings')
    TABBED_OUTPUT_FILE = os.path.join(output, 'annotations.tab')
    FA_OUTPUT_FILE = os.path.join(output, 'annotations.fusion-cdna.fa')

    mkdirp(DRAWINGS_DIRECTORY)
    # test that the sequence makes sense for a random transcript
    bpps = read_inputs(inputs, in_={COLUMNS.protocol: PROTOCOL})
    log('read {} breakpoint pairs'.format(len(bpps)))

    annotations = annotate_events(
        bpps,
        reference_genome=reference_genome, annotations=annotations,
        min_orf_size=min_orf_size, min_domain_mapping_match=min_domain_mapping_match,
        max_orf_cap=max_orf_cap,
        log=log
    )

    id_prefix = build_batch_id(prefix='annotation-', suffix='-')
    rows = []  # hold the row information for the final tsv file
    fa_sequences = {}
    # now try generating the svg
    DS = DiagramSettings(**{k: v for k, v in kwargs.items() if k in ILLUSTRATION_DEFAULTS.__dict__})

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
                    DS, ann, reference_genome=reference_genome, templates=template_metadata)

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



