"""
This is the primary module responsible for generating svg visualizations

"""
from ..annotate.genomic import IntergenicRegion
from ..annotate.variant import FusionTranscript
from ..constants import STRAND, ORIENT, CODON_SIZE, GIESMA_STAIN
from ..error import DrawingFitError, NotSpecifiedError
from ..interval import Interval
from .elements import *
from .scatter import ScatterPlot, draw_scatter
from .util import *
from colour import Color
from svgwrite import Drawing
import re
import svgwrite

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'


def draw_sv_summary_diagram(
    DS, ann, fusion_transcript=None, reference_genome=None, templates=None, ignore_absent_templates=True,
    show_template=True, user_friendly_labels=True, template_display_label_prefix='c'
):
    """
    this is the main drawing function. It decides between layouts
    where each view-level is split into one or two diagrams (side-by-side)
    dependant on whether the event is interchromosomal, within a single
    transcript, etc.

    Diagrams have four levels
        - template
        - gene
        - transcript
        - fusion transcript/translation

    Args:
        ann (Annotation): the annotation object to be illustrated
        fusion_transcript (FusionTranscript): the fusion transcript built from this annotation
        reference_genome (dict of str by str): reference sequences
        templates (list of Template): list of templates, used in drawing the template-level view
        ignore_absent_templates (bool):
            if true then will not raise an error if the template information is not given but will
            not draw the template instead
        show_template (bool): if false the template-level view is not drawn
        user_friendly_labels (bool):
            if True, genes are labelled by their aliases (where possible) and domains are labeled by their
            names (where possible)
        template_display_label_prefix (str): the character to precede the template label
    """
    templates = dict() if templates is None else templates
    canvas = Drawing(size=(DS.WIDTH, 1000))  # just set the height for now and change later
    labels = LabelMapping()  # keep labels consistent within the drawing
    y = DS.TOP_MARGIN
    x = DS.LEFT_MARGIN

    dx_label_shift = DS.LABEL_LEFT_MARGIN

    x += dx_label_shift
    drawing_width = DS.WIDTH - dx_label_shift - DS.LEFT_MARGIN - DS.RIGHT_MARGIN
    # calculate the half-width for transcripts and genes etc
    half_drawing_width = (drawing_width - DS.INNER_MARGIN - dx_label_shift) / 2
    second_drawing_shift = x + half_drawing_width + DS.INNER_MARGIN + dx_label_shift

    if show_template:
        try:
            template1 = templates[ann.transcript1.get_chr()]
            template2 = templates[ann.transcript2.get_chr()]
            if user_friendly_labels and template1.name:
                labels.set_key(template_display_label_prefix + template1.name, template1)
            if user_friendly_labels and template2.name:
                labels.set_key(template_display_label_prefix + template2.name, template2)

            h = [0]
            if template1 == template2:  # single template
                g = draw_template(
                    DS, canvas, template1, drawing_width, breakpoints=[ann.break1, ann.break2], labels=labels)
                canvas.add(g)
                g.translate(x, y)
                h.append(g.height)
            else:  # multiple templates
                g = draw_template(DS, canvas, template1, half_drawing_width, breakpoints=[ann.break1], labels=labels)
                canvas.add(g)
                g.translate(x, y)
                h.append(g.height)

                g = draw_template(DS, canvas, template2, half_drawing_width, breakpoints=[ann.break2], labels=labels)
                canvas.add(g)
                g.translate(second_drawing_shift, y)
                h.append(g.height)
            y += max(h)
        except KeyError as err:
            if not ignore_absent_templates:
                raise err

    colors = dict()
    genes = set()
    genes1 = set()
    genes2 = set()
    legend = dict()

    for g in ann.genes_overlapping_break1:
        genes1.add(g)
        colors[g] = DS.GENE1_COLOR

    for g, d in ann.genes_proximal_to_break1:
        genes1.add(g)
        colors[g] = DS.GENE1_COLOR

    for g in ann.genes_overlapping_break2:
        genes2.add(g)
        colors.setdefault(g, DS.GENE2_COLOR)

    for g, d in ann.genes_proximal_to_break2:
        genes2.add(g)
        colors.setdefault(g, DS.GENE2_COLOR)

    if ann.transcript1:
        try:
            genes1.add(ann.transcript1.gene)
            colors[ann.transcript1.gene] = DS.GENE1_COLOR_SELECTED
            for e in ann.transcript1.exons:
                colors[e] = DS.EXON1_COLOR
        except AttributeError:
            genes1.add(ann.transcript1)
            colors[ann.transcript1] = DS.GENE1_COLOR_SELECTED

    if ann.transcript2:
        try:
            genes2.add(ann.transcript2.gene)
            colors.setdefault(ann.transcript2.gene, DS.GENE2_COLOR_SELECTED)
            for e in ann.transcript2.exons:
                colors.setdefault(e, DS.EXON2_COLOR)
        except AttributeError:
            genes2.add(ann.transcript2)
            colors[ann.transcript2] = DS.GENE2_COLOR_SELECTED

    # set all the labels so that they are re-used correctly
    aliases = {}
    alias_failure = False
    for gene in sorted(genes1 | genes2, key=lambda x: (str(x.get_chr()), x.start)):
        if alias_failure:
            break
        try:
            for alias in gene.aliases:
                if alias in aliases and aliases[alias] != gene:
                    alias_failure = True
                    break
            if len(gene.aliases) == 1:
                aliases[gene.aliases[0]] = gene
            else:
                for alias in gene.aliases:  # can't label when multiple
                    aliases[alias] = None
        except AttributeError:
            pass

    alias_by_gene = {}
    for k, v in aliases.items():
        alias_by_gene[v] = k

    for gene in sorted(genes1 | genes2, key=lambda x: (str(x.get_chr()), x.start)):
        if isinstance(gene, IntergenicRegion):
            l = labels.add(gene, DS.REGION_LABEL_PREFIX)
        elif user_friendly_labels and not alias_failure and gene in alias_by_gene:
            labels[alias_by_gene[gene]] = gene
        else:
            l = labels.add(gene, DS.GENE_LABEL_PREFIX)

    gheights = [0]

    if ann.interchromosomal:
        g = draw_genes(DS, canvas, genes1, half_drawing_width, [ann.break1], colors=colors, labels=labels)
        g.translate(x, y)
        canvas.add(g)
        gheights.append(g.height)

        # second gene view
        g = draw_genes(DS, canvas, genes2, half_drawing_width, [ann.break2], colors=colors, labels=labels)
        g.translate(second_drawing_shift, y)
        canvas.add(g)
        gheights.append(g.height)
    else:
        g = draw_genes(DS,
            canvas, genes1 | genes2, drawing_width, [ann.break1, ann.break2], colors=colors, labels=labels)
        g.translate(x, y)
        canvas.add(g)
        gheights.append(g.height)

    y += max(gheights) + DS.INNER_MARGIN

    theights = [0]
    # now the transcript level drawings
    if ann.transcript1 == ann.transcript2:
        try:
            g = canvas.g(class_='transcript')
            g = draw_ustranscript(DS,
                canvas, ann.transcript1, drawing_width,
                breakpoints=[ann.break1, ann.break2],
                labels=labels,
                colors=colors,
                reference_genome=reference_genome
            )
            theights.append(g.height)
            g.translate(x, y)
            canvas.add(g)
        except AttributeError:
            pass  # Intergenic region or None
    else:  # separate drawings
        try:
            g = canvas.g(class_='transcript')
            g = draw_ustranscript(DS,
                canvas, ann.transcript1, half_drawing_width,
                breakpoints=[ann.break1],
                labels=labels,
                colors=colors,
                reference_genome=reference_genome
            )
            theights.append(g.height)
            g.translate(x, y)
            canvas.add(g)
        except AttributeError:
            pass  # Intergenic region or None

        try:
            g = canvas.g(class_='transcript')
            g = draw_ustranscript(DS,
                canvas, ann.transcript2, half_drawing_width,
                breakpoints=[ann.break2],
                labels=labels,
                colors=colors,
                reference_genome=reference_genome
            )
            theights.append(g.height)
            g.translate(second_drawing_shift, y)
            canvas.add(g)
        except AttributeError:
            pass  # Intergenic region or None

    y += max(theights)
    if max(theights) == 0:
        y -= DS.INNER_MARGIN

    # finally the fusion transcript level drawing
    if fusion_transcript:
        y += DS.INNER_MARGIN
        for exon in fusion_transcript.exons:
            old_ex = fusion_transcript.exon_mapping[exon.position]
            if old_ex in colors:
                colors[exon] = colors[old_ex]
        g = canvas.g(class_='transcript')
        g = draw_ustranscript(DS,
            canvas, fusion_transcript, drawing_width,
            colors=colors,
            labels=labels,
            reference_genome=reference_genome
        )
        g.translate(x, y)
        canvas.add(g)
        y += g.height

    y += DS.BOTTOM_MARGIN
    canvas.attribs['height'] = y
    for label, obj in labels.items():
        if label in legend:
            continue
        try:
            legend[label] = obj.to_dict()
        except AttributeError:
            legend[label] = str(obj)
    # now make the json legend
    return canvas, legend


def draw_multi_transcript_overlay(DS, gene, vmarkers=None, window_buffer=0, plots=None):
    vmarkers = [] if vmarkers is None else vmarkers
    plots = [] if plots is None else plots

    canvas = Drawing(size=(DS.WIDTH, 1000))  # just set the height for now and change later
    w = DS.WIDTH - DS.LEFT_MARGIN - DS.RIGHT_MARGIN - DS.OVERLAY_LEFT_LABEL - DS.PADDING
    labels = LabelMapping()  # keep labels consistent within the drawing

    all_exons = set()
    colors = dict()
    for tx in gene.transcripts:
        for ex in tx.exons:
            all_exons.add(ex)
            colors[ex] = DS.EXON1_COLOR if tx.is_best_transcript else DS.EXON2_COLOR

        for tr in tx.transcripts:
            for tl in tr.translations:
                for dom in tl.domains:
                    labels.set_key(dom.name, dom.name)

    st = min([max([gene.start - window_buffer, 1])] + [m.start for m in vmarkers] + [p.xmin for p in plots])
    end = max([gene.end + window_buffer] + [m.end for m in vmarkers] + [p.xmax for p in plots])

    mapping = generate_interval_mapping(
        all_exons,
        w,
        DS.EXON_INTRON_RATIO,
        DS.EXON_MIN_WIDTH,
        min_inter_width=DS.MIN_WIDTH,
        start=st, end=end
    )
    main_group = canvas.g(class_='overlay')

    x = DS.OVERLAY_LEFT_LABEL + DS.PADDING
    y = DS.MARKER_TOP_MARGIN

    for plot in plots:
        plot_group = draw_scatter(DS, canvas, plot, mapping)
        main_group.add(plot_group)
        plot_group.translate(x, y)
        y += plot.height + DS.PADDING * 2

    regular_transcripts = sorted([tx for tx in gene.transcripts if not tx.is_best_transcript], key=lambda x: x.name)
    for tx in regular_transcripts:
        g = draw_exon_track(DS, canvas, tx, mapping, colors=colors)
        main_group.add(g)
        g.translate(x, y)

        t = canvas.text(
            tx.name,
            insert=(
                x - DS.PADDING,
                y + DS.TRACK_HEIGHT / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE
            ),
            fill=DS.LABEL_COLOR,
            style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='end'),
            class_='label'
        )
        main_group.add(t)
        y += DS.PADDING + DS.TRACK_HEIGHT

    best_transcripts = sorted([tx for tx in gene.transcripts if tx.is_best_transcript], key=lambda x: x.name)
    for tx in best_transcripts:
        for txx in tx.transcripts:
            labels[tx.name] = txx

        g = draw_ustranscript(DS, canvas, tx, mapping=mapping, colors=colors, labels=labels)
        main_group.add(g)
        g.translate(x, y)

        y += DS.PADDING + g.height

    y += DS.MARKER_BOTTOM_MARGIN
    # now draw the breakpoints overtop
    for i, m in enumerate(sorted(vmarkers)):
        s = Interval.convert_ratioed_pos(mapping, m.start)
        t = Interval.convert_ratioed_pos(mapping, m.end)
        px_itvl = Interval(s.start, t.end)
        bg = draw_vmarker(DS,
            canvas, m, px_itvl.length(), y, label=labels.add(m, DS.MARKER_LABEL_PREFIX))
        bg.translate(x + px_itvl.start, 0)
        main_group.add(bg)

    main_group.translate(DS.LEFT_MARGIN, DS.TOP_MARGIN)
    y += DS.BOTTOM_MARGIN
    canvas.add(main_group)
    canvas.attribs['height'] = y
    return canvas
