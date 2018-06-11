"""
This is the primary module responsible for generating svg visualizations

"""
from svgwrite import Drawing

from .elements import draw_exon_track, draw_genes, draw_template, draw_ustranscript, draw_vmarker
from .scatter import draw_scatter
from .util import generate_interval_mapping, LabelMapping

from ..annotate.genomic import IntergenicRegion
from ..interval import Interval
from ..util import DEVNULL

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'


def draw_sv_summary_diagram(
        config, ann, reference_genome=None, templates=None, ignore_absent_templates=True,
        user_friendly_labels=True, template_display_label_prefix='',
        draw_reference_transcripts=True,
        draw_reference_genes=True,
        draw_reference_templates=True,
        draw_fusion_transcript=True,
        stack_reference_transcripts=False):
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
    if not any([draw_reference_templates, draw_reference_genes, draw_reference_transcripts, draw_fusion_transcript]):
        raise AssertionError('nothing to draw')
    fusion_transcript = ann.fusion
    templates = dict() if templates is None else templates
    canvas = Drawing(size=(config.width, 1000))  # just set the height for now and change later
    labels = LabelMapping()  # keep labels consistent within the drawing
    y = config.top_margin
    x = config.left_margin

    dx_label_shift = config.label_left_margin

    x += dx_label_shift
    drawing_width = config.width - dx_label_shift - config.left_margin - config.right_margin
    # calculate the half-width for transcripts and genes etc
    half_drawing_width = (drawing_width - config.inner_margin - dx_label_shift) / 2
    second_drawing_shift = x + half_drawing_width + config.inner_margin + dx_label_shift

    if draw_reference_templates:
        try:
            template1 = templates[ann.transcript1.get_chr()]
            template2 = templates[ann.transcript2.get_chr()]
            if user_friendly_labels and template1.name:
                labels.set_key(template_display_label_prefix + template1.name, template1)
            if user_friendly_labels and template2.name:
                labels.set_key(template_display_label_prefix + template2.name, template2)

            height = [0]
            if template1 == template2:  # single template
                svg_group = draw_template(
                    config, canvas, template1, drawing_width, breakpoints=[ann.break1, ann.break2], labels=labels)
                canvas.add(svg_group)
                svg_group.translate(x, y)
                height.append(svg_group.height)
            else:  # multiple templates
                svg_group = draw_template(config, canvas, template1, half_drawing_width, breakpoints=[ann.break1], labels=labels)
                canvas.add(svg_group)
                svg_group.translate(x, y)
                height.append(svg_group.height)

                svg_group = draw_template(config, canvas, template2, half_drawing_width, breakpoints=[ann.break2], labels=labels)
                canvas.add(svg_group)
                svg_group.translate(second_drawing_shift, y)
                height.append(svg_group.height)
            y += max(height) + config.inner_margin
        except KeyError as err:
            if not ignore_absent_templates:
                raise err

    colors = dict()
    genes1 = set()
    genes2 = set()
    legend = dict()

    for gene in ann.genes_overlapping_break1:
        genes1.add(gene)
        colors[gene] = config.gene1_color

    for gene, _ in ann.genes_proximal_to_break1:
        genes1.add(gene)
        colors[gene] = config.gene1_color

    for gene in ann.genes_overlapping_break2:
        genes2.add(gene)
        colors.setdefault(gene, config.gene2_color)

    for gene, _ in ann.genes_proximal_to_break2:
        genes2.add(gene)
        colors.setdefault(gene, config.gene2_color)

    if ann.transcript1:
        try:
            genes1.add(ann.transcript1.gene)
            colors[ann.transcript1.gene] = config.gene1_color_selected
            for exon in ann.transcript1.exons:
                colors[exon] = config.exon1_color
        except AttributeError:
            genes1.add(ann.transcript1)
            colors[ann.transcript1] = config.gene1_color_selected

    if ann.transcript2:
        same = ann.transcript1 == ann.transcript2
        try:
            genes2.add(ann.transcript2.gene)
            colors[ann.transcript2.gene] = config.gene2_color_selected if not same else config.gene1_color_selected
            for exon in ann.transcript2.exons:
                colors[exon] = config.exon2_color if not same else config.exon1_color
        except AttributeError:
            genes2.add(ann.transcript2)
            colors[ann.transcript2] = config.gene2_color_selected if not same else config.gene1_color_selected

    if draw_reference_genes:
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
        for alias, gene in aliases.items():
            alias_by_gene[gene] = alias

        for gene in sorted(genes1 | genes2, key=lambda x: (str(x.get_chr()), x.start)):
            if isinstance(gene, IntergenicRegion):
                labels.add(gene, config.region_label_prefix)
            elif user_friendly_labels and not alias_failure and gene in alias_by_gene:
                labels[alias_by_gene[gene]] = gene
            else:
                labels.add(gene, config.gene_label_prefix)
        gheights = [0]

        if ann.interchromosomal:
            svg_group = draw_genes(config, canvas, genes1, half_drawing_width, [ann.break1], colors=colors, labels=labels)
            svg_group.translate(x, y)
            canvas.add(svg_group)
            gheights.append(svg_group.height)

            # second gene view
            svg_group = draw_genes(config, canvas, genes2, half_drawing_width, [ann.break2], colors=colors, labels=labels)
            svg_group.translate(second_drawing_shift, y)
            canvas.add(svg_group)
            gheights.append(svg_group.height)
        else:
            svg_group = draw_genes(
                config, canvas, genes1 | genes2, drawing_width, [ann.break1, ann.break2], colors=colors, labels=labels)
            svg_group.translate(x, y)
            canvas.add(svg_group)
            gheights.append(svg_group.height)

        y += max(gheights) + config.inner_margin

    if draw_reference_transcripts:
        theights = []
        # now the transcript level drawings
        if any([
            ann.transcript1 == ann.transcript2,
            ann.transcript1 is None,
            ann.transcript2 is None,
            isinstance(ann.transcript1, IntergenicRegion),
            isinstance(ann.transcript2, IntergenicRegion)
        ]):
            breaks = [ann.break1, ann.break2]
            transcript = ann.transcript1
            if ann.transcript1 is None or isinstance(ann.transcript1, IntergenicRegion):
                transcript = ann.transcript2
                breaks = [ann.break2]
            elif ann.transcript2 is None or isinstance(ann.transcript2, IntergenicRegion):
                breaks = [ann.break1]

            try:
                svg_group = canvas.g(class_='transcript')
                svg_group = draw_ustranscript(
                    config,
                    canvas, transcript, drawing_width,
                    breakpoints=breaks,
                    labels=labels,
                    colors=colors,
                    reference_genome=reference_genome
                )
                theights.append(svg_group.height)
                svg_group.translate(x, y)
                canvas.add(svg_group)
            except AttributeError:
                pass  # Intergenic region or None
        else:  # separate drawings
            try:
                ratio = len(ann.transcript1.exons) / (len(ann.transcript1.exons) + len(ann.transcript2.exons))
                ratio = max(0.25, min(ratio, 0.75))  # must be between 0.25 - 0.75
            except AttributeError:
                ratio = 0.5

            try:
                svg_group = canvas.g(class_='transcript')
                svg_group = draw_ustranscript(
                    config,
                    canvas, ann.transcript1,
                    half_drawing_width * 2 * ratio if not stack_reference_transcripts else drawing_width,
                    breakpoints=[ann.break1],
                    labels=labels,
                    colors=colors,
                    reference_genome=reference_genome
                )
                svg_group.translate(x, y)
                if not stack_reference_transcripts:
                    theights.append(svg_group.height)
                else:
                    y += svg_group.height + config.inner_margin
                canvas.add(svg_group)
            except AttributeError:
                pass  # Intergenic region or None

            try:
                svg_group = canvas.g(class_='transcript')
                svg_group = draw_ustranscript(
                    config,
                    canvas, ann.transcript2,
                    half_drawing_width * 2 * (1 - ratio) if not stack_reference_transcripts else drawing_width,
                    breakpoints=[ann.break2],
                    labels=labels,
                    colors=colors,
                    reference_genome=reference_genome
                )
                theights.append(svg_group.height)
                shift = second_drawing_shift - half_drawing_width + half_drawing_width * 2 * ratio
                svg_group.translate(shift if not stack_reference_transcripts else x, y)
                canvas.add(svg_group)
            except AttributeError:
                pass  # Intergenic region or None

        if theights:
            y += max(theights) + config.inner_margin

    # finally the fusion transcript level drawing
    if fusion_transcript and draw_fusion_transcript:
        for exon in fusion_transcript.exons:
            colors[exon] = config.novel_exon_color
            try:
                old_ex = fusion_transcript.exon_mapping[exon.position]
                colors[exon] = colors[old_ex]
            except KeyError:
                pass
        svg_group = canvas.g(class_='transcript')
        svg_group = draw_ustranscript(
            config,
            canvas, fusion_transcript, drawing_width,
            colors=colors,
            labels=labels,
            reference_genome=reference_genome
        )
        svg_group.translate(x, y)
        canvas.add(svg_group)
        y += svg_group.height + config.inner_margin

    y += config.bottom_margin - config.inner_margin
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


def draw_multi_transcript_overlay(config, gene, vmarkers=None, window_buffer=0, plots=None, log=DEVNULL):
    vmarkers = [] if vmarkers is None else vmarkers
    plots = [] if plots is None else plots

    canvas = Drawing(size=(config.width, 1000))  # just set the height for now and change later
    width = config.width - config.left_margin - config.right_margin - config.overlay_left_label - config.padding
    labels = LabelMapping()  # keep labels consistent within the drawing

    all_exons = set()
    colors = dict()
    for us_tx in gene.transcripts:
        for ex in us_tx.exons:
            all_exons.add(ex)
            colors[ex] = config.exon1_color if us_tx.is_best_transcript else config.exon2_color

        for spl_tx in us_tx.transcripts:
            for translation in spl_tx.translations:
                for dom in translation.domains:
                    labels.set_key(dom.name, dom.name)
    genomic_min = min([max([gene.start - window_buffer, 1])] + [m.start for m in vmarkers] + [p.xmin for p in plots if p.xmin])
    genomic_max = max([gene.end + window_buffer] + [m.end for m in vmarkers] + [p.xmax for p in plots if p.xmax])

    mapping = generate_interval_mapping(
        all_exons,
        width,
        config.exon_intron_ratio,
        config.exon_min_width,
        min_inter_width=config.min_width,
        start=genomic_min, end=genomic_max
    )
    main_group = canvas.g(class_='overlay')

    x = config.overlay_left_label + config.padding
    y = config.marker_top_margin

    for plot in plots:
        if plot.points:
            plot_group = draw_scatter(config, canvas, plot, mapping, log=log)
            main_group.add(plot_group)
            plot_group.translate(x, y)
            y += plot.height + config.padding * 2

    regular_transcripts = sorted([us_tx for us_tx in gene.transcripts if not us_tx.is_best_transcript], key=lambda x: x.name)
    for us_tx in regular_transcripts:
        group_element = draw_exon_track(
            config, canvas, us_tx, mapping,
            colors=colors,
            genomic_min=genomic_min,
            genomic_max=genomic_max
        )
        main_group.add(group_element)
        group_element.translate(x, y)

        text_element = canvas.text(
            us_tx.name,
            insert=(
                x - config.padding,
                y + config.track_height / 2 + config.font_central_shift_ratio * config.label_font_size
            ),
            fill=config.label_color,
            style=config.font_style.format(font_size=config.label_font_size, text_anchor='end'),
            class_='label'
        )
        main_group.add(text_element)
        y += config.padding + config.track_height

    best_transcripts = sorted([us_tx for us_tx in gene.transcripts if us_tx.is_best_transcript], key=lambda x: x.name)
    for us_tx in best_transcripts:
        for spl_tx in us_tx.transcripts:
            labels[us_tx.name] = spl_tx

        group_element = draw_ustranscript(config, canvas, us_tx, mapping=mapping, colors=colors, labels=labels)
        main_group.add(group_element)
        group_element.translate(x, y)

        y += config.padding + group_element.height

    y += config.marker_bottom_margin
    # now draw the breakpoints overtop
    for marker in sorted(vmarkers):
        px_itvl = Interval(
            mapping.convert_ratioed_pos(marker.start).start,
            mapping.convert_ratioed_pos(marker.end).end)
        group_element = draw_vmarker(
            config, canvas, marker, px_itvl.length(), y, label=marker.name)
        group_element.translate(x + px_itvl.start, 0)
        main_group.add(group_element)

    main_group.translate(config.left_margin, config.top_margin)
    y += config.bottom_margin
    canvas.add(main_group)
    canvas.attribs['height'] = y
    return canvas
