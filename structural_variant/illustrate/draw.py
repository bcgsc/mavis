"""
.. todo::

    add optional scatter plots to subdiagrams for cna and expression data (possibly using matplotlib?)

.. todo::

    add mask to highlight the region not kept by the breakpoints (for 2 of more breakpoints only)
    can use the css property "pointer-events: none;" to avoid issues with interactivity

"""
import svgwrite
import re
from svgwrite import Drawing
from ..interval import Interval
from ..constants import STRAND, ORIENT, CODON_SIZE, GIESMA_STAIN
from colour import Color
from ..error import DrawingFitError, NotSpecifiedError
from ..annotate.genomic import IntergenicRegion
from ..annotate.variant import FusionTranscript
from .labels import dynamic_label_color, LabelMapping
from .scatter import ScatterPlot, draw_scatter

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'


class Tag(svgwrite.base.BaseElement):
    def __init__(DS, elementname, content='', **kwargs):
        DS.elementname = elementname
        super(Tag, DS).__init__(**kwargs)
        DS.content = content

    def get_xml(DS):
        xml = super(Tag, DS).get_xml()
        xml.text = DS.content
        return xml


def draw_legend(DS, canvas, swatches, border=True):
    main_group = canvas.g(class_='legend')
    y = DS.PADDING if border else 0
    x = DS.PADDING if border else 0
    for swatch, label in swatches:
        g = canvas.g()
        g.add(canvas.rect(
            (0, 0),
            (DS.LEGEND_SWATCH_SIZE, DS.LEGEND_SWATCH_SIZE),
            fill=swatch,
            stroke=DS.LEGEND_SWATCH_STROKE
        ))

        g.add(canvas.text(
            label,
            insert=(DS.LEGEND_SWATCH_SIZE + DS.PADDING, DS.LEGEND_SWATCH_SIZE / 2),
            fill=DS.LEGEND_FONT_COLOR,
            style=DS.FONT_STYLE.format(text_anchor='start', font_size=DS.LEGEND_FONT_SIZE),
            class_='label'
        ))
        g.translate(x, y)
        main_group.add(g)
        y += DS.LEGEND_SWATCH_SIZE + DS.PADDING

    w = max([len(l) for c, l in swatches]) * DS.LEGEND_FONT_SIZE * DS.FONT_WIDTH_HEIGHT_RATIO + \
        DS.PADDING * (3 if border else 1) + DS.LEGEND_SWATCH_SIZE

    if border:
        main_group.add(canvas.rect(
            (0, 0), (w, y), fill='none', stroke=DS.LEGEND_BORDER_STROKE,
            stroke_width=DS.LEGEND_BORDER_STROKE_WIDTH
        ))
    else:
        y -= DS.PADDING
    setattr(main_group, 'height', y)
    setattr(main_group, 'width', w)
    setattr(main_group, 'labels', None)
    setattr(main_group, 'mapping', None)
    return main_group


def draw(
    DS, ann, fusion_transcript=None, REFERENCE_GENOME=None, templates=None, ignore_absent_templates=True,
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
        REFERENCE_GENOME (dict of str by str): reference sequences
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
    print('draw(', DS, ann, fusion_transcript)
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
                REFERENCE_GENOME=REFERENCE_GENOME
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
                REFERENCE_GENOME=REFERENCE_GENOME
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
                REFERENCE_GENOME=REFERENCE_GENOME
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
            REFERENCE_GENOME=REFERENCE_GENOME
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


def _draw_exon_track(DS, canvas, transcript, mapping, colors=None, x_start=None, x_end=None, translation=None):
    """
    """
    print('_draw_exon_track', repr(transcript), transcript.exons, transcript.gene)
    colors = {} if colors is None else colors
    main_group = canvas.g(class_='exon_track')

    y = DS.TRACK_HEIGHT / 2
    exons = sorted(transcript.exons, key=lambda x: x.start)

    s = Interval.convert_ratioed_pos(mapping, exons[0].start).start if x_start is None else x_start
    t = Interval.convert_ratioed_pos(mapping, exons[-1].end).end if x_end is None else x_end

    main_group.add(
        canvas.rect(
            (s, y - DS.SCAFFOLD_HEIGHT / 2), (t - s + 1, DS.SCAFFOLD_HEIGHT),
            fill=DS.SCAFFOLD_COLOR, class_='scaffold'
        ))

    # draw the exons
    for exon in exons:
        s = Interval.convert_ratioed_pos(mapping, exon.start).start
        t = Interval.convert_ratioed_pos(mapping, exon.end).end
        pxi = Interval(s, t)
        c = colors.get(exon, DS.EXON1_COLOR)
        group = draw_exon(DS,
            canvas, exon, pxi.length(), DS.TRACK_HEIGHT, c,
            label=transcript.exon_number(exon),
            translation=translation
        )
        group.translate(pxi.start, y - DS.TRACK_HEIGHT / 2)
        main_group.add(group)

    setattr(main_group, 'height', y + DS.TRACK_HEIGHT / 2)
    setattr(main_group, 'width', t - s + 1)
    return main_group


def _draw_transcript_with_translation(
    DS, canvas, translation, labels, colors, mapping, REFERENCE_GENOME=None, x_start=None, x_end=None
):
    main_group = canvas.g()
    ust = translation.transcript.reference_object
    tr = translation.transcript

    if x_start is None:
        x_start = Interval.convert_ratioed_pos(mapping, ust.start).start
    if x_end is None:
        x_end = Interval.convert_ratioed_pos(mapping, ust.end).end

    LABEL_PREFIX = DS.TRANSCRIPT_LABEL_PREFIX
    if isinstance(ust, FusionTranscript):
        LABEL_PREFIX = DS.FUSION_LABEL_PREFIX

    # if the splicing takes up more room than the track we need to adjust for it
    y = DS.SPLICE_HEIGHT

    exon_track_group = _draw_exon_track(DS, canvas, ust, mapping, colors, translation=translation)
    exon_track_group.translate(0, y)
    exon_track_group.add(canvas.text(
        labels.add(tr, LABEL_PREFIX),
        insert=(
            0 - DS.PADDING,
            DS.TRACK_HEIGHT / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE
        ),
        fill=DS.LABEL_COLOR,
        style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='end'),
        class_='label'
    ))

    # draw the splicing pattern
    splice_group = canvas.g(class_='splicing')
    for p1, p2 in zip(tr.splicing_pattern[::2], tr.splicing_pattern[1::2]):
        a = Interval.convert_pos(mapping, p1)
        b = Interval.convert_pos(mapping, p2)
        polyline = [(a, y), (a + (b - a) / 2, y - DS.SPLICE_HEIGHT), (b, y)]
        p = canvas.polyline(polyline, fill='none')
        p.dasharray(DS.SPLICE_STROKE_DASHARRAY)
        p.stroke(DS.SPLICE_COLOR, width=DS.SPLICE_STROKE_WIDTH)
        splice_group.add(p)

    y += DS.TRACK_HEIGHT / 2

    main_group.add(splice_group)
    main_group.add(exon_track_group)
    y += DS.TRACK_HEIGHT / 2

    gp = canvas.g(class_='protein')
    y += DS.PADDING
    gp.translate(0, y)
    # translation track

    # convert the AA position to cdna position, then convert the cdna to genomic, etc
    # ==================== adding the translation track ============
    translated_genomic_regions = [
        tr.convert_cdna_to_genomic(translation.start),
        tr.convert_cdna_to_genomic(translation.end)
    ]
    translated_genomic_regions = [Interval(*sorted(translated_genomic_regions))]

    for p1, p2 in zip(tr.splicing_pattern[::2], tr.splicing_pattern[1::2]):
        try:
            spliced_out_interval = Interval(p1 + 1, p2 - 1)
            temp = []
            for region in translated_genomic_regions:
                temp.extend(region - spliced_out_interval)
            translated_genomic_regions = temp
        except AttributeError:
            pass

    s = Interval.convert_pos(mapping, translated_genomic_regions[0].start)
    t = Interval.convert_pos(mapping, translated_genomic_regions[-1].end)

    gt = canvas.g(class_='translation')
    gp.add(gt)
    h = DS.TRANSLATION_TRACK_HEIGHT

    for sec in translated_genomic_regions:
        start = Interval.convert_pos(mapping, sec.start)
        end = Interval.convert_pos(mapping, sec.end)
        gt.add(
            canvas.rect(
                (start, h / 2 - DS.TRANSLATION_TRACK_HEIGHT / 2), (end - start + 1, DS.TRANSLATION_TRACK_HEIGHT),
                fill=DS.TRANSLATION_SCAFFOLD_COLOR,
                class_='scaffold'
            ))
    gt.add(canvas.text(
        DS.TRANSLATION_END_MARKER if tr.get_strand() == STRAND.NEG else DS.TRANSLATION_START_MARKER,
        insert=(
            s - DS.TRANSLATION_MARKER_PADDING, h / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.TRANSLATION_FONT_SIZE
        ),
        fill=DS.LABEL_COLOR,
        style=DS.FONT_STYLE.format(font_size=DS.TRANSLATION_FONT_SIZE, text_anchor='end'),
        class_='label'
    ))
    gt.add(canvas.text(
        DS.TRANSLATION_START_MARKER if tr.get_strand() == STRAND.NEG else DS.TRANSLATION_END_MARKER,
        insert=(
            t + DS.TRANSLATION_MARKER_PADDING, h / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.TRANSLATION_FONT_SIZE
        ),
        fill=DS.LABEL_COLOR,
        style=DS.FONT_STYLE.format(font_size=DS.TRANSLATION_FONT_SIZE, text_anchor='start'),
        class_='label'
    ))
    gt.add(Tag('title', 'translation  cdna({}_{})  c.{}_{}  p.{}_{}'.format(
        translation.start, translation.end, 1, len(translation), 1, len(translation) // CODON_SIZE)))
    py = h
    # now draw the domain tracks
    # need to convert the domain AA positions to cds positions to genomic
    for i, d in enumerate(sorted(translation.domains, key=lambda x: x.name)):
        if not re.match(DS.DOMAIN_NAME_REGEX_FILTER, str(d.name)):
            continue
        py += DS.PADDING
        gd = canvas.g(class_='domain')
        gd.add(canvas.rect(
            (x_start, DS.DOMAIN_TRACK_HEIGHT / 2), (x_end - x_start, DS.DOMAIN_SCAFFOLD_HEIGHT),
            fill=DS.DOMAIN_SCAFFOLD_COLOR, class_='scaffold'
        ))
        fill = DS.DOMAIN_COLOR
        percent_match = None
        try:
            match, total = d.score_region_mapping(REFERENCE_GENOME)
            percent_match = int(round(match * 100 / total, 0))
            fill = DS.DOMAIN_FILL_GRADIENT[percent_match % len(DS.DOMAIN_FILL_GRADIENT) - 1]
        except (NotSpecifiedError, AttributeError):
            pass
        for region in d.regions:
            # convert the AA position to cdna position, then convert the cdna to genomic, etc
            s = translation.convert_aa_to_cdna(region.start)
            t = translation.convert_aa_to_cdna(region.end)
            s = Interval.convert_pos(mapping, tr.convert_cdna_to_genomic(s.start))
            t = Interval.convert_pos(mapping, tr.convert_cdna_to_genomic(t.end))
            if s > t:
                t, s = (s, t)
            gdr = canvas.g(class_='domain_region')
            gdr.add(canvas.rect((s, 0), (t - s + 1, DS.DOMAIN_TRACK_HEIGHT), fill=fill, class_='region'))
            gdr.add(Tag('title', 'domain {} region p.{}_{}{}'.format(
                d.name if d.name else '', region.start, region.end,
                '  matched({}%)'.format(percent_match) if percent_match is not None else '')))
            gd.add(gdr)
        gd.translate(0, py)

        f = DS.LABEL_COLOR if not DS.DYNAMIC_LABELS else dynamic_label_color(DS.DOMAIN_COLOR)
        label_group = None
        if re.match(DS.PFAM_DOMAIN, d.name):
            label_group = canvas.a(
                DS.PFAM_LINK.format(d), target='_blank')
        else:
            label_group = canvas.g()
        gd.add(label_group)
        label_group.add(canvas.text(
            labels.add(d.name, DS.DOMAIN_LABEL_PREFIX),
            insert=(
                0 - DS.PADDING,
                DS.DOMAIN_TRACK_HEIGHT / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.DOMAIN_LABEL_FONT_SIZE),
            fill=f, class_='label',
            style=DS.FONT_STYLE.format(font_size=DS.DOMAIN_LABEL_FONT_SIZE, text_anchor='end')
        ))
        gp.add(gd)
        py += DS.DOMAIN_TRACK_HEIGHT
    y += py
    main_group.add(gp)
    setattr(main_group, 'height', y)
    return main_group


def draw_ustranscript(
    DS, canvas, ust, target_width=None, breakpoints=[], labels=LabelMapping(), colors={},
    mapping=None, REFERENCE_GENOME=None, masks=None
):
    """
    builds an svg group representing the transcript. Exons are drawn in a track with the splicing
    information and domains are drawn in separate tracks below

    if there are multiple splicing variants then multiple exon tracks are drawn

    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        target_width (int): the target width of the diagram
        t (Transcript): the transcript being drawn
        exon_color (str): the color being used for the fill of the exons
        utr_color (str): the color for the fill of the UTR regions
        abrogated_splice_sites (:class:`list` of :class:`int`): list of positions to ignore as splice sites
        breakpoints (:class:`list` of :class:`Breakpoint`): the breakpoints to overlay

    Return:
        svgwrite.container.Group: the group element for the transcript diagram
                Has the added parameters of labels, height, and mapping
    """


    if ust.get_strand() not in [STRAND.POS, STRAND.NEG]:
        raise NotSpecifiedError('strand must be positive or negative to draw the ust')
    if (mapping is None and target_width is None) or (mapping is not None and target_width is not None):
        raise AttributeError('mapping and target_width arguments are required and mutually exclusive')

    if mapping is None:
        mapping = generate_interval_mapping(
            ust.exons,
            target_width,
            DS.EXON_INTRON_RATIO,
            DS.EXON_MIN_WIDTH,
            min_inter_width=DS.MIN_WIDTH
        )

    main_group = canvas.g(class_='ust')

    y = DS.BREAKPOINT_TOP_MARGIN if len(breakpoints) > 0 else 0
    x_start = Interval.convert_ratioed_pos(mapping, ust.start).start
    x_end = Interval.convert_ratioed_pos(mapping, ust.end).end

    if target_width:
        x_start = 0
        x_end = target_width

    if masks is None:
        masks = []
        try:
            if len(breakpoints) == 1:
                b = breakpoints[0]
                if b.orient == ORIENT.RIGHT:
                    masks = [Interval(ust.start, b.start - 1)]
                elif b.orient == ORIENT.LEFT:
                    masks = [Interval(b.end + 1, ust.end)]
            elif len(breakpoints) == 2:
                b1, b2 = sorted(breakpoints)
                if b1.orient == ORIENT.LEFT and b2.orient == ORIENT.RIGHT:
                    masks = [Interval(b1.end + 1, b2.start - 1)]
        except AttributeError:
            pass

    LABEL_PREFIX = DS.TRANSCRIPT_LABEL_PREFIX
    if isinstance(ust, FusionTranscript):
        LABEL_PREFIX = DS.FUSION_LABEL_PREFIX

    if len(ust.translations) == 0:
        y += DS.SPLICE_HEIGHT
        exon_track_group = _draw_exon_track(DS, canvas, ust, mapping, colors)
        exon_track_group.translate(0, y)
        exon_track_group.add(canvas.text(
            labels.add(ust, LABEL_PREFIX),
            insert=(0 - DS.PADDING, DS.TRACK_HEIGHT / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE),
            fill=DS.LABEL_COLOR,
            style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='end'),
            class_='label'
        ))
        main_group.add(exon_track_group)
        y += DS.TRACK_HEIGHT
    else:
        # draw the protein features if there are any
        for i, tl in enumerate(ust.translations):
            gp = _draw_transcript_with_translation(DS,
                canvas, tl, labels, colors, mapping, x_start=x_start, x_end=x_end
            )
            gp.translate(0, y)
            if i < len(ust.translations) - 1:
                y += DS.INNER_MARGIN
            y += gp.height
            main_group.add(gp)

    y += DS.BREAKPOINT_BOTTOM_MARGIN if len(breakpoints) > 0 else 0
    # add masks
    for mask in masks:
        pixel = Interval.convert_ratioed_pos(mapping, mask.start) | Interval.convert_ratioed_pos(mapping, mask.end)
        m = canvas.rect(
            (pixel.start, 0), (pixel.length(), y),
            class_='mask', fill=DS.MASK_FILL, opacity=DS.MASK_OPACITY, pointer_events='none'
        )
        main_group.add(m)

    # now overlay the breakpoints on top of everything
    for i, b in enumerate(breakpoints):
        pixel = Interval.convert_ratioed_pos(mapping, b.start) | Interval.convert_ratioed_pos(mapping, b.end)
        bg = draw_breakpoint(DS, canvas, b, pixel.length(), y, label=labels.add(b, DS.BREAKPOINT_LABEL_PREFIX))
        bg.translate(pixel.start, 0)
        main_group.add(bg)

    setattr(main_group, 'height', y)
    setattr(main_group, 'width', x_end - x_start)
    setattr(main_group, 'mapping', mapping)
    setattr(main_group, 'labels', labels)
    return main_group


def draw_genes(DS, canvas, genes, target_width, breakpoints=None, colors=None, labels=None, plots=None, masks=None):
    """
    draws the genes given in order of their start position trying to minimize
    the number of tracks required to avoid overlap

    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        target_width (int): the target width of the diagram
        genes (:class:`list` of :class:`Gene`): the list of genes to draw
        breakpoints (:class:`list` of :class:`Breakpoint`): the breakpoints to overlay
        colors (:class:`dict` of :class:`Gene` and :class:`str`): dictionary of the colors assigned to each Gene as
         fill

    Return:
        svgwrite.container.Group: the group element for the diagram.
            Has the added parameters of labels, height, and mapping
    """
    # mutable default argument parameters
    breakpoints = [] if breakpoints is None else breakpoints
    colors = {} if colors is None else colors
    labels = LabelMapping() if labels is None else labels
    plots = plots if plots else []

    st = max(min([g.start for g in genes] + [b.start for b in breakpoints]) - DS.GENE_MIN_BUFFER, 1)
    end = max([g.end for g in genes] + [b.end for b in breakpoints]) + DS.GENE_MIN_BUFFER
    main_group = canvas.g(class_='genes')
    mapping = generate_interval_mapping(
        [g for g in genes],
        target_width,
        DS.GENE_INTERGENIC_RATIO,
        DS.GENE_MIN_WIDTH,
        start=st, end=end,
        min_inter_width=DS.MIN_WIDTH
    )
    if masks is None:
        masks = []
        try:
            if len(breakpoints) == 1:
                b = breakpoints[0]
                if b.orient == ORIENT.RIGHT:
                    masks = [Interval(st, b.start - 1)]
                elif b.orient == ORIENT.LEFT:
                    masks = [Interval(b.end + 1, end)]
            elif len(breakpoints) == 2:
                b1, b2 = sorted(breakpoints)
                if b1.orient == ORIENT.LEFT and b2.orient == ORIENT.RIGHT:
                    masks = [Interval(b1.end, b2.start)]
        except AttributeError:
            pass

    gene_px_intervals = {}
    for i, gene in enumerate(sorted(genes, key=lambda x: x.start)):
        s = Interval.convert_ratioed_pos(mapping, gene.start)
        t = Interval.convert_ratioed_pos(mapping, gene.end)
        gene_px_intervals[Interval(s.start, t.end)] = gene
        l = labels.add(gene, DS.GENE_LABEL_PREFIX)
    tracks = split_intervals_into_tracks(gene_px_intervals)

    y = DS.BREAKPOINT_TOP_MARGIN

    main_group.add(
        canvas.rect(
            (0, y + DS.TRACK_HEIGHT / 2 - DS.SCAFFOLD_HEIGHT / 2 +
                (len(tracks) - 1) * (DS.TRACK_HEIGHT + DS.PADDING)),
            (target_width, DS.SCAFFOLD_HEIGHT),
            fill=DS.SCAFFOLD_COLOR,
            class_='scaffold'
        ))
    tracks.reverse()
    for track in tracks:  # svg works from top down
        for genepx in track:
            # draw the gene
            gene = gene_px_intervals[genepx]
            group = draw_gene(DS,
                canvas, gene, genepx.length(),
                DS.TRACK_HEIGHT,
                colors.get(gene, DS.GENE1_COLOR),
                labels.get_key(gene)
            )
            group.translate(genepx.start, y)
            main_group.add(group)
        y += DS.TRACK_HEIGHT + DS.PADDING

    y += DS.BREAKPOINT_BOTTOM_MARGIN - DS.PADDING

    # adding the masks is the final step
    for mask in masks:
        pixel = Interval.convert_ratioed_pos(mapping, mask.start) | Interval.convert_ratioed_pos(mapping, mask.end)
        m = canvas.rect(
            (pixel.start, 0), (pixel.length(), y),
            class_='mask', fill=DS.MASK_FILL, pointer_events='none', opacity=DS.MASK_OPACITY
        )
        main_group.add(m)

    # now overlay the breakpoints on top of everything
    for i, b in enumerate(sorted(breakpoints)):
        s = Interval.convert_pos(mapping, b.start)
        t = Interval.convert_pos(mapping, b.end)
        bg = draw_breakpoint(DS, canvas, b, abs(t - s) + 1, y, label=labels.add(b, DS.BREAKPOINT_LABEL_PREFIX))
        bg.translate(s, 0)
        main_group.add(bg)

    setattr(main_group, 'height', y)
    setattr(main_group, 'width', target_width)
    setattr(main_group, 'mapping', mapping)
    setattr(main_group, 'labels', labels)

    return main_group


def draw_ustranscripts_overlay(DS, gene, vmarkers=None, window_buffer=0, plots=None):
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
        g = _draw_exon_track(DS, canvas, tx, mapping, colors=colors)
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


def draw_vmarker(DS, canvas, marker, width, height, label='', color=None):
    """
    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        breakpoint (Breakpoint): the breakpoint to draw
        width (int): the pixel width
        height (int): the pixel height
    Return:
        svgwrite.container.Group: the group element for the diagram
    """
    color = DS.MARKER_COLOR if color is None else color
    width = max([DS.ABS_MIN_WIDTH, width])
    g = canvas.g(class_='marker')
    y = DS.PADDING + DS.MARKER_LABEL_FONT_SIZE / 2
    t = canvas.text(
        label,
        insert=(width / 2, y + DS.FONT_CENTRAL_SHIFT_RATIO * DS.MARKER_LABEL_FONT_SIZE),
        fill=color,
        style=DS.FONT_STYLE.format(text_anchor='middle', font_size=DS.MARKER_LABEL_FONT_SIZE),
        class_='label'
    )
    y += DS.MARKER_LABEL_FONT_SIZE / 2 + DS.PADDING

    g.add(t)
    g.add(canvas.rect((0, y), (width, height - y), stroke=color, fill='none'))
    g.add(Tag('title', 'marker {}:{}-{} {}'.format(
        marker.reference_object, marker.start, marker.end, marker.name)))
    return g


def draw_breakpoint(DS, canvas, breakpoint, width, height, label=''):
    """
    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        breakpoint (Breakpoint): the breakpoint to draw
        width (int): the pixel width
        height (int): the pixel height
    Return:
        svgwrite.container.Group: the group element for the diagram
    """
    g = canvas.g(class_='breakpoint')
    y = DS.PADDING + DS.BREAKPOINT_LABEL_FONT_SIZE / 2
    t = canvas.text(
        label,
        insert=(width / 2, y + DS.FONT_CENTRAL_SHIFT_RATIO * DS.BREAKPOINT_LABEL_FONT_SIZE),
        fill=HEX_BLACK,
        style=DS.FONT_STYLE.format(text_anchor='middle', font_size=DS.BREAKPOINT_LABEL_FONT_SIZE),
        class_='label'
    )
    y += DS.BREAKPOINT_LABEL_FONT_SIZE / 2 + DS.PADDING

    g.add(t)

    r = canvas.rect((0, y), (width, height - y), stroke=DS.BREAKPOINT_COLOR, fill='none')
    r.dasharray(DS.BREAKPOINT_STROKE_DASHARRAY)
    g.add(r)

    if breakpoint.orient == ORIENT.LEFT:
        l = canvas.line((0, y), (0, height))
        l.stroke(DS.BREAKPOINT_COLOR, width=DS.BREAKPOINT_ORIENT_STROKE_WIDTH)
        g.add(l)
    elif breakpoint.orient == ORIENT.RIGHT:
        l = canvas.line((width, y), (width, height))
        l.stroke(DS.BREAKPOINT_COLOR, width=DS.BREAKPOINT_ORIENT_STROKE_WIDTH)
        g.add(l)
    g.add(Tag('title', 'Breakpoint {}:g.{}_{}{} {}'.format(
        breakpoint.chr, breakpoint.start, breakpoint.end, breakpoint.strand, breakpoint.orient)))
    return g


def draw_exon(DS, canvas, exon, width, height, fill, label='', translation=None):
    """
    generates the svg object representing an exon

    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        exon (Exon): the exon to draw
        width (int): the pixel width
        height (int): the pixel height
        fill (str): the fill color to use for the exon

    Return:
        svgwrite.container.Group: the group element for the diagram

    .. todo::
        add markers for exons with abrogated splice sites
    """
    g = canvas.g(class_='exon')
    label = str(label)
    g.add(canvas.rect((0, 0), (width, height), fill=fill))
    t = canvas.text(
        label,
        insert=(width / 2, height / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.EXON_FONT_SIZE),
        fill=DS.LABEL_COLOR if not DS.DYNAMIC_LABELS else dynamic_label_color(fill),
        style=DS.FONT_STYLE.format(font_size=DS.EXON_FONT_SIZE, text_anchor='middle'),
        class_='label'
    )
    g.add(t)
    title = 'Exon {}  g.{}_{}'.format(
        exon.name if exon.name else '', exon.start, exon.end)
    if translation:
        cds_start = translation.convert_genomic_to_cds_notation(exon.start)
        cds_end = translation.convert_genomic_to_cds_notation(exon.end)
        if exon.get_strand() == STRAND.NEG:
            cds_start, cds_end = cds_end, cds_start
        title += '  c.{}_{}'.format(cds_start, cds_end)
        cdna_start = translation.transcript.convert_genomic_to_cdna(exon.start)
        cdna_end = translation.transcript.convert_genomic_to_cdna(exon.end)
        if cdna_end < cdna_start:
            cdna_start, cdna_end = cdna_end, cdna_start
        title += '  cdna({}_{})'.format(cdna_start, cdna_end)
    title += '  length({})'.format(len(exon))
    g.add(Tag('title', title))
    return g


def draw_template(DS, canvas, template, target_width, labels=None, colors=None, breakpoints=None):
    labels = LabelMapping() if labels is None else labels
    colors = {} if colors is None else colors
    breakpoints = [] if not breakpoints else breakpoints
    total_height = DS.TEMPLATE_TRACK_HEIGHT + DS.BREAKPOINT_TOP_MARGIN + DS.BREAKPOINT_BOTTOM_MARGIN
    group = canvas.g(class_='template')
    # 1 as input since we don't want to change the ratio here
    mapping = generate_interval_mapping(
        template.bands, target_width, 1, DS.TEMPLATE_BAND_MIN_WIDTH,
        start=template.start, end=template.end
    )
    scaffold = canvas.rect(
        (0, 0), (target_width, DS.SCAFFOLD_HEIGHT),
        fill=DS.SCAFFOLD_COLOR
    )
    group.add(scaffold)
    scaffold.translate((0, DS.BREAKPOINT_TOP_MARGIN + DS.TEMPLATE_TRACK_HEIGHT / 2 - DS.SCAFFOLD_HEIGHT / 2))
    label_group = canvas.g()
    label_group.add(canvas.text(
        labels.add(template, DS.TEMPLATE_LABEL_PREFIX),
        insert=(
            0 - DS.PADDING, DS.BREAKPOINT_TOP_MARGIN + DS.TEMPLATE_TRACK_HEIGHT / 2 +
            DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE),
        fill=DS.LABEL_COLOR,
        style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='end'),
        class_='label'
    ))
    label_group.add(Tag('title', 'template {}'.format(template.name)))
    group.add(label_group)

    for band in template.bands:
        s = Interval.convert_pos(mapping, band[0])
        t = Interval.convert_pos(mapping, band[1])

        bgroup = canvas.g(class_='cytoband')
        f = DS.TEMPLATE_BAND_FILL.get(band.data.get('giesma_stain', None), DS.TEMPLATE_DEFAULT_FILL)
        r = None
        w = t - s + 1
        if band.data.get('giesma_stain', None) == GIESMA_STAIN.ACEN:
            if band.name[0] == 'p':
                r = canvas.polyline(
                    [(0, 0), (w, DS.TEMPLATE_TRACK_HEIGHT / 2), (0, DS.TEMPLATE_TRACK_HEIGHT)],
                    fill=f, stroke=DS.TEMPLATE_BAND_STROKE, stroke_width=DS.TEMPLATE_BAND_STROKE_WIDTH
                )
            else:
                r = canvas.polyline(
                    [(w, 0), (0, DS.TEMPLATE_TRACK_HEIGHT / 2), (w, DS.TEMPLATE_TRACK_HEIGHT)],
                    fill=f, stroke=DS.TEMPLATE_BAND_STROKE, stroke_width=DS.TEMPLATE_BAND_STROKE_WIDTH
                )
        else:
            r = canvas.rect(
                (0, 0), (w, DS.TEMPLATE_TRACK_HEIGHT),
                fill=f, stroke=DS.TEMPLATE_BAND_STROKE, stroke_width=DS.TEMPLATE_BAND_STROKE_WIDTH
            )
        bgroup.add(r)
        bgroup.add(
            Tag('title', 'cytoband  {0}:y.{1}  {0}:g.{2}_{3}'.format(
                template.name, band.name, band.start, band.end)))
        bgroup.translate((s, DS.BREAKPOINT_TOP_MARGIN))
        group.add(bgroup)
    # now draw the breakpoints overtop
    for i, b in enumerate(sorted(breakpoints)):
        s = Interval.convert_pos(mapping, b.start)
        t = Interval.convert_pos(mapping, b.end)
        bg = draw_breakpoint(DS,
            canvas, b, abs(t - s) + 1, total_height, label=labels.add(b, DS.BREAKPOINT_LABEL_PREFIX))
        bg.translate(s, 0)
        group.add(bg)
    setattr(
        group, 'height', total_height)
    return group


def draw_gene(DS, canvas, gene, width, height, fill, label='', REFERENCE_GENOME=None):
    """
    generates the svg object representing a gene

    Args:
        canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
        gene (Gene): the gene to draw
        width (int): the pixel width
        height (int): the pixel height
        fill (str): the fill color to use for the gene

    Return:
        svgwrite.container.Group: the group element for the diagram
    """

    group = canvas.g(class_='gene')
    if width < DS.GENE_MIN_WIDTH:
        raise DrawingFitError('width of {} is not sufficient to draw a gene of minimum width {}'.format(
            width, DS.GENE_MIN_WIDTH))
    wrect = width - DS.GENE_ARROW_WIDTH
    if wrect < 1:
        raise DrawingFitError('width is not sufficient to draw gene')

    label_color = DS.LABEL_COLOR if not DS.DYNAMIC_LABELS else dynamic_label_color(fill)

    if gene.get_strand() == STRAND.POS:
        group.add(
            canvas.rect(
                (0, 0), (wrect, height), fill=fill
            ))
        group.add(
            canvas.polyline(
                [(wrect, 0), (wrect + DS.GENE_ARROW_WIDTH, height / 2), (wrect, height)],
                fill=fill
            ))
        group.add(
            canvas.text(
                label,
                insert=(wrect / 2, height / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE),
                fill=label_color,
                style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='middle'),
                class_='label'
            ))
    elif gene.get_strand() == STRAND.NEG:
        group.add(
            canvas.rect(
                (DS.GENE_ARROW_WIDTH, 0), (wrect, height), fill=fill
            ))
        group.add(
            canvas.polyline(
                [(DS.GENE_ARROW_WIDTH, 0), (0, height / 2), (DS.GENE_ARROW_WIDTH, height)],
                fill=fill
            ))
        group.add(
            canvas.text(
                label,
                insert=(
                    wrect / 2 + DS.GENE_ARROW_WIDTH,
                    height / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE
                ),
                fill=label_color,
                style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='middle'),
                class_='label'
            ))
    else:
        group.add(
            canvas.rect(
                (0, 0), (width, height), fill=fill
            ))
        group.add(
            canvas.text(
                label,
                insert=(width / 2, height / 2 + DS.FONT_CENTRAL_SHIFT_RATIO * DS.LABEL_FONT_SIZE),
                fill=label_color,
                style=DS.FONT_STYLE.format(font_size=DS.LABEL_FONT_SIZE, text_anchor='middle'),
                class_='label'
            ))
    aliases = ''
    try:
        if len(gene.aliases) > 0:
            aliases = ' aka {}'.format(';'.join(sorted(gene.aliases)))
    except AttributeError:
        pass
    group.add(
        Tag('title', 'Gene {} {}:g.{}_{}{}{}'.format(gene.name if gene.name else '',
            gene.chr, gene.start, gene.end, gene.get_strand(), aliases)))
    return group


def split_intervals_into_tracks(intervals):
    tracks = [[]]
    for i in sorted(intervals, key=lambda x: x[0]):
        added = False
        for t in tracks:
            overlaps = False
            for og in t:
                if Interval.overlaps(i, og):
                    overlaps = True
                    break
            if not overlaps:
                added = True
                t.append(i)
                break
        if not added:
            tracks.append([i])
    return tracks


def generate_interval_mapping(
    input_intervals, target_width, ratio, min_width,
    buffer_length=None, start=None, end=None, min_inter_width=None
):
    min_inter_width = min_width if min_inter_width is None else min_inter_width
    if all([x is not None for x in [start, end, buffer_length]]):
        raise AttributeError('buffer_length is a mutually exclusive argument with start/end')

    intervals = []
    for i in Interval.min_nonoverlapping(*input_intervals):
        if len(intervals) == 0 or abs(Interval.dist(intervals[-1], i)) > 1:
            intervals.append(i)
        else:
            intervals[-1] = intervals[-1] | i
    # break up the intervals by any intervals of length 1
    for itvl_in in input_intervals:
        if len(itvl_in) > 1:
            continue
        # try splitting all current interval
        temp = []
        for itvl in intervals:
            split = itvl - itvl_in
            if split is not None:
                temp.extend(split)
        intervals = temp
    for itvl_in in input_intervals:
        if len(itvl_in) == 1:
            intervals.append(itvl_in)
    # now split any intervals by start/end
    breaks = {}
    for i in intervals:
        # split by input intervals
        breaks[i] = set([i.start, i.end])
        for ii in input_intervals:
            if ii.start >= i.start and ii.start <= i.end:
                breaks[i].add(ii.start)
            if ii.end >= i.start and ii.end <= i.end:
                breaks[i].add(ii.end)
    temp = []
    for itvl, breakpoints in breaks.items():
        breakpoints.add(itvl.start)
        breakpoints.add(itvl.end)
        pos = sorted(breakpoints)
        if len(pos) == 1:
            temp.append(Interval(pos[0]))
        else:
            # remove all the single intervals to start?
            pos[0] -= 1
            for i in range(1, len(pos)):
                temp.append(Interval(pos[i - 1] + 1, pos[i]))
    intervals = sorted(temp, key=lambda x: x.start)

    if buffer_length is None:
        buffer_length = 0

    if start is None:
        start = max(intervals[0].start - buffer_length, 1)
    elif start <= 0:
        raise AttributeError('start must be a natural number', start)

    if end is None:
        end = intervals[-1].end + buffer_length
    elif end <= 0:
        raise AttributeError('end must be a natural number', end)

    total_length = end - start + 1
    genic_length = sum([len(i) for i in intervals])
    intergenic_length = total_length - genic_length
    intermediate_intervals = 0
    if start < intervals[0].start:
        intermediate_intervals += 1
    if end > intervals[-1].end:
        intermediate_intervals += 1

    for i in range(1, len(intervals)):
        if intervals[i].start > intervals[i - 1].end + 1:
            intermediate_intervals += 1
    width = target_width - intermediate_intervals * min_inter_width - len(intervals) * min_width  # reserved width

    if width < 0:
        raise DrawingFitError('width cannot accommodate the number of expected objects')

    intergenic_width = width // (ratio + 1) if intergenic_length > 0 else 0
    genic_width = width - intergenic_width
    intergenic_unit = lambda x: x * intergenic_width / intergenic_length
    genic_unit = lambda x: x * genic_width / genic_length

    assert(
        genic_width + intergenic_width + len(intervals) * min_width + intermediate_intervals * min_inter_width ==
        target_width)
    mapping = []

    pos = 1
    # do the intergenic region prior to the first genic region
    if start < intervals[0].start:
        ifrom = Interval(start, intervals[0].start - 1)
        s = max(intergenic_unit(len(ifrom)), 0)
        ito = Interval(pos, pos + min_inter_width + s)
        mapping.append((ifrom, ito))
        pos += ito.length()

    for i, curr in enumerate(intervals):
        if i > 0 and intervals[i - 1].end + 1 < curr.start:  # add between the intervals
            prev = intervals[i - 1]
            ifrom = Interval(prev.end + 1, curr.start - 1)
            s = max(intergenic_unit(len(ifrom)), 0)
            ito = Interval(pos, pos + min_inter_width + s)
            mapping.append((ifrom, ito))
            pos += ito.length()

        s = max(genic_unit(len(curr)), 0)
        ito = Interval(pos, pos + min_width + s)
        mapping.append((curr, ito))
        pos += ito.length()

    # now the last intergenic region will make up for the rounding error
    if end > intervals[-1].end:
        ifrom = Interval(intervals[-1].end + 1, end)
        s = max(intergenic_unit(len(ifrom)), 0)
        ito = Interval(pos, pos + min_inter_width + s)
        mapping.append((ifrom, ito))
        pos += ito.length()
    mapping[-1][1].end = target_width  # min(int(target_width), mapping[-1][1].end)
    temp = mapping
    mapping = dict()
    for ifrom, ito in temp:
        mapping[ifrom] = ito

    # assert that that mapping is correct
    for ifrom in input_intervals:
        ifrom = Interval(ifrom.start, ifrom.end)
        p1 = Interval.convert_ratioed_pos(mapping, ifrom.start)
        p2 = Interval.convert_ratioed_pos(mapping, ifrom.end)
        if ifrom in mapping and ito.end == target_width:
            continue
        n = p1 | p2
        if n.length() < min_width and abs(n.length() - min_width) > 0.01:  # precision error allowable
            raise AssertionError(
                'interval mapping should not map any intervals to less than the minimum required width. Interval {}'
                ' was mapped to a pixel interval of length {} but the minimum width is {}'.format(
                    ifrom, n.length(), min_width),
                p1, p2, mapping, input_intervals, target_width, ratio, min_inter_width)
    return mapping
