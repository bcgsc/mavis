"""
This is the primary module responsible for generating svg visualizations

"""
from ..annotate.genomic import IntergenicRegion
from ..annotate.variant import FusionTranscript
from ..constants import STRAND, ORIENT, CODON_SIZE, GIESMA_STAIN
from ..error import DrawingFitError, NotSpecifiedError
from ..interval import Interval
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

def draw_legend(DS, canvas, swatches, border=True):
    main_group = canvas.g(class_='legend')
    y = DS.padding if border else 0
    x = DS.padding if border else 0
    for swatch, label in swatches:
        g = canvas.g()
        g.add(canvas.rect(
            (0, 0),
            (DS.legend_swatch_size, DS.legend_swatch_size),
            fill=swatch,
            stroke=DS.legend_swatch_stroke
        ))

        g.add(canvas.text(
            label,
            insert=(DS.legend_swatch_size + DS.padding, DS.legend_swatch_size / 2),
            fill=DS.legend_font_color,
            style=DS.font_style.format(text_anchor='start', font_size=DS.legend_font_size),
            class_='label'
        ))
        g.translate(x, y)
        main_group.add(g)
        y += DS.legend_swatch_size + DS.padding

    w = max([len(l) for c, l in swatches]) * DS.legend_font_size * DS.font_width_height_ratio + \
        DS.padding * (3 if border else 1) + DS.legend_swatch_size

    if border:
        main_group.add(canvas.rect(
            (0, 0), (w, y), fill='none', stroke=DS.legend_border_stroke,
            stroke_width=DS.legend_border_stroke_width
        ))
    else:
        y -= DS.padding
    setattr(main_group, 'height', y)
    setattr(main_group, 'width', w)
    setattr(main_group, 'labels', None)
    setattr(main_group, 'mapping', None)
    return main_group


def draw_exon_track(DS, canvas, transcript, mapping, colors=None, x_start=None, x_end=None, translation=None):
    """
    """
    colors = {} if colors is None else colors
    main_group = canvas.g(class_='exon_track')

    y = DS.track_height / 2
    exons = sorted(transcript.exons, key=lambda x: x.start)

    s = Interval.convert_ratioed_pos(mapping, exons[0].start).start if x_start is None else x_start
    t = Interval.convert_ratioed_pos(mapping, exons[-1].end).end if x_end is None else x_end

    main_group.add(
        canvas.rect(
            (s, y - DS.scaffold_height / 2), (t - s + 1, DS.scaffold_height),
            fill=DS.scaffold_color, class_='scaffold'
        ))

    # draw the exons
    for exon in exons:
        s = Interval.convert_ratioed_pos(mapping, exon.start).start
        t = Interval.convert_ratioed_pos(mapping, exon.end).end
        pxi = Interval(s, t)
        c = colors.get(exon, DS.exon1_color)
        group = draw_exon(DS,
            canvas, exon, pxi.length(), DS.track_height, c,
            label=transcript.exon_number(exon),
            translation=translation
        )
        group.translate(pxi.start, y - DS.track_height / 2)
        main_group.add(group)

    setattr(main_group, 'height', y + DS.track_height / 2)
    setattr(main_group, 'width', t - s + 1)
    return main_group


def draw_transcript_with_translation(
    DS, canvas, translation, labels, colors, mapping, reference_genome=None, x_start=None, x_end=None
):
    main_group = canvas.g()
    ust = translation.transcript.reference_object
    tr = translation.transcript

    if x_start is None:
        x_start = Interval.convert_ratioed_pos(mapping, ust.start).start
    if x_end is None:
        x_end = Interval.convert_ratioed_pos(mapping, ust.end).end

    LABEL_PREFIX = DS.transcript_label_prefix
    if isinstance(ust, FusionTranscript):
        LABEL_PREFIX = DS.fusion_label_prefix

    # if the splicing takes up more room than the track we need to adjust for it
    y = DS.splice_height

    exon_track_group = draw_exon_track(DS, canvas, ust, mapping, colors, translation=translation)
    exon_track_group.translate(0, y)
    exon_track_group.add(canvas.text(
        labels.add(tr, LABEL_PREFIX),
        insert=(
            0 - DS.padding,
            DS.track_height / 2 + DS.font_central_shift_ratio * DS.label_font_size
        ),
        fill=DS.label_color,
        style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='end'),
        class_='label'
    ))

    # draw the splicing pattern
    splice_group = canvas.g(class_='splicing')
    for p1, p2 in zip(tr.splicing_pattern[::2], tr.splicing_pattern[1::2]):
        a = Interval.convert_pos(mapping, p1)
        b = Interval.convert_pos(mapping, p2)
        polyline = [(a, y), (a + (b - a) / 2, y - DS.splice_height), (b, y)]
        p = canvas.polyline(polyline, fill='none')
        p.dasharray(DS.splice_stroke_dasharray)
        p.stroke(DS.splice_color, width=DS.splice_stroke_width)
        splice_group.add(p)

    y += DS.track_height / 2

    main_group.add(splice_group)
    main_group.add(exon_track_group)
    y += DS.track_height / 2

    gp = canvas.g(class_='protein')
    y += DS.padding
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
    h = DS.translation_track_height

    for sec in translated_genomic_regions:
        start = Interval.convert_pos(mapping, sec.start)
        end = Interval.convert_pos(mapping, sec.end)
        gt.add(
            canvas.rect(
                (start, h / 2 - DS.translation_track_height / 2), (end - start + 1, DS.translation_track_height),
                fill=DS.translation_scaffold_color,
                class_='scaffold'
            ))
    gt.add(canvas.text(
        DS.translation_end_marker if tr.get_strand() == STRAND.NEG else DS.translation_start_marker,
        insert=(
            s - DS.translation_marker_padding, h / 2 + DS.font_central_shift_ratio * DS.translation_font_size
        ),
        fill=DS.label_color,
        style=DS.font_style.format(font_size=DS.translation_font_size, text_anchor='end'),
        class_='label'
    ))
    gt.add(canvas.text(
        DS.translation_start_marker if tr.get_strand() == STRAND.NEG else DS.translation_end_marker,
        insert=(
            t + DS.translation_marker_padding, h / 2 + DS.font_central_shift_ratio * DS.translation_font_size
        ),
        fill=DS.label_color,
        style=DS.font_style.format(font_size=DS.translation_font_size, text_anchor='start'),
        class_='label'
    ))
    gt.add(Tag('title', 'translation  cdna({}_{})  c.{}_{}  p.{}_{}'.format(
        translation.start, translation.end, 1, len(translation), 1, len(translation) // CODON_SIZE)))
    py = h
    # now draw the domain tracks
    # need to convert the domain AA positions to cds positions to genomic
    for i, d in enumerate(sorted(translation.domains, key=lambda x: x.name)):
        if not re.match(DS.domain_name_regex_filter, str(d.name)):
            continue
        py += DS.padding
        gd = canvas.g(class_='domain')
        gd.add(canvas.rect(
            (x_start, DS.domain_track_height / 2), (x_end - x_start, DS.domain_scaffold_height),
            fill=DS.domain_scaffold_color, class_='scaffold'
        ))
        fill = DS.domain_color
        percent_match = None
        try:
            match, total = d.score_region_mapping(reference_genome)
            percent_match = int(round(match * 100 / total, 0))
            fill = DS.domain_fill_gradient[percent_match % len(DS.domain_fill_gradient) - 1]
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
            gdr.add(canvas.rect((s, 0), (t - s + 1, DS.domain_track_height), fill=fill, class_='region'))
            gdr.add(Tag('title', 'domain {} region p.{}_{}{}'.format(
                d.name if d.name else '', region.start, region.end,
                '  matched({}%)'.format(percent_match) if percent_match is not None else '')))
            gd.add(gdr)
        gd.translate(0, py)

        f = DS.label_color if not DS.dynamic_labels else dynamic_label_color(DS.domain_color)
        label_group = None
        if re.match(DS.pfam_domain, d.name):
            label_group = canvas.a(
                DS.pfam_link.format(d), target='_blank')
        else:
            label_group = canvas.g()
        gd.add(label_group)
        label_group.add(canvas.text(
            labels.add(d.name, DS.domain_label_prefix),
            insert=(
                0 - DS.padding,
                DS.domain_track_height / 2 + DS.font_central_shift_ratio * DS.domain_label_font_size),
            fill=f, class_='label',
            style=DS.font_style.format(font_size=DS.domain_label_font_size, text_anchor='end')
        ))
        gp.add(gd)
        py += DS.domain_track_height
    y += py
    main_group.add(gp)
    setattr(main_group, 'height', y)
    return main_group


def draw_ustranscript(
    DS, canvas, ust, target_width=None, breakpoints=[], labels=LabelMapping(), colors={},
    mapping=None, reference_genome=None, masks=None
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
            DS.exon_intron_ratio,
            DS.exon_min_width,
            min_inter_width=DS.min_width
        )

    main_group = canvas.g(class_='ust')

    y = DS.breakpoint_top_margin if len(breakpoints) > 0 else 0
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

    LABEL_PREFIX = DS.transcript_label_prefix
    if isinstance(ust, FusionTranscript):
        LABEL_PREFIX = DS.fusion_label_prefix

    if len(ust.translations) == 0:
        y += DS.splice_height
        exon_track_group = draw_exon_track(DS, canvas, ust, mapping, colors)
        exon_track_group.translate(0, y)
        exon_track_group.add(canvas.text(
            labels.add(ust, LABEL_PREFIX),
            insert=(0 - DS.padding, DS.track_height / 2 + DS.font_central_shift_ratio * DS.label_font_size),
            fill=DS.label_color,
            style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='end'),
            class_='label'
        ))
        main_group.add(exon_track_group)
        y += DS.track_height
    else:
        # draw the protein features if there are any
        for i, tl in enumerate(ust.translations):
            gp = draw_transcript_with_translation(DS,
                canvas, tl, labels, colors, mapping, x_start=x_start, x_end=x_end
            )
            gp.translate(0, y)
            if i < len(ust.translations) - 1:
                y += DS.inner_margin
            y += gp.height
            main_group.add(gp)

    y += DS.breakpoint_bottom_margin if len(breakpoints) > 0 else 0
    # add masks
    for mask in masks:
        pixel = Interval.convert_ratioed_pos(mapping, mask.start) | Interval.convert_ratioed_pos(mapping, mask.end)
        m = canvas.rect(
            (pixel.start, 0), (pixel.length(), y),
            class_='mask', fill=DS.mask_fill, opacity=DS.mask_opacity, pointer_events='none'
        )
        main_group.add(m)

    # now overlay the breakpoints on top of everything
    for i, b in enumerate(breakpoints):
        pixel = Interval.convert_ratioed_pos(mapping, b.start) | Interval.convert_ratioed_pos(mapping, b.end)
        bg = draw_breakpoint(DS, canvas, b, pixel.length(), y, label=labels.add(b, DS.breakpoint_label_prefix))
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

    st = max(min([g.start for g in genes] + [b.start for b in breakpoints]) - DS.gene_min_buffer, 1)
    end = max([g.end for g in genes] + [b.end for b in breakpoints]) + DS.gene_min_buffer
    main_group = canvas.g(class_='genes')
    mapping = generate_interval_mapping(
        [g for g in genes],
        target_width,
        DS.gene_intergenic_ratio,
        DS.gene_min_width,
        start=st, end=end,
        min_inter_width=DS.min_width
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
        l = labels.add(gene, DS.gene_label_prefix)
    tracks = split_intervals_into_tracks(gene_px_intervals)

    y = DS.breakpoint_top_margin

    main_group.add(
        canvas.rect(
            (0, y + DS.track_height / 2 - DS.scaffold_height / 2 +
                (len(tracks) - 1) * (DS.track_height + DS.padding)),
            (target_width, DS.scaffold_height),
            fill=DS.scaffold_color,
            class_='scaffold'
        ))
    tracks.reverse()
    for track in tracks:  # svg works from top down
        for genepx in track:
            # draw the gene
            gene = gene_px_intervals[genepx]
            group = draw_gene(DS,
                canvas, gene, genepx.length(),
                DS.track_height,
                colors.get(gene, DS.gene1_color),
                labels.get_key(gene)
            )
            group.translate(genepx.start, y)
            main_group.add(group)
        y += DS.track_height + DS.padding

    y += DS.breakpoint_bottom_margin - DS.padding

    # adding the masks is the final step
    for mask in masks:
        pixel = Interval.convert_ratioed_pos(mapping, mask.start) | Interval.convert_ratioed_pos(mapping, mask.end)
        m = canvas.rect(
            (pixel.start, 0), (pixel.length(), y),
            class_='mask', fill=DS.mask_fill, pointer_events='none', opacity=DS.mask_opacity
        )
        main_group.add(m)

    # now overlay the breakpoints on top of everything
    for i, b in enumerate(sorted(breakpoints)):
        s = Interval.convert_ratioed_pos(mapping, b.start).start
        t = Interval.convert_ratioed_pos(mapping, b.end).end
        bg = draw_breakpoint(DS, canvas, b, abs(t - s) + 1, y, label=labels.add(b, DS.breakpoint_label_prefix))
        bg.translate(s, 0)
        main_group.add(bg)

    setattr(main_group, 'height', y)
    setattr(main_group, 'width', target_width)
    setattr(main_group, 'mapping', mapping)
    setattr(main_group, 'labels', labels)

    return main_group


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
    color = DS.marker_color if color is None else color
    width = max([DS.abs_min_width, width])
    g = canvas.g(class_='marker')
    y = DS.padding + DS.marker_label_font_size / 2
    t = canvas.text(
        label,
        insert=(width / 2, y + DS.font_central_shift_ratio * DS.marker_label_font_size),
        fill=color,
        style=DS.font_style.format(text_anchor='middle', font_size=DS.marker_label_font_size),
        class_='label'
    )
    y += DS.marker_label_font_size / 2 + DS.padding

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
    y = DS.padding + DS.breakpoint_label_font_size / 2
    t = canvas.text(
        label,
        insert=(width / 2, y + DS.font_central_shift_ratio * DS.breakpoint_label_font_size),
        fill=HEX_BLACK,
        style=DS.font_style.format(text_anchor='middle', font_size=DS.breakpoint_label_font_size),
        class_='label'
    )
    y += DS.breakpoint_label_font_size / 2 + DS.padding

    g.add(t)

    r = canvas.rect((0, y), (width, height - y), stroke=DS.breakpoint_color, fill='none')
    r.dasharray(DS.breakpoint_stroke_dasharray)
    g.add(r)

    if breakpoint.orient == ORIENT.LEFT:
        l = canvas.line((0, y), (0, height))
        l.stroke(DS.breakpoint_color, width=DS.breakpoint_orient_stroke_width)
        g.add(l)
    elif breakpoint.orient == ORIENT.RIGHT:
        l = canvas.line((width, y), (width, height))
        l.stroke(DS.breakpoint_color, width=DS.breakpoint_orient_stroke_width)
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
        insert=(width / 2, height / 2 + DS.font_central_shift_ratio * DS.exon_font_size),
        fill=DS.label_color if not DS.dynamic_labels else dynamic_label_color(fill),
        style=DS.font_style.format(font_size=DS.exon_font_size, text_anchor='middle'),
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
        try:
            cdna_start = translation.transcript.convert_genomic_to_cdna(exon.start)
            cdna_end = translation.transcript.convert_genomic_to_cdna(exon.end)
            if cdna_end < cdna_start:
                cdna_start, cdna_end = cdna_end, cdna_start
            title += '  cdna({}_{})'.format(cdna_start, cdna_end)
        except IndexError:
            title += '  cdna(N/A)'
    title += '  length({})'.format(len(exon))
    g.add(Tag('title', title))
    return g


def draw_template(DS, canvas, template, target_width, labels=None, colors=None, breakpoints=None):
    labels = LabelMapping() if labels is None else labels
    colors = {} if colors is None else colors
    breakpoints = [] if not breakpoints else breakpoints
    total_height = DS.template_track_height + DS.breakpoint_top_margin + DS.breakpoint_bottom_margin
    group = canvas.g(class_='template')
    # 1 as input since we don't want to change the ratio here
    mapping = generate_interval_mapping(
        template.bands, target_width, 1, DS.template_band_min_width,
        start=template.start, end=template.end
    )
    scaffold = canvas.rect(
        (0, 0), (target_width, DS.scaffold_height),
        fill=DS.scaffold_color
    )
    group.add(scaffold)
    scaffold.translate((0, DS.breakpoint_top_margin + DS.template_track_height / 2 - DS.scaffold_height / 2))
    label_group = canvas.g()
    label_group.add(canvas.text(
        labels.add(template, DS.template_label_prefix),
        insert=(
            0 - DS.padding, DS.breakpoint_top_margin + DS.template_track_height / 2 +
            DS.font_central_shift_ratio * DS.label_font_size),
        fill=DS.label_color,
        style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='end'),
        class_='label'
    ))
    label_group.add(Tag('title', 'template {}'.format(template.name)))
    group.add(label_group)

    for band in template.bands:
        s = Interval.convert_pos(mapping, band[0])
        t = Interval.convert_pos(mapping, band[1])

        bgroup = canvas.g(class_='cytoband')
        f = DS.template_band_fill.get(band.data.get('giesma_stain', None), DS.template_default_fill)
        r = None
        w = t - s + 1
        if band.data.get('giesma_stain', None) == GIESMA_STAIN.ACEN:
            if band.name[0] == 'p':
                r = canvas.polyline(
                    [(0, 0), (w, DS.template_track_height / 2), (0, DS.template_track_height)],
                    fill=f, stroke=DS.template_band_stroke, stroke_width=DS.template_band_stroke_width
                )
            else:
                r = canvas.polyline(
                    [(w, 0), (0, DS.template_track_height / 2), (w, DS.template_track_height)],
                    fill=f, stroke=DS.template_band_stroke, stroke_width=DS.template_band_stroke_width
                )
        else:
            r = canvas.rect(
                (0, 0), (w, DS.template_track_height),
                fill=f, stroke=DS.template_band_stroke, stroke_width=DS.template_band_stroke_width
            )
        bgroup.add(r)
        bgroup.add(
            Tag('title', 'cytoband  {0}:y.{1}  {0}:g.{2}_{3}'.format(
                template.name, band.name, band.start, band.end)))
        bgroup.translate((s, DS.breakpoint_top_margin))
        group.add(bgroup)
    # now draw the breakpoints overtop
    for i, b in enumerate(sorted(breakpoints)):
        s = Interval.convert_pos(mapping, b.start)
        t = Interval.convert_pos(mapping, b.end)
        bg = draw_breakpoint(DS,
            canvas, b, abs(t - s) + 1, total_height, label=labels.add(b, DS.breakpoint_label_prefix))
        bg.translate(s, 0)
        group.add(bg)
    setattr(
        group, 'height', total_height)
    return group


def draw_gene(DS, canvas, gene, width, height, fill, label='', reference_genome=None):
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
    if width < DS.gene_min_width:
        raise DrawingFitError('width of {} is not sufficient to draw a gene of minimum width {}'.format(
            width, DS.gene_min_width))
    wrect = width - DS.gene_arrow_width
    if wrect < 1:
        raise DrawingFitError('width is not sufficient to draw gene')

    label_color = DS.label_color if not DS.dynamic_labels else dynamic_label_color(fill)

    if gene.get_strand() == STRAND.POS:
        group.add(
            canvas.rect(
                (0, 0), (wrect, height), fill=fill
            ))
        group.add(
            canvas.polyline(
                [(wrect, 0), (wrect + DS.gene_arrow_width, height / 2), (wrect, height)],
                fill=fill
            ))
        group.add(
            canvas.text(
                label,
                insert=(wrect / 2, height / 2 + DS.font_central_shift_ratio * DS.label_font_size),
                fill=label_color,
                style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='middle'),
                class_='label'
            ))
    elif gene.get_strand() == STRAND.NEG:
        group.add(
            canvas.rect(
                (DS.gene_arrow_width, 0), (wrect, height), fill=fill
            ))
        group.add(
            canvas.polyline(
                [(DS.gene_arrow_width, 0), (0, height / 2), (DS.gene_arrow_width, height)],
                fill=fill
            ))
        group.add(
            canvas.text(
                label,
                insert=(
                    wrect / 2 + DS.gene_arrow_width,
                    height / 2 + DS.font_central_shift_ratio * DS.label_font_size
                ),
                fill=label_color,
                style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='middle'),
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
                insert=(width / 2, height / 2 + DS.font_central_shift_ratio * DS.label_font_size),
                fill=label_color,
                style=DS.font_style.format(font_size=DS.label_font_size, text_anchor='middle'),
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
