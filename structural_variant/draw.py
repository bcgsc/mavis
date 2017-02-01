"""
.. todo::

    add optional scatter plots to subdiagrams for cna and expression data (possibly using matplotlib?)

"""
import svgwrite
import re
from svgwrite import Drawing
from .interval import Interval
from .constants import STRAND, ORIENT, CODON_SIZE, GIESMA_STAIN
from colour import Color
from .error import DrawingFitError, NotSpecifiedError, DiscontinuousMappingError
from .annotate.genomic import IntergenicRegion
from .annotate.variant import FusionTranscript

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'


class Tag(svgwrite.base.BaseElement):
    def __init__(self, elementname, content='', **kwargs):
        self.elementname = elementname
        super(Tag, self).__init__(**kwargs)
        self.content = content

    def get_xml(self):
        xml = super(Tag, self).get_xml()
        xml.text = self.content
        return xml


class ScatterPlot:
    """
    holds settings that will go into matplotlib after conversion using the mapping system
    """
    def __init__(
        self, points, y_axis_label,
        ymax=None, ymin=None, xmin=None, xmax=None, hmarkers=None, height=100, point_radius=2,
        title='', yticks=None
    ):
        self.hmarkers = hmarkers if hmarkers is not None else []
        self.yticks = yticks if yticks is not None else []
        self.ymin = ymin
        self.ymax = ymax
        self.points = points
        if self.ymin is None:
            self.ymin = min([y.start for x, y in points] + yticks)
        if self.ymax is None:
            self.ymax = max([y.end for x, y in points] + yticks)
        self.xmin = xmin
        self.xmax = xmax
        if self.xmin is None:
            self.xmin = min([x.start for x, y in points])
        if self.xmax is None:
            self.xmax = max([x.end for x, y in points])
        self.y_axis_label = y_axis_label
        self.height = 100
        self.point_radius = 2
        self.title = title


class LabelMapping:
    def __init__(self, **kwargs):
        self._mapping = dict()
        self._reverse_mapping = dict()
        for k, v in kwargs.items():
            self[k] = v

    def __setitem__(self, key, value):
        if key in self._mapping:
            raise KeyError('duplicate key: keys must be unique', key)
        if value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self._mapping[key] = value
        self._reverse_mapping[value] = key

    def items(self):
        return self._mapping.items()

    def __getitem__(self, key):
        return self._mapping[key]

    def __len__(self):
        return len(self._mapping.keys())

    def get_key(self, value):
        return self._reverse_mapping[value]

    def set_key(self, key, value):
        if key in self._mapping:
            current_value = self._mapping[key]
            if value == current_value:
                return
            elif value in self._reverse_mapping:
                raise KeyError('duplicate value: values must be unique', value)
            del self._mapping[key]
            del self._reverse_mapping[current_value]
        elif value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self[key] = value

    def add(self, value, prefix=''):
        if value in self._reverse_mapping:
            return self._reverse_mapping[value]
        i = 1
        while True:
            key = '{}{}'.format(prefix, i)
            if key not in self._mapping:
                self[key] = value
                break
            i += 1
        return self._reverse_mapping[value]


class Diagram:
    """
    class which holds the settings for drawing a fusion digram
    """
    def __init__(
        self,
        WIDTH=1000,
    ):
        self.MIN_WIDTH = 10  # no element (exon, gene, etc can be less than this wide)
        self.TRACK_LINE_HEIGHT = 4
        self.LEFT_MARGIN = 20
        self.RIGHT_MARGIN = 20
        self.TOP_MARGIN = 20
        self.BOTTOM_MARGIN = 20
        self.INNER_MARGIN = 20
        self.PADDING = 5
        self.SCAFFOLD_HEIGHT = 3
        self.SCAFFOLD_COLOR = '#000000'
        self.TRACK_HEIGHT = 50
        self.WIDTH = WIDTH
        # removing unsupported attr: 'alignment-baseline:central;dominant-baseline:central;' \
        self.FONT_STYLE = 'font-size:{font_size}px;font-weight:bold;alignment-baseline:baseline;' \
            'text-anchor:{text_anchor};font-family: consolas, courier new, monospace'
        # ratio for courier new which is wider than consolas, used for estimating width
        self.FONT_WIDTH_HEIGHT_RATIO = 1229 / 2048
        self.FONT_CENTRAL_SHIFT_RATIO = 0.3

        self.GENE1_COLOR_SELECTED = '#518DC5'
        self.GENE2_COLOR_SELECTED = '#4C9677'
        self.GENE1_COLOR = '#657E91'
        self.GENE2_COLOR = '#325556'
        self.GENE_DEFAULT_COLOR = self.GENE1_COLOR
        self.GENE_MIN_BUFFER = 1000
        self.GENE_ARROW_WIDTH = 20
        self.GENE_INTERGENIC_RATIO = 5
        self.GENE_MIN_WIDTH = 40 + self.GENE_ARROW_WIDTH
        self.GENE_LABEL_PREFIX = 'G'

        self.LABEL_COLOR = HEX_BLACK
        self.LABEL_FONT_SIZE = 28
        self.DYNAMIC_LABELS = True
        self.LABEL_LEFT_MARGIN = self.LABEL_FONT_SIZE * self.FONT_WIDTH_HEIGHT_RATIO * 4

        self.DOMAIN_COLOR = '#ccccb3'
        self.DOMAIN_TRACK_HEIGHT = 30
        self.DOMAIN_SCAFFOLD_HEIGHT = 1
        self.DOMAIN_SCAFFOLD_COLOR = HEX_BLACK
        self.DOMAIN_LABEL_PREFIX = 'D'
        self.DOMAIN_LABEL_FONT_SIZE = 20
        self.DOMAIN_MISMATCH_COLOR = '#B2182B'
        self.DOMAIN_FILL_GRADIENT = [
            c.hex for c in Color(self.DOMAIN_MISMATCH_COLOR).range_to(Color(self.DOMAIN_COLOR), 10)]
        self.DOMAIN_NAME_REGEX_FILTER = '.*'

        self.SPLICE_HEIGHT = self.TRACK_HEIGHT / 2
        self.SPLICE_STROKE_DASHARRAY = [2, 2]
        self.SPLICE_STROKE_WIDTH = 2
        self.SPLICE_COLOR = HEX_BLACK

        self.BREAKPOINT_STROKE_DASHARRAY = [3, 3]
        self.BREAKPOINT_ORIENT_STROKE_WIDTH = 2
        self.BREAKPOINT_COLOR = HEX_BLACK
        self.BREAKPOINT_LABEL_FONT_SIZE = 20
        self.BREAKPOINT_BOTTOM_MARGIN = 20
        self.BREAKPOINT_TOP_MARGIN = self.PADDING * 2 + self.BREAKPOINT_LABEL_FONT_SIZE + self.BREAKPOINT_BOTTOM_MARGIN
        self.BREAKPOINT_LABEL_PREFIX = 'B'

        self.MARKER_LABEL_FONT_SIZE = self.BREAKPOINT_LABEL_FONT_SIZE
        self.MARKER_LABEL_PREFIX = 'M'
        self.MARKER_TOP_MARGIN = self.BREAKPOINT_TOP_MARGIN
        self.MARKER_BOTTOM_MARGIN = self.BREAKPOINT_BOTTOM_MARGIN
        self.MARKER_COLOR = self.BREAKPOINT_COLOR

        self.EXON_TEAR_TOOTH_WIDTH = 2
        self.EXON_MIN_WIDTH = self.MIN_WIDTH + self.EXON_TEAR_TOOTH_WIDTH * 2
        self.EXON_TEAR_TOOTH_HEIGHT = 2
        self.EXON_INTRON_RATIO = 20
        self.EXON1_COLOR = self.GENE1_COLOR_SELECTED
        self.EXON2_COLOR = self.GENE2_COLOR_SELECTED
        self.EXON_FONT_SIZE = 20

        self.TRANSCRIPT_LABEL_PREFIX = 'T'
        self.FUSION_LABEL_PREFIX = 'F'

        self.TRANSLATION_FONT_SIZE = 14
        self.TRANSLATION_SCAFFOLD_COLOR = self.SCAFFOLD_COLOR
        self.TRANSLATION_TRACK_HEIGHT = self.TRANSLATION_FONT_SIZE
        self.TRANSLATION_START_MARKER = 'M'
        self.TRANSLATION_END_MARKER = '*'
        self.TRANSLATION_MARKER_PADDING = 4

        self.LEGEND_SWATCH_SIZE = 50
        self.LEGEND_FONT_SIZE = 20
        self.LEGEND_SWATCH_STROKE = HEX_BLACK
        self.LEGEND_FONT_COLOR = HEX_BLACK
        self.LEGEND_BORDER_STROKE = HEX_BLACK
        self.LEGEND_BORDER_STROKE_WIDTH = 1

        self.TEMPLATE_BAND_STROKE_WIDTH = 0.5
        temp = [c.hex for c in Color(HEX_WHITE).range_to(Color(HEX_BLACK), 5)]
        self.TEMPLATE_BAND_FILL = {
            GIESMA_STAIN.ACEN: '#800000',
            GIESMA_STAIN.GPOS25: temp[1],
            GIESMA_STAIN.GPOS50: temp[2],
            GIESMA_STAIN.GPOS75: temp[3],
            GIESMA_STAIN.GPOS100: temp[4],
            GIESMA_STAIN.GNEG: HEX_WHITE
        }
        self.TEMPLATE_BAND_STROKE = HEX_BLACK
        self.TEMPLATE_TRACK_HEIGHT = max([
            self.TRACK_HEIGHT / 3,
            self.LABEL_FONT_SIZE - self.BREAKPOINT_BOTTOM_MARGIN -
            self.BREAKPOINT_TOP_MARGIN + self.BREAKPOINT_LABEL_FONT_SIZE])
        self.TEMPLATE_DEFAULT_FILL = HEX_WHITE
        self.TEMPLATE_BAND_MIN_WIDTH = 2
        self.TEMPLATE_LABEL_PREFIX = 'C'

        self.REGION_LABEL_PREFIX = 'R'
        self.OVERLAY_LEFT_LABEL = 16 * self.FONT_WIDTH_HEIGHT_RATIO * self.EXON_FONT_SIZE

        self.SCATTER_AXIS_FONT_SIZE = 12
        self.SCATTER_ERROR_BAR_STROKE_WIDTH = 1
        self.SCATTER_MARKER_RADIUS = 2
        self.SCATTER_YAXIS_TICK_SIZE = self.PADDING

    def draw_legend(self, canvas, swatches, border=True):
        main_group = canvas.g(class_='legend')
        y = self.PADDING if border else 0
        x = self.PADDING if border else 0
        for swatch, label in swatches:
            g = canvas.g()
            g.add(canvas.rect(
                (0, 0),
                (self.LEGEND_SWATCH_SIZE, self.LEGEND_SWATCH_SIZE),
                fill=swatch,
                stroke=self.LEGEND_SWATCH_STROKE
            ))

            g.add(canvas.text(
                label,
                insert=(self.LEGEND_SWATCH_SIZE + self.PADDING, self.LEGEND_SWATCH_SIZE / 2),
                fill=self.LEGEND_FONT_COLOR,
                style=self.FONT_STYLE.format(text_anchor='start', font_size=self.LEGEND_FONT_SIZE),
                class_='label'
            ))
            g.translate(x, y)
            main_group.add(g)
            y += self.LEGEND_SWATCH_SIZE + self.PADDING

        w = max([len(l) for c, l in swatches]) * self.LEGEND_FONT_SIZE * self.FONT_WIDTH_HEIGHT_RATIO + \
            self.PADDING * (3 if border else 1) + self.LEGEND_SWATCH_SIZE

        if border:
            main_group.add(canvas.rect(
                (0, 0), (w, y), fill='none', stroke=self.LEGEND_BORDER_STROKE,
                stroke_width=self.LEGEND_BORDER_STROKE_WIDTH
            ))
        else:
            y -= self.PADDING
        setattr(main_group, 'height', y)
        setattr(main_group, 'width', w)
        setattr(main_group, 'labels', None)
        setattr(main_group, 'mapping', None)
        return main_group

    def draw(
        self,
        ann,
        fusion_transcript=None,
        REFERENCE_GENOME=None,
        templates=None,
        ignore_absent_templates=True,
        draw_template=True,
        user_friendly_labels=True,
        template_display_label_prefix='c'
    ):
        """
        this is the main drawing function. It decides between the 3 basic layouts

        ::

            1. Breakpoints are in the same transcript

                +---------------------------+
                | Genomic Level             |
                +---------------------------+
                | Transcript Level          |
                +---------------------------+
                | Fusion Level              |
                +---------------------------+

            2. Breakpoints are on the same template/chromosome but different genes

                +---------------------------+
                | Genomic Level             |
                +-------------+-------------+
                | Transcript1 | Transcript2 |
                +-------------+-------------+
                | Fusion Level              |
                +---------------------------+

            3. Breakpoints are on different templates/chromosomes

                +-------------+-------------+
                | Gene1       | Gene2       |
                +-------------+-------------+
                | Transcript1 | Transcript2 |
                +-------------+-------------+
                | Fusion Level              |
                +---------------------------+

        """
        templates = dict() if templates is None else templates
        canvas = Drawing(size=(self.WIDTH, 1000))  # just set the height for now and change later
        labels = LabelMapping()  # keep labels consistent within the drawing
        y = self.TOP_MARGIN
        x = self.LEFT_MARGIN

        dx_label_shift = self.LABEL_LEFT_MARGIN

        x += dx_label_shift
        drawing_width = self.WIDTH - dx_label_shift - self.LEFT_MARGIN - self.RIGHT_MARGIN
        # calculate the half-width for transcripts and genes etc
        half_drawing_width = (drawing_width - self.INNER_MARGIN - dx_label_shift) / 2
        second_drawing_shift = x + half_drawing_width + self.INNER_MARGIN + dx_label_shift

        if draw_template:
            try:
                template1 = templates[ann.transcript1.get_chr()]
                template2 = templates[ann.transcript2.get_chr()]
                if user_friendly_labels and \
                        (len(template_display_label_prefix) + max([len(template1.name), len(template2.name)])) \
                        * self.LEGEND_FONT_SIZE * self.FONT_WIDTH_HEIGHT_RATIO <= self.LABEL_LEFT_MARGIN and \
                        template1.name and template2.name:
                    labels[template_display_label_prefix + template1.name] = template1
                    labels[template_display_label_prefix + template2.name] = template2

                h = [0]
                if template1 == template2:  # single template
                    g = self.draw_template(
                        canvas, template1, drawing_width,
                        breakpoints=[ann.break1, ann.break2], labels=labels)
                    canvas.add(g)
                    g.translate(x, y)
                    h.append(g.height)
                else:  # multiple templates
                    g = self.draw_template(
                        canvas, template1, half_drawing_width, breakpoints=[ann.break1], labels=labels)
                    canvas.add(g)
                    g.translate(x, y)
                    h.append(g.height)

                    g = self.draw_template(
                        canvas, template2, half_drawing_width, breakpoints=[ann.break2], labels=labels)
                    canvas.add(g)
                    g.translate(second_drawing_shift, y)
                    h.append(g.height)
                y += max(h)
            except KeyError as err:
                if not ignore_absent_templates:
                    raise err
                print(repr(err))

        colors = dict()
        genes = set()
        genes1 = set()
        genes2 = set()
        legend = dict()

        for g, d in ann.genes_proximal_to_break1:
            genes1.add(g)
            colors[g] = self.GENE1_COLOR

        for g, d in ann.genes_proximal_to_break2:
            genes2.add(g)
            colors.setdefault(g, self.GENE2_COLOR)

        for g in ann.genes_overlapping_break1:
            genes1.add(g)
            colors[g] = self.GENE1_COLOR

        for g in ann.genes_overlapping_break2:
            genes2.add(g)
            colors.setdefault(g, self.GENE2_COLOR)

        if ann.transcript1:
            try:
                genes1.add(ann.transcript1.gene)
                colors[ann.transcript1.gene] = self.GENE1_COLOR_SELECTED
                for e in ann.transcript1.exons:
                    colors[e] = self.EXON1_COLOR
            except AttributeError:
                genes1.add(ann.transcript1)
                colors[ann.transcript1] = self.GENE1_COLOR_SELECTED

        if ann.transcript2:
            try:
                genes2.add(ann.transcript2.gene)
                colors.setdefault(ann.transcript2.gene, self.GENE2_COLOR_SELECTED)
                for e in ann.transcript2.exons:
                    colors.setdefault(e, self.EXON2_COLOR)
            except AttributeError:
                genes2.add(ann.transcript2)
                colors[ann.transcript2] = self.GENE2_COLOR_SELECTED

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
                l = labels.add(gene, self.REGION_LABEL_PREFIX)
            elif user_friendly_labels and not alias_failure and gene in alias_by_gene:
                labels[alias_by_gene[gene]] = gene
            else:
                l = labels.add(gene, self.GENE_LABEL_PREFIX)

        gheights = [0]

        if ann.interchromosomal:
            # first gene view
            g = self.draw_genes(canvas, genes1, half_drawing_width, [ann.break1], colors=colors, labels=labels)
            g.translate(x, y)
            canvas.add(g)
            gheights.append(g.height)

            # second gene view
            g = self.draw_genes(canvas, genes2, half_drawing_width, [ann.break2], colors=colors, labels=labels)
            g.translate(second_drawing_shift, y)
            canvas.add(g)
            gheights.append(g.height)
        else:
            g = self.draw_genes(canvas, genes1 | genes2, drawing_width, [ann.break1, ann.break2], colors=colors, labels=labels)
            g.translate(x, y)
            canvas.add(g)
            gheights.append(g.height)

        y += max(gheights) + self.INNER_MARGIN

        theights = [0]
        # now the transcript level drawings
        if ann.transcript1 == ann.transcript2:
            try:
                g = canvas.g(class_='transcript')
                g = self.draw_ustranscript(
                    canvas,
                    ann.transcript1,
                    drawing_width,
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
                g = self.draw_ustranscript(
                    canvas,
                    ann.transcript1,
                    half_drawing_width,
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
                g = self.draw_ustranscript(
                    canvas,
                    ann.transcript2,
                    half_drawing_width,
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
            y -= self.INNER_MARGIN

        # finally the fusion transcript level drawing
        if fusion_transcript:
            y += self.INNER_MARGIN
            for ex, old_ex in fusion_transcript.exon_mapping.items():
                if old_ex in colors:
                    colors[ex] = colors[old_ex]
            g = canvas.g(class_='transcript')
            g = self.draw_ustranscript(
                canvas,
                fusion_transcript,
                drawing_width,
                colors=colors,
                labels=labels,
                REFERENCE_GENOME=REFERENCE_GENOME
            )
            g.translate(x, y)
            canvas.add(g)
            y += g.height

        y += self.BOTTOM_MARGIN
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

    def _draw_exon_track(self, canvas, transcript, mapping, colors=None):
        """
        """
        colors = {} if colors is None else colors
        main_group = canvas.g(class_='exon_track')

        y = self.TRACK_HEIGHT / 2
        exons = sorted(transcript.exons, key=lambda x: x.start)

        s = Interval.convert_pos(mapping, exons[0].start)
        t = Interval.convert_pos(mapping, exons[-1].end)

        main_group.add(
            canvas.rect(
                (s, y - self.SCAFFOLD_HEIGHT / 2),
                (t - s + 1, self.SCAFFOLD_HEIGHT),
                fill=self.SCAFFOLD_COLOR,
                class_='scaffold'
            ))

        # draw the exons
        for exon in exons:
            s = Interval.convert_pos(mapping, exon.start)
            t = Interval.convert_pos(mapping, exon.end)
            pxi = Interval(s, t)
            c = colors.get(exon, self.EXON1_COLOR)
            group = self.draw_exon(canvas, exon, pxi.length(), self.TRACK_HEIGHT, c, label=transcript.exon_number(exon))
            group.translate(pxi.start, y - self.TRACK_HEIGHT / 2)
            main_group.add(group)

        setattr(main_group, 'height', y + self.TRACK_HEIGHT / 2)
        setattr(main_group, 'width', t - s + 1)
        return main_group

    def draw_scatter(self, canvas, plot, xmapping):
        """
        given a xmapping, draw the scatter plot svg group
        """
        # generate the y coordinate mapping
        plot_group = canvas.g(class_='scatter_plot')

        yratio = plot.height / (abs(plot.ymax - plot.ymin))
        print('yratio', yratio)
        ypx = []
        xpx = []
        for xp, yp in plot.points:
            try:
                temp = Interval.convert_ratioed_pos(xmapping, xp.start)
                xp = Interval.convert_ratioed_pos(xmapping, xp.end)
                xp = xp | temp
                xpx.append(xp)
                yp = plot.height - abs(yp - plot.ymin) * yratio
                ypx.append(yp)
            except DiscontinuousMappingError:
                pass

        for xp, yp in zip(xpx, ypx):
            if xp.length() > self.SCATTER_MARKER_RADIUS:
                plot_group.add(canvas.line(
                    (xp.start, yp),
                    (xp.end, yp),
                    stroke=HEX_BLACK,
                    stroke_width=self.SCATTER_ERROR_BAR_STROKE_WIDTH
                ))
            plot_group.add(canvas.circle(
                center=(xp.center, yp),
                fill=HEX_BLACK,
                r=self.SCATTER_MARKER_RADIUS
            ))

        for py in plot.hmarkers:
            py = plot.height - abs(py - plot.ymin) * yratio
            plot_group.add(
                canvas.line(
                    start=(min(xpx).start, py),
                    end=(max(xpx).end, py),
                    stroke='blue'
                )
            )
        # draw left y axis
        plot_group.add(canvas.line(
            start=(0, 0), end=(0, plot.height), stroke=HEX_BLACK
        ))
        # draw start and end markers on the y axis
        for y in plot.yticks:
            py = plot.height - abs(y - plot.ymin) * yratio
            plot_group.add(
                canvas.line(
                    start=(0 - self.SCATTER_YAXIS_TICK_SIZE, py),
                    end=(0, py),
                    stroke=HEX_BLACK
                )
            )

        x = 0 - self.PADDING - self.SCATTER_AXIS_FONT_SIZE - self.SCATTER_YAXIS_TICK_SIZE
        y = plot.height / 2 #+ len(plot.y_axis_label) * self.SCATTER_AXIS_FONT_SIZE
        yaxis = canvas.text(
            plot.y_axis_label,
            insert=(x, y),
            fill=self.LABEL_COLOR,
            style=self.FONT_STYLE.format(font_size=self.SCATTER_AXIS_FONT_SIZE, text_anchor='start'),
            class_='y_axis_label'
        )
        print('yaxis', yaxis.tostring())
        plot_group.add(yaxis)
        cx = len(plot.y_axis_label) * self.FONT_WIDTH_HEIGHT_RATIO * self.SCATTER_AXIS_FONT_SIZE / 2
        yaxis.rotate(270, (x + cx, y))
        print('yaxis', yaxis.tostring())
        yaxis.translate(0, 0)
        print('yaxis', yaxis.tostring())

        y = plot.height
        setattr(plot_group, 'height', y)
        return plot_group

    def draw_ustranscript(
        self, canvas, ust,
        target_width=None,
        breakpoints=[],
        labels=LabelMapping(),
        colors={},
        mapping=None,
        REFERENCE_GENOME=None
    ):
        """
        builds an svg group representing the transcript. Exons are drawn in a track with the splicing
        information and domains are drawn in separate tracks below

        if there are mutltiple splicing variants then mutliple exon tracks are drawn

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
            mapping = self._generate_interval_mapping(
                ust.exons,
                target_width,
                self.EXON_INTRON_RATIO,
                self.EXON_MIN_WIDTH,
                min_inter_width=self.MIN_WIDTH
            )

        main_group = canvas.g(class_='ust')
        BTM = self.BREAKPOINT_TOP_MARGIN if len(breakpoints) > 0 else 0
        BBM = self.BREAKPOINT_BOTTOM_MARGIN if len(breakpoints) > 0 else 0

        y = 0

        LABEL_PREFIX = self.TRANSCRIPT_LABEL_PREFIX
        if isinstance(ust, FusionTranscript):
            LABEL_PREFIX = self.FUSION_LABEL_PREFIX

        if len(ust.translations) == 0:
            y += self.SPLICE_HEIGHT + BTM
            exon_track_group = self._draw_exon_track(canvas, ust, mapping, colors)
            exon_track_group.translate(0, y)
            exon_track_group.add(canvas.text(
                labels.add(ust, LABEL_PREFIX),
                insert=(0 - self.PADDING, self.TRACK_HEIGHT / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
                fill=self.LABEL_COLOR,
                style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='end'),
                class_='label'
            ))
            main_group.add(exon_track_group)
            y += self.TRACK_HEIGHT
            y += BBM
        else:
            # draw the protein features if there are any
            for i, tl in enumerate(ust.translations):
                tr = tl.transcript
                # if the splicing takes up more room than the track we need to adjust for it
                y += self.SPLICE_HEIGHT + BTM

                exon_track_group = self._draw_exon_track(canvas, ust, mapping, colors)
                exon_track_group.translate(0, y)
                exon_track_group.add(canvas.text(
                    labels.add(tr, LABEL_PREFIX),
                    insert=(0 - self.PADDING, self.TRACK_HEIGHT / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
                    fill=self.LABEL_COLOR,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='end'),
                    class_='label'
                ))

                # draw the splicing pattern
                splice_group = canvas.g(class_='splicing')
                for p1, p2 in zip(tr.splicing_pattern[::2], tr.splicing_pattern[1::2]):
                    a = Interval.convert_pos(mapping, p1)
                    b = Interval.convert_pos(mapping, p2)
                    polyline = [(a, y), (a + (b - a) / 2, y - self.SPLICE_HEIGHT), (b, y)]
                    p = canvas.polyline(polyline, fill='none')
                    p.dasharray(self.SPLICE_STROKE_DASHARRAY)
                    p.stroke(self.SPLICE_COLOR, width=self.SPLICE_STROKE_WIDTH)
                    splice_group.add(p)

                y += self.TRACK_HEIGHT / 2

                main_group.add(splice_group)
                main_group.add(exon_track_group)
                y += self.TRACK_HEIGHT / 2

                gp = canvas.g(class_='protein')
                # translation track
                y += self.PADDING
                # convert the AA position to cdna position, then convert the cdna to genomic, etc
                # ==================== adding the translation track ============
                translated_genomic_regions = [
                    tr.convert_cdna_to_genomic(tl.start), tr.convert_cdna_to_genomic(tl.end)]
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

                #if len(translated_genomic_regions) == 0:
                #    continue
                s = Interval.convert_pos(mapping, translated_genomic_regions[0].start)
                t = Interval.convert_pos(mapping, translated_genomic_regions[-1].end)

                gt = canvas.g(class_='translation')
                gp.add(gt)
                h = self.TRANSLATION_TRACK_HEIGHT

                for sec in translated_genomic_regions:
                    start = Interval.convert_pos(mapping, sec.start)
                    end = Interval.convert_pos(mapping, sec.end)
                    gt.add(canvas.rect(
                        (start, h / 2 - self.TRANSLATION_TRACK_HEIGHT / 2),
                        (end - start + 1, self.TRANSLATION_TRACK_HEIGHT),
                        fill=self.TRANSLATION_SCAFFOLD_COLOR,
                        class_='scaffold'
                        ))
                gt.add(canvas.text(
                    self.TRANSLATION_END_MARKER if tr.get_strand() == STRAND.NEG else self.TRANSLATION_START_MARKER,
                    insert=(s - self.TRANSLATION_MARKER_PADDING, h / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.TRANSLATION_FONT_SIZE),
                    fill=self.LABEL_COLOR,
                    style=self.FONT_STYLE.format(font_size=self.TRANSLATION_FONT_SIZE, text_anchor='end'),
                    class_='label'
                ))
                gt.add(canvas.text(
                    self.TRANSLATION_START_MARKER if tr.get_strand() == STRAND.NEG else self.TRANSLATION_END_MARKER,
                    insert=(t + self.TRANSLATION_MARKER_PADDING, h / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.TRANSLATION_FONT_SIZE),
                    fill=self.LABEL_COLOR,
                    style=self.FONT_STYLE.format(font_size=self.TRANSLATION_FONT_SIZE, text_anchor='start'),
                    class_='label'
                ))
                gt.add(Tag('title', 'cds: {}_{}; aa length: {}'.format(tl.start, tl.end, len(tl) // CODON_SIZE)))
                py = h
                # now draw the domain tracks
                # need to convert the domain AA positions to cds positions to genomic
                for i, d in enumerate(sorted(tl.domains, key=lambda x: x.name)):
                    if not re.match(self.DOMAIN_NAME_REGEX_FILTER, str(d.name)):
                        continue
                    py += self.PADDING
                    gd = canvas.g(class_='domain')
                    gd.add(canvas.rect(
                        (0, self.DOMAIN_TRACK_HEIGHT / 2),
                        (target_width, self.DOMAIN_SCAFFOLD_HEIGHT),
                        fill=self.DOMAIN_SCAFFOLD_COLOR,
                        class_='scaffold'
                    ))
                    fill = self.DOMAIN_COLOR
                    percent_match = None
                    try:
                        match, total = d.score_region_mapping(REFERENCE_GENOME)
                        percent_match = int(round(match * 100 / total, 0))
                        fill = self.DOMAIN_FILL_GRADIENT[percent_match % len(self.DOMAIN_FILL_GRADIENT) - 1]
                    except (NotSpecifiedError, AttributeError):
                        pass
                    for region in d.regions:
                        # convert the AA position to cdna position, then convert the cdna to genomic, etc
                        s = tl.convert_aa_to_cdna(region.start)
                        t = tl.convert_aa_to_cdna(region.end)
                        s = tr.convert_cdna_to_genomic(s.start)
                        t = tr.convert_cdna_to_genomic(t.end)
                        s = Interval.convert_pos(mapping, s)
                        t = Interval.convert_pos(mapping, t)
                        if s > t:
                            t, s = (s, t)
                        gdr = canvas.g(class_='domain_region')
                        gdr.add(canvas.rect(
                            (s, 0), (t - s + 1, self.DOMAIN_TRACK_HEIGHT),
                            fill=fill, class_='region'))
                        gdr.add(Tag('title', 'Domain {} region {}_{}aa{}'.format(
                            d.name if d.name else '', region.start, region.end,
                            ' matched {}%'.format(percent_match) if percent_match is not None else '')))
                        gd.add(gdr)
                    gd.translate(0, py)

                    gd.add(canvas.text(
                        labels.add(d.name, self.DOMAIN_LABEL_PREFIX),
                        insert=(
                            0 - self.PADDING,
                            self.DOMAIN_TRACK_HEIGHT / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.DOMAIN_LABEL_FONT_SIZE),
                        fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(self.DOMAIN_COLOR),
                        style=self.FONT_STYLE.format(font_size=self.DOMAIN_LABEL_FONT_SIZE, text_anchor='end'),
                        class_='label'
                    ))
                    gp.add(gd)
                    py += self.DOMAIN_TRACK_HEIGHT
                gp.translate(0, y)
                if i < len(ust.translations):
                    y += self.INNER_MARGIN
                y += py + BBM
                main_group.add(gp)

        # now overlay the breakpoints on top of everything
        for i, b in enumerate(breakpoints):
            s = Interval.convert_ratioed_pos(mapping, b.start)
            t = Interval.convert_ratioed_pos(mapping, b.end)
            px_itvl = Interval(s.start, t.end)
            bg = self.draw_breakpoint(canvas, b, px_itvl.length(), y, label=labels.add(b, self.BREAKPOINT_LABEL_PREFIX))
            bg.translate(px_itvl.start, 0)
            main_group.add(bg)

        setattr(main_group, 'height', y)
        setattr(main_group, 'width', target_width)
        setattr(main_group, 'mapping', mapping)
        setattr(main_group, 'labels', labels)
        return main_group

    def draw_genes(self, canvas, genes, target_width, breakpoints=None, colors=None, labels=None):
        """
        draws the genes given in order of their start position trying to minimize
        the number of tracks required to avoid overlap

        Args:
            canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
            target_width (int): the target width of the diagram
            genes (:class:`list` of :class:`Gene`): the list of genes to draw
            breakpoints (:class:`list` of :class:`Breakpoint`): the breakpoints to overlay
            colors (:class:`dict` of :class:`Gene` and :class:`str`): dictionary of the colors assigned to each Gene as fill

        Return:
            svgwrite.container.Group: the group element for the diagram.
                Has the added parameters of labels, height, and mapping
        """
        # mutable default argument parameters
        breakpoints = [] if breakpoints is None else breakpoints
        colors = {} if colors is None else colors
        labels = LabelMapping() if labels is None else labels

        st = max(min([g.start for g in genes] + [b.start for b in breakpoints]) - self.GENE_MIN_BUFFER, 1)
        end = max([g.end for g in genes] + [b.end for b in breakpoints]) + self.GENE_MIN_BUFFER
        main_group = canvas.g(class_='genes')
        mapping = self._generate_interval_mapping(
            [g for g in genes],
            target_width,
            self.GENE_INTERGENIC_RATIO,
            self.GENE_MIN_WIDTH,
            start=st, end=end,
            min_inter_width=self.MIN_WIDTH
        )
        gene_px_intervals = {}
        for i, gene in enumerate(sorted(genes, key=lambda x: x.start)):
            s = Interval.convert_ratioed_pos(mapping, gene.start)
            t = Interval.convert_ratioed_pos(mapping, gene.end)
            gene_px_intervals[Interval(s.start, t.end)] = gene
            l = labels.add(gene, self.GENE_LABEL_PREFIX)
        tracks = Diagram._split_intervals_into_tracks(gene_px_intervals)

        y = self.BREAKPOINT_TOP_MARGIN

        main_group.add(
            canvas.rect(
                (0, y + self.TRACK_HEIGHT / 2 - self.SCAFFOLD_HEIGHT / 2
                    + (len(tracks) - 1) * (self.TRACK_HEIGHT + self.PADDING)),
                (target_width, self.SCAFFOLD_HEIGHT),
                fill=self.SCAFFOLD_COLOR,
                class_='scaffold'
            ))
        tracks.reverse()
        for track in tracks:  # svg works from top down
            for genepx in track:
                # draw the gene
                gene = gene_px_intervals[genepx]
                group = self.draw_gene(
                    canvas, gene, genepx.length(),
                    self.TRACK_HEIGHT,
                    colors.get(gene, self.GENE1_COLOR),
                    labels.get_key(gene)
                )
                group.translate(genepx.start, y)
                main_group.add(group)
            y += self.TRACK_HEIGHT + self.PADDING

        y += self.BREAKPOINT_BOTTOM_MARGIN - self.PADDING
        # now overlay the breakpoints on top of everything
        for i, b in enumerate(sorted(breakpoints)):
            s = Interval.convert_pos(mapping, b.start)
            t = Interval.convert_pos(mapping, b.end)
            bg = self.draw_breakpoint(canvas, b, abs(t - s) + 1, y, label=labels.add(b, self.BREAKPOINT_LABEL_PREFIX))
            bg.translate(s, 0)
            main_group.add(bg)

        setattr(main_group, 'height', y)
        setattr(main_group, 'width', target_width)
        setattr(main_group, 'mapping', mapping)
        setattr(main_group, 'labels', labels)

        return main_group

    def draw_area_plot(self, canvas, data, height, fill, xmapping=None, **kwargs):
        ymin = kwargs.pop('ymin', min([y for x, y in data]))
        ymax = kwargs.pop('ymax', max([y for x, y in data]))
        yratio = height / (ymax - ymin + 1)

        if kwargs:
            raise AttributeError('invalid arguments', kwargs)

        temp = []
        for x, y in data:
            x = Interval.convert_pos(xmapping, x) if xmapping is not None else x
            temp.append((x, (ymax - y) * yratio))
        data = temp

        g = canvas.g(class_='area_plot')

        data = sorted(data)
        polyline = [(data[0][0], ymax * yratio)] + data + [(data[-1][0], ymax * yratio)]
        p = canvas.polyline(polyline, fill=fill)
        g.add(p)

        h = max([y for x, y in data])
        setattr(g, 'height', h)
        return g

    def read_config(self):
        pass

    def draw_ustranscripts_overlay(self, gene, vmarkers=None, window_buffer=0, plots=None):
        vmarkers = [] if vmarkers is None else vmarkers
        plots = [] if plots is None else plots

        canvas = Drawing(size=(self.WIDTH, 1000))  # just set the height for now and change later
        w = self.WIDTH - self.LEFT_MARGIN - self.RIGHT_MARGIN - self.OVERLAY_LEFT_LABEL - self.PADDING
        labels = LabelMapping()  # keep labels consistent within the drawing

        all_exons = set()
        colors = dict()
        for tx in gene.transcripts:
            for ex in tx.exons:
                all_exons.add(ex)
                colors[ex] = self.EXON1_COLOR if tx.is_best_transcript else self.EXON2_COLOR

        st = min([max([gene.start - window_buffer, 1])] + [m.start for m in vmarkers] + [p.xmin for p in plots])
        end = max([gene.end + window_buffer] + [m.end for m in vmarkers] + [p.xmax for p in plots])

        mapping = Diagram._generate_interval_mapping(
            all_exons,
            w,
            self.EXON_INTRON_RATIO,
            self.MIN_WIDTH,
            start=st, end=end
        )
        main_group = canvas.g(class_='overlay')

        x = self.OVERLAY_LEFT_LABEL + self.PADDING
        y = self.MARKER_TOP_MARGIN

        for plot in plots:
            plot_group = self.draw_scatter(canvas, plot, mapping)
            main_group.add(plot_group)
            plot_group.translate(x, y)
            y += plot.height + self.PADDING

        for tx in gene.transcripts:
            g = Diagram._draw_exon_track(self, canvas, tx, mapping, colors=colors)
            main_group.add(g)
            g.translate(x, y)

            t = canvas.text(
                tx.name,
                insert=(x - self.PADDING, y + self.TRACK_HEIGHT / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.EXON_FONT_SIZE),
                fill=self.LABEL_COLOR,
                style=self.FONT_STYLE.format(font_size=self.EXON_FONT_SIZE, text_anchor='end'),
                class_='label'
            )
            main_group.add(t)
            y += self.PADDING + self.TRACK_HEIGHT

        y += self.MARKER_BOTTOM_MARGIN
        # now draw the breakpoints overtop
        for i, m in enumerate(sorted(vmarkers)):
            s = Interval.convert_ratioed_pos(mapping, m.start)
            t = Interval.convert_ratioed_pos(mapping, m.end)
            px_itvl = Interval(s.start, t.end)
            bg = self.draw_vmarker(
                canvas, m, px_itvl.length(), y, label=labels.add(m, self.MARKER_LABEL_PREFIX))
            bg.translate(x + px_itvl.start, 0)
            main_group.add(bg)

        main_group.translate(self.LEFT_MARGIN, self.TOP_MARGIN)
        y += self.BOTTOM_MARGIN
        canvas.add(main_group)
        canvas.attribs['height'] = y
        return canvas

    def draw_vmarker(self, canvas, marker, width, height, label='', color=None):
        """
        Args:
            canvas (svgwrite.drawing.Drawing): the main svgwrite object used to create new svg elements
            breakpoint (Breakpoint): the breakpoint to draw
            width (int): the pixel width
            height (int): the pixel height
        Return:
            svgwrite.container.Group: the group element for the diagram
        """
        color = self.MARKER_COLOR if color is None else color
        g = canvas.g(class_='marker')
        y = self.PADDING + self.MARKER_LABEL_FONT_SIZE / 2
        t = canvas.text(
            label,
            insert=(width / 2, y + self.FONT_CENTRAL_SHIFT_RATIO * self.MARKER_LABEL_FONT_SIZE),
            fill=color,
            style=self.FONT_STYLE.format(text_anchor='middle', font_size=self.MARKER_LABEL_FONT_SIZE),
            class_='label'
        )
        y += self.MARKER_LABEL_FONT_SIZE / 2 + self.PADDING

        g.add(t)

        r = canvas.rect(
            (0, y),
            (width, height - y),
            stroke=color,
            fill='none'
        )
        g.add(r)

        g.add(Tag('title', 'marker {}:{}-{} {}'.format(
            marker.reference_object, marker.start, marker.end, marker.name)))
        return g

    def draw_breakpoint(self, canvas, breakpoint, width, height, label=''):
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
        y = self.PADDING + self.BREAKPOINT_LABEL_FONT_SIZE / 2
        t = canvas.text(
            label,
            insert=(width / 2, y + self.FONT_CENTRAL_SHIFT_RATIO * self.BREAKPOINT_LABEL_FONT_SIZE),
            fill=HEX_BLACK,
            style=self.FONT_STYLE.format(text_anchor='middle', font_size=self.BREAKPOINT_LABEL_FONT_SIZE),
            class_='label'
        )
        y += self.BREAKPOINT_LABEL_FONT_SIZE / 2 + self.PADDING

        g.add(t)

        r = canvas.rect(
            (0, y),
            (width, height - y),
            stroke=self.BREAKPOINT_COLOR,
            fill='none'
        )
        r.dasharray(self.BREAKPOINT_STROKE_DASHARRAY)
        g.add(r)

        if breakpoint.orient == ORIENT.LEFT:
            l = canvas.line((0, y), (0, height))
            l.stroke(self.BREAKPOINT_COLOR, width=self.BREAKPOINT_ORIENT_STROKE_WIDTH)
            g.add(l)
        elif breakpoint.orient == ORIENT.RIGHT:
            l = canvas.line((width, y), (width, height))
            l.stroke(self.BREAKPOINT_COLOR, width=self.BREAKPOINT_ORIENT_STROKE_WIDTH)
            g.add(l)
        g.add(Tag('title', 'Breakpoint {}:{}-{} strand={} orient={}'.format(
            breakpoint.chr, breakpoint.start, breakpoint.end, breakpoint.strand, breakpoint.orient)))
        return g

    def draw_exon(self, canvas, exon, width, height, fill, label='', tear_left=False, tear_right=False):
        """
        generates the svg object representing an exon

        ::

            intact exon

            +-----+
            |     |
            +-----+

            exon "torn" on the right side (abrogated 3' splice site if on the positive strand)

            +----->
            |     >
            +----->

            exon "torn" on the left side (abrogated 5' splice site if on the positive strand)

            <-----+
            <     |
            <-----+


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
        if tear_right and tear_left:
            raise NotImplementedError('have not added support for tearing exons yet')
        elif tear_left:
            raise NotImplementedError('have not added support for tearing exons yet')
        elif tear_right:
            raise NotImplementedError('have not added support for tearing exons yet')
        else:
            g.add(canvas.rect((0, 0), (width, height), fill=fill))
            t = canvas.text(
                label,
                insert=(width / 2, height / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.EXON_FONT_SIZE),
                fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(fill),
                style=self.FONT_STYLE.format(font_size=self.EXON_FONT_SIZE, text_anchor='middle'),
                class_='label'
            )
            g.add(t)
            g.add(Tag('title', 'Exon {} {}_{} L={}'.format(exon.name if exon.name else '', exon.start, exon.end, len(exon))))
        return g

    def draw_template(self, canvas, template, target_width, labels=None, colors=None, breakpoints=None):
        labels = LabelMapping() if labels is None else labels
        colors = {} if colors is None else colors
        breakpoints = [] if not breakpoints else breakpoints
        total_height = self.TEMPLATE_TRACK_HEIGHT + self.BREAKPOINT_TOP_MARGIN + self.BREAKPOINT_BOTTOM_MARGIN
        group = canvas.g(class_='template')
        mapping = self._generate_interval_mapping(
            template.bands,
            target_width,
            1,  # do not alter ratio
            self.TEMPLATE_BAND_MIN_WIDTH,
            start=template.start, end=template.end
        )
        scaffold = canvas.rect(
            (0, 0),
            (target_width, self.SCAFFOLD_HEIGHT),
            fill=self.SCAFFOLD_COLOR
        )
        group.add(scaffold)
        scaffold.translate((0, self.BREAKPOINT_TOP_MARGIN + self.TEMPLATE_TRACK_HEIGHT / 2 - self.SCAFFOLD_HEIGHT / 2))
        group.add(canvas.text(
            labels.add(template, self.TEMPLATE_LABEL_PREFIX),
            insert=(0 - self.PADDING, self.BREAKPOINT_TOP_MARGIN + self.TEMPLATE_TRACK_HEIGHT / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
            fill=self.LABEL_COLOR,
            style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='end'),
            class_='label'
        ))

        for band in template.bands:
            s = Interval.convert_pos(mapping, band[0])
            t = Interval.convert_pos(mapping, band[1])

            bgroup = canvas.g(class_='cytoband')
            f = self.TEMPLATE_BAND_FILL.get(band.data.get('giesma_stain', None), self.TEMPLATE_DEFAULT_FILL)
            r = None
            w = t - s + 1
            if band.data.get('giesma_stain', None) == GIESMA_STAIN.ACEN:
                if band.name[0] == 'p':
                    r = canvas.polyline(
                        [(0, 0), (w, self.TEMPLATE_TRACK_HEIGHT / 2), (0, self.TEMPLATE_TRACK_HEIGHT)],
                        fill=f, stroke=self.TEMPLATE_BAND_STROKE, stroke_width=self.TEMPLATE_BAND_STROKE_WIDTH
                    )
                else:
                    r = canvas.polyline(
                        [(w, 0), (0, self.TEMPLATE_TRACK_HEIGHT / 2), (w, self.TEMPLATE_TRACK_HEIGHT)],
                        fill=f, stroke=self.TEMPLATE_BAND_STROKE, stroke_width=self.TEMPLATE_BAND_STROKE_WIDTH
                    )
            else:
                r = canvas.rect(
                    (0, 0),
                    (w, self.TEMPLATE_TRACK_HEIGHT),
                    fill=f, stroke=self.TEMPLATE_BAND_STROKE, stroke_width=self.TEMPLATE_BAND_STROKE_WIDTH
                )
            bgroup.add(r)
            bgroup.add(Tag('title', 'template {}: band {} {}-{} L={}M'.format(
                template.name, band.name, band.start, band.end, round(len(band)/1000000, 0))))
            bgroup.translate((s, self.BREAKPOINT_TOP_MARGIN))
            group.add(bgroup)
        # now draw the breakpoints overtop
        for i, b in enumerate(sorted(breakpoints)):
            s = Interval.convert_pos(mapping, b.start)
            t = Interval.convert_pos(mapping, b.end)
            bg = self.draw_breakpoint(
                canvas, b, abs(t - s) + 1, total_height, label=labels.add(b, self.BREAKPOINT_LABEL_PREFIX))
            bg.translate(s, 0)
            group.add(bg)
        setattr(
            group, 'height', total_height)
        return group

    def draw_gene(self, canvas, gene, width, height, fill, label='', REFERENCE_GENOME=None):
        """
        generates the svg object representing a gene

        ::

            gene on the positive/forward strand

            +-----\\
            +-----/

            gene on the negative/reverse strand

            /-----+
            \-----+

            gene with a non-specified strand

            +-----+
            +-----+

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
        if width < self.GENE_MIN_WIDTH:
            raise DrawingFitError('width of {} is not sufficient to draw a gene of minimum width {}'.format(width, self.GENE_MIN_WIDTH))
        wrect = width - self.GENE_ARROW_WIDTH
        if wrect < 1:
            raise DrawingFitError('width is not sufficient to draw gene')

        label_color = self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(fill)

        if gene.get_strand() == STRAND.POS:
            group.add(
                canvas.rect(
                    (0, 0),
                    (wrect, height),
                    fill=fill
                ))
            group.add(
                canvas.polyline(
                    [
                        (wrect, 0),
                        (wrect + self.GENE_ARROW_WIDTH, height / 2),
                        (wrect, height)
                    ],
                    fill=fill
                ))
            group.add(
                canvas.text(
                    label,
                    insert=(wrect / 2, height / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
                    fill=label_color,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                    class_='label'
                ))
        elif gene.get_strand() == STRAND.NEG:
            group.add(
                canvas.rect(
                    (self.GENE_ARROW_WIDTH, 0),
                    (wrect, height),
                    fill=fill
                ))
            group.add(
                canvas.polyline(
                    [
                        (self.GENE_ARROW_WIDTH, 0),
                        (0, height / 2),
                        (self.GENE_ARROW_WIDTH, height)
                    ], fill=fill
                ))
            group.add(
                canvas.text(
                    label,
                    insert=(wrect / 2 + self.GENE_ARROW_WIDTH, height / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
                    fill=label_color,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                    class_='label'
                ))
        else:
            group.add(
                canvas.rect(
                    (0, 0),
                    (width, height),
                    fill=fill
                ))
            group.add(
                canvas.text(
                    label,
                    insert=(width / 2, height / 2 + self.FONT_CENTRAL_SHIFT_RATIO * self.LABEL_FONT_SIZE),
                    fill=label_color,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                    class_='label'
                ))
        aliases = ''
        try:
            if len(gene.aliases) > 0:
                aliases = ' aka {}'.format(';'.join(sorted(gene.aliases)))
        except AttributeError:
            pass
        group.add(
            Tag('title', 'Gene {} {}:{}_{}{}{}'.format(gene.name if gene.name else '',
                gene.chr, gene.start, gene.end, gene.get_strand(), aliases)))
        return group

    @classmethod
    def dynamic_label_color(cls, color):
        """
        calculates the luminance of a color and determines if a black or white label will be more contrasting
        """
        f = Color(color)
        if f.get_luminance() < 0.5:
            return HEX_WHITE
        else:
            return HEX_BLACK

    @classmethod
    def _split_intervals_into_tracks(cls, intervals):
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

    @classmethod
    def _generate_interval_mapping(
        cls, input_intervals, target_width, ratio, min_width,
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
        for itvl in input_intervals:
            p1 = Interval.convert_ratioed_pos(mapping, itvl.start)
            p2 = Interval.convert_ratioed_pos(mapping, itvl.end)
            n = p1 | p2
            if n.length() < min_width:
                raise AssertionError(
                    'interval mapping should not map any intervals to less than the minimum required width. Interval {}'
                    ' was mapped to a pixel interval of length {} but the minimum width is {}'.format(
                        itvl, n.length(), min_width), p1, p2, mapping, input_intervals, target_width, ratio, min_inter_width)
        return mapping
