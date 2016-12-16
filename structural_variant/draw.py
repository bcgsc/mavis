from svgwrite import Drawing
from structural_variant.interval import Interval
from structural_variant.constants import STRAND, ORIENT
from colour import Color

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'


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
        BREAKPOINT_COLOR=HEX_BLACK,
        BREAKPOINT_LABEL_FONT_SIZE=20,
        GENE1_COLOR_SELECTED='#4C9677',
        GENE2_COLOR_SELECTED='#518DC5',
        GENE1_COLOR='#325556',
        GENE2_COLOR='#657E91',
        GENE_DEFAULT_COLOR='#325556',
        LABEL_COLOR=HEX_BLACK,
        LABEL_FONT_SIZE=28,
        EXON1_COLOR='#4C9677',
        EXON2_COLOR='#518DC5',
        EXON1_UTR_COLOR='#7bddc1',
        EXON2_UTR_COLOR='#7dc3d8',
        DOMAIN_COLOR='#b8d3ba',
        DOMAIN_SCAFFOLD_COLOR=HEX_BLACK,
        SPLICE_COLOR=HEX_BLACK,
        WIDTH=1000,
        PADDING=5,
        DYNAMIC_LABELS=True,
        LEGEND_FONT_SIZE=20,
        LEGEND_BORDER_STROKE=HEX_BLACK,
        DOMAIN_LABEL_FONT_SIZE=20
    ):
        self.MIN_WIDTH = 2  # no element (exon, gene, etc can be less than this wide)
        self.TRACK_LINE_HEIGHT = 4
        self.LEFT_MARGIN = 20
        self.RIGHT_MARGIN = 20
        self.TOP_MARGIN = 20
        self.BOTTOM_MARGIN = 20
        self.INNER_MARGIN = 20
        self.PADDING = PADDING
        self.SCAFFOLD_HEIGHT = 3
        self.SCAFFOLD_COLOR = '#000000'
        self.TRACK_HEIGHT = 50
        self.WIDTH = WIDTH
        self.FONT_STYLE = 'font-size:{font_size}px;font-weight:bold;alignment-baseline:central;' \
            'text-anchor:{text_anchor};font-family: consolas, courier new, monospace'
        # ratio for courier new which is wider than consolas, used for estimating width
        self.FONT_WIDTH_HEIGHT_RATIO = 1229 / 2048

        self.GENE1_COLOR_SELECTED = GENE1_COLOR_SELECTED
        self.GENE2_COLOR_SELECTED = GENE2_COLOR_SELECTED
        self.GENE1_COLOR = GENE1_COLOR
        self.GENE2_COLOR = GENE2_COLOR
        self.GENE_DEFAULT_COLOR = GENE_DEFAULT_COLOR
        self.GENE_MIN_BUFFER = 200
        self.GENE_ARROW_WIDTH = 20
        self.GENE_INTERGENIC_RATIO = 5
        self.GENE_MIN_WIDTH = self.MIN_WIDTH + self.GENE_ARROW_WIDTH
        self.GENE_LABEL_PREFIX = 'G'

        self.LABEL_COLOR = LABEL_COLOR
        self.LABEL_FONT_SIZE = LABEL_FONT_SIZE
        self.DYNAMIC_LABELS = DYNAMIC_LABELS

        self.DOMAIN_COLOR = DOMAIN_COLOR
        self.DOMAIN_TRACK_HEIGHT = 30
        self.DOMAIN_SCAFFOLD_HEIGHT = 1
        self.DOMAIN_SCAFFOLD_COLOR = DOMAIN_SCAFFOLD_COLOR
        self.DOMAIN_LABEL_PREFIX = 'D'
        self.DOMAIN_LABEL_FONT_SIZE = DOMAIN_LABEL_FONT_SIZE

        self.SPLICE_HEIGHT = self.TRACK_HEIGHT
        self.SPLICE_STROKE_DASHARRAY = [2, 2]
        self.SPLICE_STROKE_WIDTH = 2
        self.SPLICE_COLOR = SPLICE_COLOR

        self.BREAKPOINT_STROKE_DASHARRAY = [3, 3]
        self.BREAKPOINT_ORIENT_STROKE_WIDTH = 2
        self.BREAKPOINT_COLOR = BREAKPOINT_COLOR
        self.BREAKPOINT_LABEL_FONT_SIZE = BREAKPOINT_LABEL_FONT_SIZE
        self.BREAKPOINT_BOTTOM_MARGIN = 20
        self.BREAKPOINT_TOP_MARGIN = self.PADDING * 2 + self.BREAKPOINT_LABEL_FONT_SIZE + self.BREAKPOINT_BOTTOM_MARGIN
        self.BREAKPOINT_LABEL_PREFIX = 'B'

        self.EXON_TEAR_TOOTH_WIDTH = 2
        self.EXON_MIN_WIDTH = self.MIN_WIDTH + self.EXON_TEAR_TOOTH_WIDTH * 2
        self.EXON_TEAR_TOOTH_HEIGHT = 2
        self.EXON_INTRON_RATIO = 5
        self.EXON1_COLOR = EXON1_COLOR
        self.EXON2_COLOR = EXON2_COLOR
        self.EXON1_UTR_COLOR = EXON1_UTR_COLOR
        self.EXON2_UTR_COLOR = EXON2_UTR_COLOR

        self.TRANSLATION_FONT_SIZE = 14
        self.TRANSLATION_SCAFFOLD_COLOR = self.SCAFFOLD_COLOR
        self.TRANSLATION_SCAFFOLD_HEIGHT = self.TRANSLATION_FONT_SIZE
        self.TRANSLATION_START_MARKER = 'M'
        self.TRANSLATION_END_MARKER = '*'

        self.LEGEND_SWATCH_SIZE = 50
        self.LEGEND_FONT_SIZE = LEGEND_FONT_SIZE
        self.LEGEND_SWATCH_STROKE = HEX_BLACK
        self.LEGEND_FONT_COLOR = HEX_BLACK
        self.LEGEND_BORDER_STROKE = LEGEND_BORDER_STROKE
        self.LEGEND_BORDER_STROKE_WIDTH = 1

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

    def draw(self, ann, fusion_transcript=None):
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

        TODO: fusion level
        """
        canvas = Drawing(height=1000, width=self.WIDTH)  # just set the height for now and change later
        labels = LabelMapping()  # keep labels consistent within the drawing
        y = self.TOP_MARGIN
        x = self.LEFT_MARGIN

        dx_label_shift = self.PADDING + self.LABEL_FONT_SIZE * 2 * self.FONT_WIDTH_HEIGHT_RATIO

        x += dx_label_shift
        drawing_width = self.WIDTH - dx_label_shift - self.LEFT_MARGIN - self.RIGHT_MARGIN

        if ann.interchromosomal:  # two gene tracks
            raise NotImplementedError('does not support two gene tracks yet')
        elif (ann.transcript1 or ann.transcript1) and ann.transcript1 == ann.transcript2:  # single gene track
            colors = {}
            genes = set() | ann.encompassed_genes

            for g in ann.encompassed_genes:
                genes.add(g)
                colors[g] = self.GENE2_COLOR

            for g, d in ann.nearest_gene_break1 | ann.nearest_gene_break2:
                genes.add(g)

            for gene in ann.genes_at_break1 | ann.genes_at_break2:
                genes.add(gene)
                colors[gene] = self.GENE1_COLOR

            if ann.transcript1:
                genes.add(ann.transcript1.gene)
                colors[ann.transcript1.gene] = self.GENE1_COLOR_SELECTED

            for gene in sorted(ann.genes_at_break1, key=lambda x: x.start):
                labels.add(gene, self.GENE_LABEL_PREFIX)
            g = self.draw_genes(canvas, genes, drawing_width, breakpoints=[ann.break1, ann.break2], colors=colors, labels=labels)
            g.translate(x, y)
            canvas.add(g)
            y += g.height + self.INNER_MARGIN

            for exon in ann.transcript1.exons:
                colors[exon] = self.EXON1_COLOR

            # now draw the transcript track
            g = self.draw_transcript(
                canvas,
                ann.transcript1,
                drawing_width,
                breakpoints=[ann.break1, ann.break2],
                labels=labels,
                colors=colors
            )
            g.translate(x, y)
            canvas.add(g)
            y += g.height
            # draw the genes in order
        elif ann.transcript1 or ann.transcript1:
            print('one gene diagram; two transcript diagrams')
            w = self.WIDTH - self.LEFT_MARGIN - self.RIGHT_MARGIN
            colors = {}
            genes = set() | ann.encompassed_genes

            for g, d in ann.nearest_gene_break1:
                genes.add(g)
                colors[g] = self.GENE1_COLOR

            for g, d in ann.nearest_gene_break2:
                genes.add(g)
                colors[g] = self.GENE2_COLOR

            for gene in ann.genes_at_break1:
                genes.add(gene)
                colors[gene] = self.GENE1_COLOR

            for gene in ann.genes_at_break2:
                genes.add(gene)
                colors[gene] = self.GENE2_COLOR

            if ann.transcript1:
                genes.add(ann.transcript1.gene)
                colors[ann.transcript1.gene] = self.GENE1_COLOR_SELECTED
                for e in ann.transcript1.exons:
                    colors[e] = self.EXON1_COLOR
            if ann.transcript2:
                genes.add(ann.transcript2.gene)
                colors[ann.transcript2.gene] = self.GENE2_COLOR_SELECTED
                for e in ann.transcript2.exons:
                    colors[e] = self.EXON2_COLOR

            for gene in sorted(ann.genes_at_break1, key=lambda x: x.start):
                labels.add(gene, self.GENE_LABEL_PREFIX)

            g = self.draw_genes(canvas, genes, drawing_width, breakpoints=[ann.break1, ann.break2], colors=colors, labels=labels)
            g.translate(x, y)
            canvas.add(g)
            y += g.height + self.INNER_MARGIN

            tw = (w - self.INNER_MARGIN) / 2

            h = [0]

            twidth = (drawing_width - self.INNER_MARGIN - dx_label_shift) / 2

            if ann.transcript1:
                g = canvas.g(class_='transcript')
                g = self.draw_transcript(
                    canvas,
                    ann.transcript1,
                    twidth,
                    breakpoints=[ann.break1],
                    labels=labels,
                    colors=colors
                )
                h.append(g.height)
                g.translate(x, y)
                canvas.add(g)

            if ann.transcript2:
                g = canvas.g(class_='transcript')
                g = self.draw_transcript(
                    canvas,
                    ann.transcript2,
                    twidth,
                    breakpoints=[ann.break2],
                    labels=labels,
                    colors=colors
                )
                h.append(g.height)
                g.translate(x + twidth + self.INNER_MARGIN + dx_label_shift, y)
                canvas.add(g)

            y += max(h)

            if fusion_transcript:
                y += self.PADDING
                for ex, old_ex in fusion_transcript.exon_mapping.items():
                    if old_ex in colors:
                        colors[ex] = colors[old_ex]
                g = canvas.g(class_='transcript')
                g = self.draw_transcript(
                    canvas,
                    fusion_transcript,
                    drawing_width,
                    colors=colors,
                    labels=labels
                )
                g.translate(x, y)
                canvas.add(g)
                y += g.height

        y += self.BOTTOM_MARGIN
        canvas.height = y
        return canvas

        # add transcript tracks if applicable (and domains)
        # add fusion track if applicable

    def _draw_exon_track(self, canvas, transcript, mapping, colors=None, labels=None):
        """
        """
        labels = LabelMapping() if labels is None else labels
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
            group = self.draw_exon(canvas, exon, len(pxi), self.TRACK_HEIGHT, c, label=transcript.exon_number(exon))
            group.translate(pxi.start, y - self.TRACK_HEIGHT / 2)
            main_group.add(group)

        setattr(main_group, 'height', y + self.TRACK_HEIGHT / 2)
        setattr(main_group, 'width', t - s + 1)
        return main_group

    def draw_transcript(
        self, canvas, transcript, target_width=None, breakpoints=[], labels=LabelMapping(), colors={}, mapping=None
    ):
        """
        builds an svg group representing the transcript. Exons are drawn in a track with the splicing
        information and domains are drawn in separate tracks below

        if there are mutltiple splicing variants then mutliple exon tracks are drawn

        Args:
            canvas (svgwrite.Drawing): the main svgwrite object used to create new svg elements
            target_width (int): the target width of the diagram
            t (structural_variant.annotate.Transcript): the transcript being drawn
            exon_color (str): the color being used for the fill of the exons
            utr_color (str): the color for the fill of the UTR regions
            abrogated_splice_sites (list of int): list of positions to ignore as splice sites
            breakpoints (iterable of Breakpoint): the breakpoints to overlay

        Return:
            svgwrite.container.Group: the group element for the transcript diagram
                    Has the added parameters of labels, height, and mapping
        """
        if transcript.strand not in [STRAND.POS, STRAND.NEG]:
            raise AttributeError('strand must be positive or negative to draw the transcript')
        if (mapping is None and target_width is None) or (mapping is not None and target_width is not None):
            raise AttributeError('mapping and target_width arguments are required and mutually exclusive')

        if mapping is None:
            target_width -= self.PADDING + self.LABEL_FONT_SIZE * 2 * self.FONT_WIDTH_HEIGHT_RATIO
            mapping = self._generate_interval_mapping(
                transcript.exons,
                target_width,
                self.EXON_INTRON_RATIO,
                self.EXON_MIN_WIDTH
            )

        main_group = canvas.g()

        y = self.BREAKPOINT_TOP_MARGIN

        if len(transcript.translations) == 0:
            y += self.TRACK_HEIGHT / 2
            exon_track_group = self._draw_exon_track(canvas, transcript, mapping, colors, labels)
            main_group.add(exon_track_group)
            y += self.TRACK_HEIGHT / 2
        else:
            # draw the protein features if there are any
            for tl in transcript.translations:
                # if the splicing takes up more room than the track we need to adjust for it
                y += max(self.SPLICE_HEIGHT, self.TRACK_HEIGHT / 2) - self.TRACK_HEIGHT / 2

                exon_track_group = self._draw_exon_track(canvas, transcript, mapping, colors, labels)
                exon_track_group.translate(0, y)

                y += self.TRACK_HEIGHT / 2
                # draw the splicing pattern
                polyline = []
                print('splicing_pattern', tl.splicing_pattern)
                for i in range(0, len(tl.splicing_pattern), 2):
                    a = tl.splicing_pattern[i]
                    b = tl.splicing_pattern[i + 1]
                    polyline.extend(
                        [(a, y), (a + (b - a) / 2, y - self.SPLICE_HEIGHT), (b, y)])

                p = canvas.polyline(polyline, fill='none', class_='splicing')
                p.dasharray(self.SPLICE_STROKE_DASHARRAY)
                p.stroke(self.SPLICE_COLOR, width=self.SPLICE_STROKE_WIDTH)
                main_group.add(p)
                main_group.add(exon_track_group)
                y += self.TRACK_HEIGHT / 2

                gp = canvas.g(class_='protein')
                # translation track
                y += self.PADDING
                # convert the AA position to cdna position, then convert the cdna to genomic, etc
                print(tl, tl.start, tl.end)
                print(transcript.convert_cdna_to_genomic(tl.start), transcript.convert_cdna_to_genomic(tl.end))
                s = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(tl.start))
                t = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(tl.end))
                print(s, t)
                reverse = False
                if s > t:
                    t, s = (s, t)
                    reverse = True
                gt = canvas.g(class_='translation')
                gp.add(gt)
                h = max(self.TRANSLATION_SCAFFOLD_HEIGHT, self.TRANSLATION_FONT_SIZE)
                gt.add(canvas.rect(
                    (s, h / 2 - self.TRANSLATION_SCAFFOLD_HEIGHT / 2),
                    (t - s + 1, self.TRANSLATION_SCAFFOLD_HEIGHT),
                    fill=self.TRANSLATION_SCAFFOLD_COLOR,
                    class_='scaffold'
                    ))
                gt.add(canvas.text(
                    self.TRANSLATION_END_MARKER if reverse else self.TRANSLATION_START_MARKER,
                    insert=(s - 1, h / 2),
                    fill=self.LABEL_COLOR,
                    style=self.FONT_STYLE.format(font_size=self.TRANSLATION_FONT_SIZE, text_anchor='end'),
                    class_='label'
                ))
                gt.add(canvas.text(
                    self.TRANSLATION_START_MARKER if reverse else self.TRANSLATION_END_MARKER,
                    insert=(t + 1, h / 2),
                    fill=self.LABEL_COLOR,
                    style=self.FONT_STYLE.format(font_size=self.TRANSLATION_FONT_SIZE, text_anchor='start'),
                    class_='label'
                ))
                py = h
                # now draw the domain tracks
                # need to convert the domain AA positions to cds positions to genomic
                for i, d in enumerate(sorted(tl.domains, key=lambda x: x.name)):
                    py += self.PADDING
                    gd = canvas.g(class_='domain')
                    gd.add(canvas.rect(
                        (0, self.DOMAIN_TRACK_HEIGHT / 2),
                        (target_width, self.DOMAIN_SCAFFOLD_HEIGHT),
                        fill=self.DOMAIN_SCAFFOLD_COLOR,
                        class_='scaffold'
                    ))
                    for region in d.regions:
                        # convert the AA position to cdna position, then convert the cdna to genomic, etc
                        s = tl.convert_aa_to_cdna(region[0])
                        t = tl.convert_aa_to_cdna(region[1])
                        temp = s | t
                        s = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(temp.start))
                        t = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(temp.end))
                        if s > t:
                            t, s = (s, t)
                        gd.add(canvas.rect(
                            (s, 0), (t - s + 1, self.DOMAIN_TRACK_HEIGHT),
                            fill=self.DOMAIN_COLOR, class_='region'))
                    gd.translate(0, py)

                    gd.add(canvas.text(
                        labels.add(d, self.DOMAIN_LABEL_PREFIX),
                        insert=(0 - self.PADDING, self.DOMAIN_TRACK_HEIGHT / 2),
                        fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(self.DOMAIN_COLOR),
                        style=self.FONT_STYLE.format(font_size=self.DOMAIN_LABEL_FONT_SIZE, text_anchor='end'),
                        class_='label'
                    ))
                    gp.add(gd)
                    py += self.DOMAIN_TRACK_HEIGHT
                gp.translate(0, y)
                y += py
                main_group.add(gp)

        y += self.BREAKPOINT_BOTTOM_MARGIN
        # now overlay the breakpoints on top of everything
        for i, b in enumerate(breakpoints):
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

    def draw_genes(self, canvas, genes, target_width, breakpoints=None, colors=None, labels=None):
        """
        draws the genes given in order of their start position trying to minimize
        the number of tracks required to avoid overlap

        Args:
            canvas (svgwrite.Drawing): the main svgwrite object used to create new svg elements
            target_width (int): the target width of the diagram
            genes (iterable of Gene): the list of genes to draw
            breakpoints (iterable of Breakpoint): the breakpoints to overlay
            colors (Dict of Gene to str): dictionary of the colors assigned to each Gene as fill

        Return:
            svgwrite.container.Group: the group element for the diagram.
                Has the added parameters of labels, height, and mapping
        """
        # mutable default argument parameters
        breakpoints = [] if breakpoints is None else breakpoints
        colors = {} if colors is None else colors
        labels = LabelMapping() if labels is None else labels

        main_group = canvas.g()
        mapping = self._generate_interval_mapping(
            [g for g in genes] + breakpoints,
            target_width,
            self.GENE_INTERGENIC_RATIO,
            self.GENE_MIN_WIDTH,
            buffer=self.GENE_MIN_BUFFER
        )
        gene_px_intervals = {}
        for i, gene in enumerate(sorted(genes, key=lambda x: x.start)):
            s = Interval.convert_pos(mapping, gene.start)
            t = Interval.convert_pos(mapping, gene.end)
            gene_px_intervals[Interval(s, t)] = gene
            labels.add(gene, self.GENE_LABEL_PREFIX)
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
                    canvas, gene, len(genepx),
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

    def read_config(self):
        pass

    def draw_breakpoint(self, canvas, breakpoint, width, height, label=''):
        """
        Args:
            canvas (svgwrite.Drawing): the main svgwrite object used to create new svg elements
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
            insert=(width / 2, y),
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
            canvas (svgwrite.Drawing): the main svgwrite object used to create new svg elements
            exon (Exon): the exon to draw
            width (int): the pixel width
            height (int): the pixel height
            fill (str): the fill color to use for the exon

        Return:
            svgwrite.container.Group: the group element for the diagram
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
                insert=(width / 2, height / 2),
                fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(fill),
                style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                class_='label'
            )
            g.add(t)
        return g

    def draw_gene(self, canvas, gene, width, height, fill, label=''):
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
            canvas (svgwrite.Drawing): the main svgwrite object used to create new svg elements
            gene (Gene): the gene to draw
            width (int): the pixel width
            height (int): the pixel height
            fill (str): the fill color to use for the gene

        Return:
            svgwrite.container.Group: the group element for the diagram
        """

        group = canvas.g(class_='gene')

        wrect = width - self.GENE_ARROW_WIDTH
        if wrect < 1:
            raise AttributeError('width is not sufficient to draw gene')

        label_color = self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(fill)

        if gene.strand == STRAND.POS:
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
                    insert=(wrect / 2, height / 2),
                    fill=label_color,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                    class_='label'
                ))
        elif gene.strand == STRAND.NEG:
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
                    insert=(wrect / 2 + self.GENE_ARROW_WIDTH, height / 2),
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
                    insert=(width / 2, height / 2),
                    fill=label_color,
                    style=self.FONT_STYLE.format(font_size=self.LABEL_FONT_SIZE, text_anchor='middle'),
                    class_='label'
                ))
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
    def _generate_interval_mapping(cls, input_intervals, target_width, ratio, min_width, buffer=0):
        intervals = []

        for i in Interval.min_nonoverlapping(*input_intervals):
            if len(intervals) == 0 or abs(Interval.dist(intervals[-1], i)) > 1:
                intervals.append(i)
            else:
                intervals[-1] = intervals[-1] | i

        start = max(intervals[0].start - buffer, 1)
        end = intervals[-1].end + buffer
        genic_length = sum([len(i) for i in intervals])
        intergenic_length = (end - start) + 1 - genic_length

        width = target_width - (len(intervals) * 2 + 1) * min_width  # reserved width

        if width < 0:
            raise AttributeError('width cannot accommodate the number of expected objects')

        intergenic_unit = width / (genic_length * ratio + intergenic_length)
        genic_unit = intergenic_unit * ratio

        if (intergenic_unit * intergenic_length + genic_unit * genic_length) > target_width:
            raise AssertionError('assert failed', intergenic_unit * intergenic_length + genic_unit * genic_length, target_width)

        mapping = []

        pos = 1
        # do the intergenic region prior to the first genic region
        if buffer > 0:
            ifrom = Interval(start, intervals[0].start - 1)
            ito = Interval(pos, pos + min_width + max(len(ifrom) * intergenic_unit - 1, 0))
            mapping.append((ifrom, ito))
            pos += len(ito)

        for i, curr in enumerate(intervals):
            if i > 0:
                prev = intervals[i - 1]
                ifrom = Interval(prev.end + 1, curr.start - 1)
                ito = Interval(pos, pos + min_width + max(len(ifrom) * intergenic_unit - 1, 0))
                mapping.append((ifrom, ito))
                pos += len(ito)

            ito = Interval(pos, pos + min_width + max(len(curr) * genic_unit - 1, 0))
            mapping.append((curr, ito))
            pos += len(ito)
        # now the last intergenic region will make up for the rounding error
        if buffer > 0:
            ifrom = Interval(intervals[-1].end + 1, end)
            ito = Interval(pos, pos + min_width + max(len(ifrom) * intergenic_unit - 1, 0))
            mapping.append((ifrom, ito))
            pos += len(ito)
        mapping[-1][1].end = target_width
        temp = mapping
        mapping = dict()
        for ifrom, ito in temp:
            mapping[ifrom] = ito
        return mapping
