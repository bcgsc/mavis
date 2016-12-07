from svgwrite import Drawing
from structural_variant.interval import Interval
from string import ascii_uppercase
from structural_variant.constants import STRAND, ORIENT
from colour import Color

# draw gene level view
# draw gene box
HEX_WHITE = '#FFFFFF'
HEX_BLACK = '#000000'

class Diagram:
    """
    class which holds the settings for drawing a fusion digram
    """
    def __init__(
        self,
        GENE1_COLOR_SELECTED='#4C9677',
        GENE2_COLOR_SELECTED='#518DC5',
        GENE1_COLOR='#325556',
        GENE2_COLOR='#657E91',
        LABEL_COLOR='#000000',
        LABEL_FONT_SIZE=28,
        EXON1_COLOR='#4C9677',
        EXON2_COLOR='#518DC5',
        EXON1_UTR_COLOR='#7bddc1',
        EXON2_UTR_COLOR='#7dc3d8',
        DOMAIN_COLOR='#b8d3ba',
        WIDTH=1000,
        PADDING=5,
        DYNAMIC_LABELS=True
    ):
        self.MIN_WIDTH = 2  # no element (exon, gene, etc can be less than this wide)
        self.TRACK_LINE_HEIGHT = 4
        self.LEFT_MARGIN = 20
        self.RIGHT_MARGIN = 20
        self.TOP_MARGIN = 20
        self.BOTTOM_MARGIN = 20
        self.INNER_MARGIN = 20
        self.PADDING = PADDING
        self.LINE_WIDTH = 3
        self.LINE_COLOR = '#000000'
        self.TRACK_HEIGHT = 50
        self.WIDTH = WIDTH

        self.GENE1_COLOR_SELECTED = GENE1_COLOR_SELECTED
        self.GENE2_COLOR_SELECTED = GENE2_COLOR_SELECTED
        self.GENE1_COLOR = GENE1_COLOR
        self.GENE2_COLOR = GENE2_COLOR
        self.GENE_MIN_BUFFER = 200
        self.GENE_ARROW_WIDTH = 20
        self.GENE_INTERGENIC_RATIO = 5
        self.GENE_MIN_WIDTH = self.MIN_WIDTH + self.GENE_ARROW_WIDTH

        self.LABEL_COLOR = LABEL_COLOR
        self.LABEL_FONT_SIZE = LABEL_FONT_SIZE
        self.DYNAMIC_LABELS = DYNAMIC_LABELS

        self.DOMAIN_COLOR = DOMAIN_COLOR
        self.DOMAIN_TRACK_HEIGHT = 30

        self.SPLICE_HEIGHT = self.TRACK_HEIGHT
        self.SPLICE_STROKE_DASHARRAY = [2, 2]
        self.SPLICE_STROKE_WIDTH = 2

        self.BREAKPOINT_STROKE_DASHARRAY = [3, 3]
        self.BREAKPOINT_ORIENT_STROKE_WIDTH = 2

        self.EXON_TEAR_TOOTH_WIDTH = 2
        self.EXON_MIN_WIDTH = self.MIN_WIDTH + self.EXON_TEAR_TOOTH_WIDTH * 2
        self.EXON_TEAR_TOOTH_HEIGHT = 2
        self.EXON_INTRON_RATIO = 5
        self.EXON1_COLOR = EXON1_COLOR
        self.EXON2_COLOR = EXON2_COLOR
        self.EXON1_UTR_COLOR = EXON1_UTR_COLOR
        self.EXON2_UTR_COLOR = EXON2_UTR_COLOR

    def draw_legend(self):
        pass

    def draw(self, ann):
        # decide one or two gene tracks
        left_gene_track = None
        right_gene_track = None

        canvas = Drawing(height=100, width=self.WIDTH)  # just set the height for now and change later

        if ann.interchromosomal:  # two gene tracks
            pass
            w = self.WIDTH - self.LEFT_MARGIN - self.RIGHT_MARGIN - self.INNER_MARGIN
            w = w // 2
        else:  # single gene track
            w = self.WIDTH - self.LEFT_MARGIN - self.RIGHT_MARGIN
            genes = ann.genes_at_break1 | ann.genes_at_break2 | ann.encompassed_genes
            if ann.transcript1 and hasattr(ann.transcript1, 'gene'):
                genes.add(ann.transcript1.gene)
            if ann.transcript2 and hasattr(ann.transcript2, 'gene'):
                genes.add(ann.transcript2.gene)


            # draw the genes in order


        # add transcript tracks if applicable (and domains)
        # add fusion track if applicable

    def draw_transcript(self, canvas, target_width, transcript, exon_color, utr_color, abrogated_splice_sites=[], breakpoints=[]):
        """
        builds an svg group representing the transcript. Exons are drawn in a track with the splicing
        information and domains are drawn in separate tracks below

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
        """
        labels = {}
        if transcript.strand not in [STRAND.POS, STRAND.NEG]:
            raise AttributeError('strand must be positive or negative to draw the transcript')
        exons = sorted(transcript.exons, key=lambda x: x.start)
        main_group = canvas.g()
        x = self.PADDING * 2 + self.LABEL_FONT_SIZE * 2
        target_width -= x

        mapping = self._generate_interval_mapping(
            target_width,
            transcript.exons,
            self.EXON_INTRON_RATIO,
            self.EXON_MIN_WIDTH
        )


        y = max(self.SPLICE_HEIGHT, self.TRACK_HEIGHT / 2) + self.PADDING

        # draw the splicing lines
        splice_sites = []
        for i, exon in enumerate(exons):
            if i > 0 and exon.start not in abrogated_splice_sites:
                splice_sites.append(Interval.convert_pos(mapping, exon.start))
            if i < len(transcript.exons) - 1 and exon.end not in abrogated_splice_sites:
                splice_sites.append(Interval.convert_pos(mapping, exon.end))

        if len(splice_sites) % 2 != 0:
            raise AttributeError('splice sites must be an even multiple')

        polyline = []

        for i in range(0, len(splice_sites), 2):
            a = splice_sites[i]
            b = splice_sites[i + 1]
            polyline.extend(
                [(x + a, y), (x + a + (b - a) / 2, y - self.SPLICE_HEIGHT), (x + b, y)])

        p = canvas.polyline(polyline, fill='none', class_='splicing')
        p.dasharray(self.SPLICE_STROKE_DASHARRAY)
        p.stroke(self.LINE_COLOR, width=self.SPLICE_STROKE_WIDTH)
        main_group.add(p)

        main_group.add(
            canvas.rect(
                (x, y - self.LINE_WIDTH / 2),
                (target_width, self.LINE_WIDTH),
                fill=self.LINE_COLOR,
                class_='scaffold'
                ))

        # draw the exons
        for i, exon in enumerate(exons):
            s = Interval.convert_pos(mapping, exon.start)
            t = Interval.convert_pos(mapping, exon.end)
            pxi = Interval(s, t)
            utr = []
            for u in transcript.genomic_utr_regions():
                temp = Interval(Interval.convert_pos(mapping, u.start), Interval.convert_pos(mapping, u.end))
                if Interval.overlaps(u, exon):
                    temp = temp & pxi
                    utr.append(Interval(temp.start - s, temp.end - s))
            group = self.draw_exon(
                canvas, exon, len(pxi), self.TRACK_HEIGHT, exon_color,
                utr=utr,
                utr_fill=utr_color,
                label=i + 1 if transcript.strand == STRAND.POS else len(exons) - i)
            group.translate(x + pxi.start, y - self.TRACK_HEIGHT / 2)
            main_group.add(group)

        y += self.TRACK_HEIGHT / 2 + self.PADDING

        # now draw the domain tracks
        # need to convert the domain AA positions to cds positions to genomic
        for i, d in enumerate(sorted(transcript.domains, key=lambda x: x.name)):
            g = canvas.g(class_='domain')
            g.add(canvas.rect(
                (x, self.DOMAIN_TRACK_HEIGHT / 2),
                (target_width, self.LINE_WIDTH),
                fill=self.LINE_COLOR,
                class_='scaffold'
                ))
            for region in d.regions:
                # convert the AA position to cdna position, then convert the cdna to genomic, etc
                s = transcript.convert_aa_to_cdna(region[0])
                t = transcript.convert_aa_to_cdna(region[1])
                temp = s | t
                s = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(temp.start))
                t = Interval.convert_pos(mapping, transcript.convert_cdna_to_genomic(temp.end))
                if s > t:
                    t, s = (s, t)
                g.add(canvas.rect(
                    (s, 0), (t - s + 1, self.DOMAIN_TRACK_HEIGHT),
                    fill=self.DOMAIN_COLOR, class_='region'))
            g.translate(0, y)
            g.add(canvas.text(
                'D{}'.format(i + 1),
                insert=(self.PADDING, self.DOMAIN_TRACK_HEIGHT / 2),
                fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(self.DOMAIN_COLOR),
                font_size=self.LABEL_FONT_SIZE,
                alignment_baseline='central',
                class_='label'
            ))
            main_group.add(g)

            y += self.DOMAIN_TRACK_HEIGHT + self.PADDING


        # now overlay the breakpoints on top of everything
        for b in breakpoints:
            s = Interval.convert_pos(mapping, b.start)
            t = Interval.convert_pos(mapping, b.end)
            bg = self.draw_breakpoint(canvas, b, abs(t - s) + 1, y)
            bg.translate(x + s, 0)
            main_group.add(bg)

        setattr(main_group, 'height', y)
        setattr(main_group, 'mapping', mapping)
        setattr(main_group, 'labels', labels)
        return main_group

    def draw_genes(self, canvas, target_width, genes, breakpoints=[], colors={}):
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
            svgwrite.container.Group: the group element for the diagram
        """
        main_group = canvas.g()
        mapping = self._generate_interval_mapping(
            target_width,
            [g for g in genes] + breakpoints,
            self.GENE_INTERGENIC_RATIO,
            self.GENE_MIN_WIDTH,
            buffer=self.GENE_MIN_BUFFER
        )
        print('draw_gene_subdiagram', canvas, target_width, genes, colors)
        gene_px_intervals = {}
        labels = {}
        for i, gene in enumerate(sorted(genes, key=lambda x: x.start)):
            s = Interval.convert_pos(mapping, gene.start)
            t = Interval.convert_pos(mapping, gene.end)
            gene_px_intervals[Interval(s, t)] = gene
            labels[gene] = 'G{}'.format(i + 1)
        tracks = Diagram._split_intervals_into_tracks(gene_px_intervals)

        y = self.PADDING

        main_group.add(
            canvas.rect(
                (0, y + self.TRACK_HEIGHT / 2 - self.LINE_WIDTH / 2 + (len(tracks) - 1) * (self.TRACK_HEIGHT + self.PADDING)),
                (target_width, self.LINE_WIDTH),
                fill=self.LINE_COLOR,
                class_='scaffold'
                ))
        tracks.reverse()
        for track in tracks:  # svg works from top down
            for genepx in track:
                # draw the gene
                gene = gene_px_intervals[genepx]
                group = self.draw_gene(
                    canvas, gene, len(genepx), self.TRACK_HEIGHT, colors.get(gene, self.GENE1_COLOR), labels[gene])
                group.translate(genepx.start, y)
                main_group.add(group)
            y += self.TRACK_HEIGHT + self.PADDING
        # now overlay the breakpoints on top of everything
        for b in breakpoints:
            s = Interval.convert_pos(mapping, b.start)
            t = Interval.convert_pos(mapping, b.end)
            bg = self.draw_breakpoint(canvas, b, abs(t - s) + 1, y)
            bg.translate(s, 0)
            main_group.add(bg)

        temp = labels
        labels = {}
        for k, v in temp.items():
            labels[v] = k
        setattr(main_group, 'height', y)
        setattr(main_group, 'mapping', mapping)
        setattr(main_group, 'labels', labels)

        return main_group

    def read_config(self):
        pass

    def draw_breakpoint(self, canvas, breakpoint, width, height):
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
        r = canvas.rect(
            (0, 0),
            (width, height),
            stroke=self.LINE_COLOR,
            fill='none'
        )
        r.dasharray(self.BREAKPOINT_STROKE_DASHARRAY)
        g.add(r)

        if breakpoint.orient == ORIENT.LEFT:
            l = canvas.line((0, 0), (0, height))
            l.stroke(self.LINE_COLOR, width=self.BREAKPOINT_ORIENT_STROKE_WIDTH)
            g.add(l)
        elif breakpoint.orient == ORIENT.RIGHT:
            l = canvas.line((width, 0), (width, height))
            l.stroke(self.LINE_COLOR, width=self.BREAKPOINT_ORIENT_STROKE_WIDTH)
            g.add(l)
        return g

    def draw_exon(self, canvas, exon, width, height, fill, utr=[], utr_fill='#000000', label='', tear_left=False, tear_right=False):
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
            for u in utr:
                if u.start < 0 or u.end > width:
                    print(u, width, height)
                    raise AttributeError('utr outside exon region')
                g.add(canvas.rect((u.start, 0), (len(u), height), class_='UTR', fill=utr_fill))
            t = canvas.text(
                label,
                insert=(width / 2, height / 2),
                fill=self.LABEL_COLOR if not self.DYNAMIC_LABELS else Diagram.dynamic_label_color(fill),
                font_size=self.LABEL_FONT_SIZE,
                text_anchor='middle',
                alignment_baseline='central',
                class_='label'
            )
            g.add(t)
        return g

    def draw_gene(self, canvas, gene, width, height, fill, label=''):
        """
        generates the svg object representing a gene

        ::

            gene on the positive/forward strand

            +-----\
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
                    font_size=self.LABEL_FONT_SIZE,
                    text_anchor='middle',
                    alignment_baseline='central',
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
                    font_size=self.LABEL_FONT_SIZE,
                    text_anchor='middle',
                    alignment_baseline='central',
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
                    font_size=self.LABEL_FONT_SIZE,
                    text_anchor='middle',
                    alignment_baseline='central',
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
    def _generate_interval_mapping(cls, target_width, input_intervals, ratio, min_width, buffer=0):
        intervals = []

        for i in Interval.min_nonoverlapping(*input_intervals):
            if len(intervals) == 0 or abs(Interval.dist(intervals[-1], i)) > 1:
                intervals.append(i)
            else:
                intervals[-1] = intervals[-1] | i

        start = intervals[0].start - buffer
        end = intervals[1].end + buffer
        genic_length = sum([len(i) for i in intervals])
        intergenic_length = end - start + 1 - genic_length

        width = target_width - (len(intervals) * 2 + 1) * min_width  # reserved width

        if width < 0:
            raise AttributeError('width cannot accommodate the number of expected objects')

        intergenic_unit = width / (genic_length * ratio + intergenic_length)
        genic_unit = intergenic_unit * ratio

        mapping = []

        pos = 1
        # do the intergenic region prior to the first genic region
        if buffer > 0:
            ifrom = Interval(start, intervals[0].start - 1)
            ito = Interval(pos, pos + min_width - 1 + max(len(ifrom) * intergenic_unit - 1, 0))
            mapping.append((ifrom, ito))
            pos += len(ito)

        for i, curr in enumerate(intervals):
            if i > 0:
                prev = intervals[i - 1]
                ifrom = Interval(prev.end + 1, curr.start - 1)
                ito = Interval(pos, pos + min_width - 1 + max(len(ifrom) * intergenic_unit - 1, 0))
                mapping.append((ifrom, ito))
                pos += len(ito)

            ito = Interval(pos, pos + min_width - 1 + max(len(curr) * genic_unit - 1, 0))
            mapping.append((curr, ito))
            pos += len(ito)

        # now the last intergenic region will make up for the rounding error
        if buffer > 0:
            ifrom = Interval(intervals[-1].end + 1, end)
            ito = Interval(pos, pos + min_width - 1 + max(len(ifrom) * intergenic_unit - 1, 0))
            mapping.append((ifrom, ito))
            pos += len(ito)

        mapping[-1][1].end = target_width - 1

        temp = mapping
        mapping = dict()
        for ifrom, ito in temp:
            mapping[ifrom] = ito
        return mapping
