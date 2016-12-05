from structural_variant.annotate import Gene
from svgwrite import Drawing
from structural_variant.interval import Interval
from structural_variant.constants import STRAND

# draw gene level view
# draw gene box

GENE_FORWARD_SVG="""
<g class='.gene'>
    <path d='M{0} 0 L{1} {2}'>
    </path>
    <rect>
    </rect>
</g>"""


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
        LABEL_FILL='#FFFFFF',
        LABEL_STROKE='#000000',
        TRANSCRIPT1_COLOR='#4C9677',
        TRANSCRIPT2_COLOR='#518DC5',
        TRANSCRIPT1_UTR_COLOR='#7bddc1',
        TRANSCRIPT2_UTR_COLOR='#7dc3d8',
        DOMAIN_COLOR='#b8d3ba',
        WIDTH=1000,
        PADDING=5,
        GENE_MIN_BUFFER=200
    ):
        self.GENE1_COLOR_SELECTED = GENE1_COLOR_SELECTED
        self.GENE2_COLOR_SELECTED = GENE2_COLOR_SELECTED
        self.GENE1_COLOR = GENE1_COLOR
        self.GENE2_COLOR = GENE2_COLOR
        self.LABEL_FILL = LABEL_FILL
        self.LABEL_STROKE = LABEL_STROKE
        self.TRANSCRIPT1_COLOR = TRANSCRIPT1_COLOR
        self.TRANSCRIPT2_COLOR = TRANSCRIPT2_COLOR
        self.TRANSCRIPT1_UTR_COLOR = TRANSCRIPT1_UTR_COLOR
        self.TRANSCRIPT2_UTR_COLOR = TRANSCRIPT2_UTR_COLOR
        self.DOMAIN_COLOR = DOMAIN_COLOR
        self.TRACK_HEIGHT = 50
        self.LINE_WIDTH = 4
        self.LINE_COLOR = '#000000'
        self.DOMAIN_TRACK_HEIGHT = 20
        self.WIDTH = WIDTH
        self.GENE_INTERGENIC_RATIO = 10
        self.EXON_INTRON_RATIO = 10
        self.MIN_WIDTH = 2 # no element (exon, gene, etc can be less than this wide)
        self.GENE_MIN_BUFFER = GENE_MIN_BUFFER
        self.TRACK_LINE_HEIGHT = 4
        self.LEFT_MARGIN = 20
        self.RIGHT_MARGIN = 20
        self.TOP_MARGIN = 20
        self.BOTTOM_MARGIN = 20
        self.INNER_MARGIN = 20
        self.GENE_ARROW_WIDTH = 20
        self.PADDING = PADDING

    def draw_legend(self):
        pass

    def draw(self, ann):
        # decide one or two gene tracks
        left_gene_track = None
        right_gene_track = None

        canvas = Drawing(height=100, width=self.WIDTH) # just set the height for now and change later

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

    def draw_transcript(self, canvas, target_width, t, exon_color, utr_color):
        main_group = canvas.g()
        mapping = self._generate_exon_mapping(target_width, t.exons)

        main_group.add(
            canvas.rect(
                (0, self.TRACK_HEIGHT / 2 - self.LINE_WIDTH / 2),
                (target_width, self.LINE_WIDTH),
                fill=self.LINE_COLOR
                ))
        y = 0
        # draw the exons
        for exon in sorted(t.exons, key=lambda x: x.start):
            s = Interval.convert_pos(mapping, exon.start)
            t = Interval.convert_pos(mapping, exon.end)

            group = canvas.g(class_='exon')
            main_group.add(group)

            group.add(
                canvas.rect(
                    (s, y),
                    (t - s + 1, self.TRACK_HEIGHT),
                    fill=exon_color
                ))
        #y += self.TRACK_HEIGHT + self.PADDING

        # now draw the domain tracks
        # need to convert the domain AA positions to cds positions to genomic

        return main_group


    def draw_gene_subdiagram(self, canvas, target_width, genes, breakpoint=None, colors={}):
        main_group = canvas.g()
        mapping = self._generate_gene_mapping(target_width, genes)
        print('draw_gene_subdiagram', canvas, target_width, genes, colors)
        gene_px_intervals = {}
        for gene in genes:
            s = Interval.convert_pos(mapping, gene.start)
            t = Interval.convert_pos(mapping, gene.end)
            gene_px_intervals[Interval(s, t)] = gene
        tracks = Diagram._split_intervals_into_tracks(gene_px_intervals)

        y = 0

        main_group.add(
            canvas.rect(
                (0, self.TRACK_HEIGHT / 2 - self.LINE_WIDTH / 2 + (len(tracks) - 1) * (self.TRACK_HEIGHT + self.PADDING)),
                (target_width, self.LINE_WIDTH),
                fill=self.LINE_COLOR
                ))
        tracks.reverse()
        for track in tracks:  # svg works from top down
            for genepx in track:
                # draw the gene
                wrect = genepx.end - self.GENE_ARROW_WIDTH - genepx.start + 1
                gene = gene_px_intervals[genepx]
                group = canvas.g(class_='gene')
                main_group.add(group)

                if gene.strand == STRAND.POS:
                    group.add(
                        canvas.rect(
                            (genepx.start, y),
                            (wrect, self.TRACK_HEIGHT),
                            fill=colors.get(gene, self.GENE1_COLOR)
                        )
                    )
                    group.add(
                        canvas.polyline(
                            [
                                (genepx.start + wrect, y),
                                (genepx.start + wrect + self.GENE_ARROW_WIDTH, y + self.TRACK_HEIGHT / 2),
                                (genepx.start + wrect, y + self.TRACK_HEIGHT)
                            ],
                            fill=colors.get(gene, self.GENE1_COLOR)
                        )
                    )
                elif gene.strand == STRAND.NEG:
                    group.add(
                        canvas.rect(
                            (genepx.start + self.GENE_ARROW_WIDTH, y),
                            (wrect, self.TRACK_HEIGHT),
                            fill=colors.get(gene, self.GENE1_COLOR)
                        )
                    )
                    group.add(
                        canvas.polyline(
                            [
                                (genepx.start + self.GENE_ARROW_WIDTH, y),
                                (genepx.start, y + self.TRACK_HEIGHT / 2),
                                (genepx.start + self.GENE_ARROW_WIDTH, y + self.TRACK_HEIGHT)
                            ],
                            fill=colors.get(gene, self.GENE1_COLOR)
                        )
                    )
                else:
                    raise AttributeError('gene must specify positive or negative strand to be drawn')
            y += self.TRACK_HEIGHT + self.PADDING
        return main_group

    def read_config(self):
        pass

    def _generate_exon_mapping(self, target_width, exons):
        return self._generate_interval_mapping(
            target_width,
            exons,
            self.EXON_INTRON_RATIO,
            self.MIN_WIDTH
        )

    def _generate_gene_mapping(self, target_width, genes):
        return self._generate_interval_mapping(
            target_width,
            genes,
            self.GENE_INTERGENIC_RATIO,
            self.MIN_WIDTH + self.GENE_ARROW_WIDTH,
            buffer=self.GENE_MIN_BUFFER
        )

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
