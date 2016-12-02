from structural_variant.annotate import Gene
import svgwrite
from structural_variant.interval import Interval

# draw gene level view
# draw gene box


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
        GENE_LABEL_FILL='#FFFFFF',
        GENE_LABEL_STROKE='#000000',
        TRANSCRIPT1_COLOR='#4C9677',
        TRANSCRIPT2_COLOR='#518DC5',
        TRANSCRIPT1_UTR_COLOR='#7bddc1',
        TRANSCRIPT2_UTR_COLOR='#7dc3d8',
        DOMAIN_COLOR='#b8d3ba',
        WIDTH=1000,
    ):
        self.GENE1_COLOR_SELECTED = GENE1_COLOR_SELECTED
        self.GENE2_COLOR_SELECTED = GENE2_COLOR_SELECTED
        self.GENE1_COLOR = GENE1_COLOR
        self.GENE2_COLOR = GENE2_COLOR
        self.GENE_LABEL_FILL = GENE_LABEL_FILL
        self.GENE_LABEL_STROKE = GENE_LABEL_STROKE
        self.TRANSCRIPT1_COLOR = TRANSCRIPT1_COLOR
        self.TRANSCRIPT2_COLOR = TRANSCRIPT2_COLOR
        self.TRANSCRIPT1_UTR_COLOR = TRANSCRIPT1_UTR_COLOR
        self.TRANSCRIPT2_UTR_COLOR = TRANSCRIPT2_UTR_COLOR
        self.DOMAIN_COLOR = DOMAIN_COLOR
        self.TRACK_HEIGHT = 50
        self.DOMAIN_TRACK_HEIGHT = 20
        self.WIDTH = WIDTH
        self.GENE_INTERGENIC_RATIO = 10
        self.EXON_INTRON_RATIO = 10
        self.MIN_WIDTH = 2 # no element (exon, gene, etc can be less than this wide)
        self.MIN_INTERGENIC_BUFFER = 200
        self.TRACK_LINE_HEIGHT = 4
        self.LEFT_MARGIN = 20
        self.RIGHT_MARGIN = 20
        self.TOP_MARGIN = 20
        self.BOTTOM_MARGIN = 20
        self.INNER_MARGIN = 20
    
    def draw_legend(self):
        pass

    def draw(self, ann):
        # decide one or two gene tracks
        left_gene_track = None
        right_gene_track = None

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
            
            gene_track, mapping = self._draw_gene_track(w, genes)

        # add transcript tracks if applicable (and domains)
        # add fusion track if applicable


    def read_config(self):
        pass
   
    def _generate_gene_mapping(self, target_width, genes):
        intervals = []
        
        for i in Interval.min_nonoverlapping(*genes):
            if len(intervals) == 0 or abs(Interval.dist(intervals[-1], i)) > 1:
                intervals.append(i)
            else:
                intervals[-1] = intervals[-1] | i
        print(intervals)
        start = intervals[0].start - self.MIN_INTERGENIC_BUFFER
        end = intervals[1].end + self.MIN_INTERGENIC_BUFFER
        genic_length = sum([len(i) for i in intervals])
        intergenic_length = end - start + 1 - genic_length
        
        width = target_width - (len(intervals) * 2 + 1) * self.MIN_WIDTH  # reserved width

        if width < 0:
            raise AttributeError('width cannot accommodate the number of expected objects')

        intergenic_unit = width / (genic_length * self.GENE_INTERGENIC_RATIO + intergenic_length)
        genic_unit = intergenic_unit * self.GENE_INTERGENIC_RATIO

        mapping = {}
        
        pos = 1
        # do the intergenic region prior to the first genic region
        ifrom = Interval(start, intervals[0].start - 1)
        ito = Interval(pos, pos + self.MIN_WIDTH - 1 + max(len(ifrom) * intergenic_unit - 1, 0))
        mapping[ifrom] = ito
        pos += len(ito)

        for i, curr in enumerate(intervals):
            print('curr', curr)
            if i > 0:
                prev = intervals[i - 1]
                print('prev', prev)
                ifrom = Interval(prev.end + 1, curr.start - 1)
                ito = Interval(pos, pos + self.MIN_WIDTH - 1 + max(len(ifrom) * intergenic_unit - 1, 0))
                mapping[ifrom] = ito
                pos += len(ito)
            
            ito = Interval(pos, pos + self.MIN_WIDTH - 1 + max(len(curr) * genic_unit - 1, 0))
            mapping[curr] = ito
            pos += len(ito)
        
        # now the last intergenic region will make up for the rounding error
        ifrom = Interval(intervals[-1].end + 1, end)
        ito = Interval(pos, target_width - 1)
        mapping[ifrom] = ito
        pos += len(ito)
        print(sorted(mapping.items()))
        print('pos', pos)
        print('target_width', target_width)
        self.assertTrue(False)
        return mapping


