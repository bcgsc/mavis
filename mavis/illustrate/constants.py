from colour import Color
from ..constants import GIESMA_STAIN


class DiagramSettings:
    def __init__(
        self, WIDTH=1000,
        SCAFFOLD_COLOR = '#000000',
        GENE1_COLOR_SELECTED='#518DC5',
        GENE2_COLOR_SELECTED = '#4C9677',
        GENE1_COLOR = '#657E91',
        GENE2_COLOR = '#325556',
        LABEL_COLOR = '#000000',
        DOMAIN_COLOR = '#ccccb3',
        DOMAIN_MISMATCH_COLOR = '#B2182B',
        SPLICE_COLOR = '#000000',
        BREAKPOINT_COLOR = '#000000',
        MASK_FILL = '#FFFFFF',
        MASK_OPACITY=0.7,
        DOMAIN_NAME_REGEX_FILTER='.*'
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
        self.SCAFFOLD_COLOR = SCAFFOLD_COLOR
        self.TRACK_HEIGHT = 50
        self.WIDTH = WIDTH
        # removing unsupported attr: 'alignment-baseline:central;dominant-baseline:central;' \
        self.FONT_STYLE = 'font-size:{font_size}px;font-weight:bold;alignment-baseline:baseline;' \
            'text-anchor:{text_anchor};font-family: consolas, courier new, monospace'
        # ratio for courier new which is wider than consolas, used for estimating width
        self.FONT_WIDTH_HEIGHT_RATIO = 1229 / 2048
        self.FONT_CENTRAL_SHIFT_RATIO = 0.3
        self.ABS_MIN_WIDTH = 0.01

        self.GENE1_COLOR_SELECTED = GENE1_COLOR_SELECTED
        self.GENE2_COLOR_SELECTED = GENE2_COLOR_SELECTED
        self.GENE1_COLOR = GENE1_COLOR
        self.GENE2_COLOR = GENE2_COLOR
        self.GENE_DEFAULT_COLOR = self.GENE1_COLOR
        self.GENE_MIN_BUFFER = 1000
        self.GENE_ARROW_WIDTH = 20
        self.GENE_INTERGENIC_RATIO = 5
        self.GENE_MIN_WIDTH = 40 + self.GENE_ARROW_WIDTH
        self.GENE_LABEL_PREFIX = 'G'

        self.LABEL_COLOR = LABEL_COLOR
        self.LABEL_FONT_SIZE = 20
        self.DYNAMIC_LABELS = True
        self.LABEL_LEFT_MARGIN = self.LABEL_FONT_SIZE * self.FONT_WIDTH_HEIGHT_RATIO * 4

        self.DOMAIN_COLOR = DOMAIN_COLOR
        self.DOMAIN_TRACK_HEIGHT = 30
        self.DOMAIN_SCAFFOLD_HEIGHT = 1
        self.DOMAIN_SCAFFOLD_COLOR = SCAFFOLD_COLOR
        self.DOMAIN_LABEL_PREFIX = 'D'
        self.DOMAIN_LABEL_FONT_SIZE = 20
        self.DOMAIN_MISMATCH_COLOR = DOMAIN_MISMATCH_COLOR
        self.DOMAIN_FILL_GRADIENT = [
            c.hex for c in Color(self.DOMAIN_MISMATCH_COLOR).range_to(Color(self.DOMAIN_COLOR), 10)]
        self.DOMAIN_NAME_REGEX_FILTER = DOMAIN_NAME_REGEX_FILTER
        self.PFAM_DOMAIN = '^PF\d+$'
        self.PFAM_LINK = 'http://pfam.xfam.org/family/{.name}'

        self.SPLICE_HEIGHT = self.TRACK_HEIGHT / 2
        self.SPLICE_STROKE_DASHARRAY = [2, 2]
        self.SPLICE_STROKE_WIDTH = 2
        self.SPLICE_COLOR = SPLICE_COLOR

        self.BREAKPOINT_STROKE_DASHARRAY = [3, 3]
        self.BREAKPOINT_ORIENT_STROKE_WIDTH = 2
        self.BREAKPOINT_COLOR = BREAKPOINT_COLOR
        self.BREAKPOINT_LABEL_FONT_SIZE = 20
        self.BREAKPOINT_BOTTOM_MARGIN = 20
        self.BREAKPOINT_TOP_MARGIN = self.PADDING * 2 + self.BREAKPOINT_LABEL_FONT_SIZE + self.BREAKPOINT_BOTTOM_MARGIN
        self.BREAKPOINT_LABEL_PREFIX = 'B'

        self.MARKER_LABEL_FONT_SIZE = self.BREAKPOINT_LABEL_FONT_SIZE
        self.MARKER_LABEL_PREFIX = 'M'
        self.MARKER_TOP_MARGIN = self.BREAKPOINT_TOP_MARGIN
        self.MARKER_BOTTOM_MARGIN = self.BREAKPOINT_BOTTOM_MARGIN
        self.MARKER_COLOR = self.BREAKPOINT_COLOR

        self.TRANSCRIPT_HYPERLINK = 'http://dec2013.archive.ensembl.org/Homo_sapiens/Transcript/Summary?' \
            'db=core;t={.name}'
        self.EXON_FONT_SIZE = 20
        self.EXON_TEAR_TOOTH_WIDTH = 2
        self.EXON_MIN_WIDTH = max([
            self.MIN_WIDTH + self.EXON_TEAR_TOOTH_WIDTH * 2,
            self.EXON_FONT_SIZE * 2 * self.FONT_WIDTH_HEIGHT_RATIO
        ])
        self.EXON_TEAR_TOOTH_HEIGHT = 2
        self.EXON_INTRON_RATIO = 20
        self.EXON1_COLOR = self.GENE1_COLOR_SELECTED
        self.EXON2_COLOR = self.GENE2_COLOR_SELECTED

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
        self.LEGEND_SWATCH_STROKE = '#000000'
        self.LEGEND_FONT_COLOR = '#000000'
        self.LEGEND_BORDER_STROKE = '#000000'
        self.LEGEND_BORDER_STROKE_WIDTH = 1

        self.TEMPLATE_BAND_STROKE_WIDTH = 0.5
        temp = [c.hex for c in Color('#FFFFFF').range_to(Color('#000000'), 5)]
        self.TEMPLATE_BAND_FILL = {
            GIESMA_STAIN.ACEN: '#800000',
            GIESMA_STAIN.GPOS25: temp[1],
            GIESMA_STAIN.GPOS50: temp[2],
            GIESMA_STAIN.GPOS75: temp[3],
            GIESMA_STAIN.GPOS100: temp[4],
            GIESMA_STAIN.GNEG: '#FFFFFF'
        }
        self.TEMPLATE_BAND_STROKE = '#000000'
        self.TEMPLATE_TRACK_HEIGHT = max([
            self.TRACK_HEIGHT / 3,
            self.LABEL_FONT_SIZE - self.BREAKPOINT_BOTTOM_MARGIN -
            self.BREAKPOINT_TOP_MARGIN + self.BREAKPOINT_LABEL_FONT_SIZE])
        self.TEMPLATE_DEFAULT_FILL = '#FFFFFF'
        self.TEMPLATE_BAND_MIN_WIDTH = 2
        self.TEMPLATE_LABEL_PREFIX = 'C'

        self.REGION_LABEL_PREFIX = 'R'
        self.OVERLAY_LEFT_LABEL = 16 * self.FONT_WIDTH_HEIGHT_RATIO * self.EXON_FONT_SIZE

        self.SCATTER_AXIS_FONT_SIZE = 12
        self.SCATTER_ERROR_BAR_STROKE_WIDTH = 1
        self.SCATTER_MARKER_RADIUS = 2
        self.SCATTER_YAXIS_TICK_SIZE = self.PADDING
        self.SCATTER_YTICK_FONT_SIZE = 10

        self.MASK_FILL = MASK_FILL
        self.MASK_OPACITY = MASK_OPACITY
