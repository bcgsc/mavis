from colour import Color
from ..constants import GIEMSA_STAIN, float_fraction
from ..util import WeakMavisNamespace

DEFAULTS = WeakMavisNamespace()
"""
- :term:`breakpoint_color`
- :term:`domain_color`
- :term:`domain_mismatch_color`
- :term:`domain_name_regex_filter`
- :term:`domain_scaffold_color`
- :term:`drawing_width_iter_increase`
- :term:`gene1_color_selected`
- :term:`gene1_color`
- :term:`gene2_color_selected`
- :term:`gene2_color`
- :term:`label_color`
- :term:`mask_fill`
- :term:`mask_opacity`
- :term:`max_drawing_retries`
- :term:`novel_exon_color`
- :term:`scaffold_color`
- :term:`splice_color`
- :term:`width`
"""
DEFAULTS.add(
    'width', 1000, defn='The drawing width in pixels')
DEFAULTS.add(
    'domain_name_regex_filter', r'^PF\d+$',
    defn='The regular expression used to select domains to be displayed (filtered by name)')
DEFAULTS.add(
    'max_drawing_retries', 5,
    defn='The maximum number of retries for attempting a drawing. Each iteration the width is extended. If it '
    'is still insufficient after this number a gene-level only drawing will be output')
DEFAULTS.add('scaffold_color', '#000000', defn='The color used for the gene/transcripts scaffolds')
DEFAULTS.add('gene1_color_selected', '#518dc5', defn='The color of the first gene')
DEFAULTS.add('gene2_color_selected', '#4c9677', defn='The color of the second gene')
DEFAULTS.add('gene1_color', '#657e91', defn='The color of genes near the first gene')
DEFAULTS.add('gene2_color', '#325556', defn='The color of genes near the second gene')
DEFAULTS.add('label_color', '#000000', defn='The label color')
DEFAULTS.add('domain_color', '#ccccb3', defn='Domain fill color')
DEFAULTS.add('domain_mismatch_color', '#b2182b', defn='Domain fill color on 0%% match')
DEFAULTS.add('novel_exon_color', '#5D3F6A', defn='Novel Exon fill color')
DEFAULTS.add('splice_color', '#000000', defn='Splicing lines color')
DEFAULTS.add('breakpoint_color', '#000000', defn='Breakpoint outline color')
DEFAULTS.add('mask_fill', '#ffffff', defn='Color of mask (for deleted region etc.)')
DEFAULTS.add('mask_opacity', 0.7, defn='opacity of the mask layer', cast_type=float_fraction)
DEFAULTS.add('domain_scaffold_color', '#000000', defn='The color of the domain scaffold')
DEFAULTS.add('drawing_width_iter_increase', 500, defn='The amount (in  pixels) by which to increase the drawing width upon failure to fit')
DEFAULTS.add('exon_min_focus_size', 10, defn='minimum size of an exon for it to be granted a label or min exon width')


class DiagramSettings:
    """
    holds settings related to colors/sizes for the drawing
    """
    def __init__(
        self, **kwargs
    ):
        inputs = {}
        inputs.update(DEFAULTS.items())
        inputs.update(kwargs)
        for arg, val in inputs.items():
            if arg not in DEFAULTS:
                raise KeyError('unrecognized argument', arg)
            setattr(self, arg, val)
        self.min_width = 10  # no element (exon, gene, etc can be less than this wide)
        self.track_line_height = 4
        self.left_margin = 20
        self.right_margin = 20
        self.top_margin = 20
        self.bottom_margin = 20
        self.inner_margin = 20
        self.padding = 5
        self.scaffold_height = 3
        self.track_height = 50

        self.ins_increase = 6
        # removing unsupported attr: 'alignment-baseline:central;dominant-baseline:central;' \
        self.font_style = 'font-size:{font_size}px;font-weight:bold;alignment-baseline:baseline;' \
            'text-anchor:{text_anchor};font-family: consolas, courier new, monospace'
        # ratio for courier new which is wider than consolas, used for estimating width
        self.font_width_height_ratio = 1229 / 2048
        self.font_central_shift_ratio = 0.3
        self.non_focus_min_width = 2

        self.gene_default_color = self.gene1_color
        self.gene_min_buffer = 1000
        self.gene_arrow_width = 20
        self.gene_intergenic_ratio = 5
        self.gene_min_width = 40 + self.gene_arrow_width
        self.gene_label_prefix = 'G'

        self.label_font_size = 20
        self.dynamic_labels = True
        self.label_left_margin = self.label_font_size * self.font_width_height_ratio * 4

        self.domain_track_height = 30
        self.domain_scaffold_height = 1
        self.domain_label_prefix = 'D'
        self.domain_label_font_size = 20
        self.domain_fill_gradient = [
            c.hex for c in Color(self.domain_mismatch_color).range_to(Color(self.domain_color), 10)]
        self.domain_links = {
            r'^PF\d+$': 'http://pfam.xfam.org/family/{.name}'
        }

        self.splice_height = self.track_height / 2
        self.splice_stroke_dasharray = [2, 2]
        self.splice_stroke_width = 2

        self.breakpoint_stroke_dasharray = [3, 3]
        self.breakpoint_orient_stroke_width = 2
        self.breakpoint_label_font_size = 20
        self.breakpoint_bottom_margin = 20
        self.breakpoint_top_margin = self.padding * 2 + self.breakpoint_label_font_size + self.breakpoint_bottom_margin
        self.breakpoint_label_prefix = 'B'

        self.marker_label_font_size = self.breakpoint_label_font_size
        self.marker_label_prefix = 'M'
        self.marker_top_margin = self.breakpoint_top_margin
        self.marker_bottom_margin = self.breakpoint_bottom_margin
        self.marker_color = self.breakpoint_color

        self.transcript_hyperlink = 'http://dec2013.archive.ensembl.org/Homo_sapiens/Transcript/Summary?' \
            'db=core;t={.name}'
        self.exon_font_size = 20
        self.exon_tear_tooth_width = 2
        self.exon_min_width = max([
            self.min_width + self.exon_tear_tooth_width * 2,
            self.exon_font_size * 2 * self.font_width_height_ratio
        ])
        self.exon_tear_tooth_height = 2
        self.exon_intron_ratio = 20
        self.exon1_color = self.gene1_color_selected
        self.exon2_color = self.gene2_color_selected

        self.transcript_label_prefix = 'T'
        self.fusion_label_prefix = 'F'

        self.translation_font_size = 14
        self.translation_scaffold_color = self.scaffold_color
        self.translation_track_height = self.translation_font_size
        self.translation_start_marker = 'M'
        self.translation_end_marker = '*'
        self.translation_marker_padding = 4

        self.legend_swatch_size = 50
        self.legend_font_size = 20
        self.legend_swatch_stroke = '#000000'
        self.legend_font_color = '#000000'
        self.legend_border_stroke = '#000000'
        self.legend_border_stroke_width = 1

        self.template_band_stroke_width = 0.5
        temp = [c.hex for c in Color('#ffffff').range_to(Color('#000000'), 7)]
        self.template_band_fill = {
            GIEMSA_STAIN.ACEN: '#800000',
            GIEMSA_STAIN.GPOS25: temp[1],
            GIEMSA_STAIN.GPOS33: temp[2],
            GIEMSA_STAIN.GPOS50: temp[3],
            GIEMSA_STAIN.GPOS66: temp[4],
            GIEMSA_STAIN.GPOS75: temp[5],
            GIEMSA_STAIN.GPOS100: temp[6],
            GIEMSA_STAIN.GNEG: '#ffffff'
        }
        self.template_band_stroke = '#000000'
        self.template_track_height = max([
            self.track_height / 3,
            self.label_font_size - self.breakpoint_bottom_margin - self.breakpoint_top_margin + self.breakpoint_label_font_size
        ])
        self.template_default_fill = '#ffffff'
        self.template_band_min_width = 2
        self.template_label_prefix = 'C'

        self.region_label_prefix = 'R'
        self.overlay_left_label = 16 * self.font_width_height_ratio * self.exon_font_size

        self.scatter_axis_font_size = 12
        self.scatter_error_bar_stroke_width = 1
        self.scatter_marker_radius = 2
        self.scatter_yaxis_tick_size = self.padding
        self.scatter_ytick_font_size = 10
