import random

import pytest
from svgwrite import Drawing

from mavis.annotate import fusion, genomic, protein, variant
from mavis.annotate.base import BioInterval
from mavis.annotate.file_io import load_templates
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, PROTOCOL, STRAND, SVTYPE
from mavis.illustrate.constants import DEFAULTS, DiagramSettings
from mavis.illustrate.diagram import (
    HEX_BLACK,
    HEX_WHITE,
    draw_multi_transcript_overlay,
    draw_sv_summary_diagram,
    generate_interval_mapping,
)
from mavis.illustrate.elements import (
    draw_genes,
    draw_legend,
    draw_template,
    draw_ustranscript,
)
from mavis.illustrate.scatter import ScatterPlot
from mavis.illustrate.util import dynamic_label_color, split_intervals_into_tracks
from mavis.interval import Interval

from ...util import get_data
from ..mock import OUTPUT_SVG, MockObject, MockString, build_transcript

TEMPLATE_METADATA = None
DEFAULTS.domain_name_regex_filter = r'.*'


def setUpModule():
    global TEMPLATE_METADATA, EXAMPLE_ANNOTATIONS
    TEMPLATE_METADATA = load_templates(get_data('cytoBand.txt'))


@pytest.fixture
def canvas():
    return Drawing(height=100, width=1000)


class TestDraw:
    def test_generate_interval_mapping_outside_range_error(self):
        temp = [
            Interval(48556470, 48556646),
            Interval(48573290, 48573665),
            Interval(48575056, 48575078),
        ]
        mapping = generate_interval_mapping(
            input_intervals=temp,
            target_width=431.39453125,
            ratio=20,
            min_width=14,
            buffer_length=None,
            end=None,
            start=None,
            min_inter_width=10,
        )
        st = min([x.start for x in temp])
        end = min([x.end for x in temp])
        Interval.convert_pos(mapping, st)
        Interval.convert_pos(mapping, end)

    def test_generate_gene_mapping_err(self, canvas):
        #  _generate_interval_mapping [genomic.IntergenicRegion(11:77361962_77361962+)] 1181.39453125 5 30 None 77356962 77366962)
        ir = genomic.IntergenicRegion('11', 5000, 5000, STRAND.POS)
        tgt_width = 1000
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        d.gene_min_buffer = 10
        # (self, canvas, gene, width, height, fill, label='', reference_genome=None)
        draw_genes(d, canvas, [ir], tgt_width, [])

        # _generate_interval_mapping ['Interval(29684391, 29684391)', 'Interval(29663998, 29696515)'] 1181.39453125 5 60 None 29662998 29697515
        # def generate_interval_mapping(cls, input_intervals, target_width, ratio, min_width, buffer_length=None, start=None, end=None, min_inter_width=None)
        itvls = [Interval(29684391, 29684391), Interval(29663998, 29696515)]
        generate_interval_mapping(itvls, 1181.39, 5, 60, None, 29662998, 29697515)

    def test_split_intervals_into_tracks(self):
        # ----======---------
        # ------======--------
        # -----===============
        t = split_intervals_into_tracks([(1, 3), (3, 7), (2, 2), (4, 5), (3, 10)])
        assert len(t) == 3
        assert t[0] == [(1, 3), (4, 5)]
        assert t[1] == [(2, 2), (3, 7)]
        assert t[2] == [(3, 10)]

    def test_draw_genes(self, canvas):
        x = genomic.Gene('1', 1000, 2000, strand=STRAND.POS)
        y = genomic.Gene('1', 5000, 7000, strand=STRAND.NEG)
        z = genomic.Gene('1', 1500, 2500, strand=STRAND.POS)

        d = DiagramSettings(domain_name_regex_filter=r'.*')
        breakpoints = [Breakpoint('1', 1100, 1200, orient=ORIENT.RIGHT)]
        g = draw_genes(
            d,
            canvas,
            [x, y, z],
            500,
            breakpoints,
            {x: d.gene1_color, y: d.gene2_color_selected, z: d.gene2_color},
        )

        # test the class structure
        assert len(g.elements) == 6
        assert g.elements[0].attribs.get('class', '') == 'scaffold'
        for i in range(1, 4):
            assert g.elements[i].attribs.get('class', '') == 'gene'
        assert g.elements[4].attribs.get('class', '') == 'mask'
        assert g.elements[5].attribs.get('class', '') == 'breakpoint'
        assert (
            g.height
            == d.track_height * 2 + d.padding + d.breakpoint_bottom_margin + d.breakpoint_top_margin
        )
        canvas.add(g)
        assert 4 == len(g.labels)
        assert g.labels['G1'] == x
        assert g.labels['G2'] == z
        assert g.labels['G3'] == y
        assert g.labels['B1'] == breakpoints[0]

    def test_draw_ustranscript(self, canvas):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        # domains = [protein.Domain()]
        d1 = protein.Domain('first', [(55, 61), (71, 73)])
        d2 = protein.Domain('second', [(10, 20), (30, 34)])

        t = build_transcript(
            gene=None,
            cds_start=50,
            cds_end=249,
            exons=[(1, 99), (200, 299), (400, 499)],
            strand=STRAND.NEG,
            domains=[d2, d1],
        )
        b = Breakpoint('1', 350, 410, orient=ORIENT.LEFT)
        g = draw_ustranscript(d, canvas, t, 500, colors={t.exons[1]: '#FFFF00'}, breakpoints=[b])
        canvas.add(g)
        # canvas.saveas('test_draw_ustranscript.svg')
        assert len(canvas.elements) == 2
        assert len(g.elements) == 3
        for el, cls in zip(g.elements[0].elements, ['splicing', 'exon_track', 'protein']):
            assert el.attribs.get('class', '') == cls

        for el, cls in zip(
            g.elements[0].elements[1].elements, ['scaffold', 'exon', 'exon', 'exon']
        ):
            assert el.attribs.get('class', '') == cls

        for el, cls in zip(g.elements[0].elements[2].elements, ['translation', 'domain', 'domain']):
            assert el.attribs.get('class', '') == cls

        assert g.height == sum(
            [
                d.track_height,
                d.splice_height,
                2 * d.padding,
                d.domain_track_height * 2,
                d.translation_track_height,
                d.padding,
                d.breakpoint_top_margin,
                d.breakpoint_bottom_margin,
            ]
        )
        assert g.labels['D1'] == d1.name
        assert g.labels['D2'] == d2.name

    def test_draw_consec_exons(self, canvas):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        # domains = [protein.Domain()]
        t = build_transcript(
            gene=None,
            cds_start=50,
            cds_end=249,
            exons=[(1, 99), (200, 299), (300, 350), (400, 499)],
            strand=STRAND.POS,
            domains=[],
        )
        b = Breakpoint('1', 350, 410, orient=ORIENT.LEFT)
        g = draw_ustranscript(d, canvas, t, 500, colors={t.exons[1]: '#FFFF00'}, breakpoints=[b])
        canvas.add(g)
        if OUTPUT_SVG:
            canvas.saveas('test_draw_consec_exons.svg')

        # canvas.saveas('test_draw_ustranscript.svg')
        assert len(canvas.elements) == 2
        assert len(g.elements) == 3
        # check that only 2 splicing marks were created
        assert len(g.elements[0].elements[0].elements) == 2
        # get the second exon
        ex2 = g.elements[0].elements[1].elements[2].elements[0]
        print(ex2)
        assert pytest.approx(ex2.attribs.get('width')) == 120.7783426339
        # get the third exon
        ex3 = g.elements[0].elements[1].elements[3].elements[0]
        print(ex3)
        assert pytest.approx(ex3.attribs.get('width')) == 96.52494419642852

    def test_dynamic_label_color(self):
        assert dynamic_label_color(HEX_BLACK) == HEX_WHITE
        assert dynamic_label_color(HEX_WHITE) == HEX_BLACK

    def test_draw_legend(self, canvas):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        swatches = [
            ('#000000', 'black'),
            ('#FF0000', 'red'),
            ('#0000FF', 'blue'),
            ('#00FF00', 'green'),
            ('#FFFF00', 'yellow'),
        ]
        g = draw_legend(d, canvas, swatches)
        canvas.add(g)

        assert g.attribs.get('class', '') == 'legend'
        assert g.height == d.legend_swatch_size * len(swatches) + d.padding * (
            len(swatches) - 1 + 2
        )

        assert len(g.elements) == 6
        assert (
            g.width
            == 6 * d.legend_font_size * d.font_width_height_ratio
            + d.padding * 3
            + d.legend_swatch_size
        )

    def test_draw_layout_single_transcript(self):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        d1 = protein.Domain('first', [(55, 61), (71, 73)])
        d2 = protein.Domain('second', [(10, 20), (30, 34)])
        g1 = genomic.Gene('1', 150, 1000, strand=STRAND.POS)
        t = build_transcript(g1, [(200, 299), (400, 499), (700, 899)], 50, 249, [d2, d1])
        b1 = Breakpoint('1', 350, orient=ORIENT.RIGHT)
        b2 = Breakpoint('1', 600, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='')
        ann = variant.Annotation(
            bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP, protocol=PROTOCOL.GENOME
        )
        ann.add_gene(genomic.Gene('1', 1500, 1950, strand=STRAND.POS))

        reference_genome = {'1': MockObject(seq=MockString('A'))}
        ft = variant.FusionTranscript.build(ann, reference_genome)
        ann.fusion = ft
        canvas, legend = draw_sv_summary_diagram(d, ann)
        assert len(canvas.elements) == 4  # defs counts as element
        expected_height = (
            d.top_margin
            + d.bottom_margin
            + d.track_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.inner_margin
            + d.track_height
            + d.splice_height
            + d.padding
            + d.translation_track_height
            + d.padding * 2
            + d.domain_track_height * 2
            + d.inner_margin
            + d.track_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.splice_height
        )
        if OUTPUT_SVG:
            canvas.saveas('test_draw_layout_single_transcript.svg')
        assert canvas.attribs['height'] == expected_height

    def test_draw_layout_single_genomic(self):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        d1 = protein.Domain('first', [(55, 61), (71, 73)])
        d2 = protein.Domain('second', [(10, 20), (30, 34)])
        g1 = genomic.Gene('1', 150, 1000, strand=STRAND.POS)
        g2 = genomic.Gene('1', 5000, 7500, strand=STRAND.POS)
        t1 = build_transcript(
            gene=g1,
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1],
        )
        t2 = build_transcript(
            gene=g2,
            cds_start=20,
            cds_end=500,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[],
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('1', 6500, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_seq='')
        ann = variant.Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME
        )
        ann.add_gene(genomic.Gene('1', 1500, 1950, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('1', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {'1': MockObject(seq=MockString('A'))}

        ft = variant.FusionTranscript.build(ann, reference_genome)
        ann.fusion = ft
        assert ft.exon_mapping[ft.exons[0].position] == t1.exons[0]
        assert ft.exon_mapping[ft.exons[1].position] == t2.exons[2]
        assert ft.exon_mapping[ft.exons[2].position] == t2.exons[3]

        canvas, legend = draw_sv_summary_diagram(d, ann)
        assert len(canvas.elements) == 5  # defs counts as element

        expected_height = (
            d.top_margin
            + d.bottom_margin
            + d.track_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.inner_margin
            + d.track_height
            + d.splice_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.padding
            + d.translation_track_height
            + d.padding * 2
            + d.domain_track_height * 2
            + d.inner_margin
            + d.track_height
            + d.splice_height
        )
        assert canvas.attribs['height'] == expected_height
        if OUTPUT_SVG:
            canvas.saveas('test_draw_layout_single_genomic.svg')

    def test_draw_layout_translocation(self):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        d1 = protein.Domain('first', [(55, 61), (71, 73)])
        d2 = protein.Domain('second', [(10, 20), (30, 34)])
        g1 = genomic.Gene('1', 150, 1000, strand=STRAND.POS)
        g2 = genomic.Gene('2', 5000, 7500, strand=STRAND.NEG)
        t1 = build_transcript(
            gene=g1,
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1],
        )
        t2 = build_transcript(
            gene=g2,
            cds_start=120,
            cds_end=700,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[],
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('2', 6520, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='')
        ann = variant.Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.ITRANS, protocol=PROTOCOL.GENOME
        )
        # genes 1
        ann.add_gene(genomic.Gene('1', 1500, 1950, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('1', 3700, 4400, strand=STRAND.NEG))
        # genes 2
        ann.add_gene(genomic.Gene('2', 1500, 1950, strand=STRAND.NEG))
        ann.add_gene(genomic.Gene('2', 5500, 9000, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('2', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {
            '1': MockObject(seq=MockString('A')),
            '2': MockObject(seq=MockString('A')),
        }

        ft = variant.FusionTranscript.build(ann, reference_genome)
        ann.fusion = ft
        canvas, legend = draw_sv_summary_diagram(d, ann)
        assert len(canvas.elements) == 6  # defs counts as element
        expected_height = (
            d.top_margin
            + d.bottom_margin
            + d.track_height * 2
            + d.padding
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.inner_margin
            + d.track_height
            + d.splice_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.padding
            + d.translation_track_height
            + d.padding * 2
            + d.domain_track_height * 2
            + d.inner_margin
            + d.track_height
            + d.splice_height
        )
        assert canvas.attribs['height'] == expected_height

    def test_draw_template(self):
        # def draw_template(self, canvas, template, target_width, height, labels=None, colors=None):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        canvas = Drawing(size=(1000, 50))
        t = genomic.Template(
            '1',
            1,
            100000,
            bands=[BioInterval(None, 1, 8000, 'p1'), BioInterval(None, 10000, 15000, 'p2')],
        )
        g = draw_template(d, canvas, t, 1000)
        canvas.add(g)
        canvas.attribs['height'] = g.height
        canvas = Drawing(size=(1000, 50))

        g = draw_template(d, canvas, TEMPLATE_METADATA['1'], 1000)
        assert (
            g.height
            == d.breakpoint_top_margin + d.breakpoint_bottom_margin + d.template_track_height
        )
        canvas.add(g)
        canvas.attribs['height'] = g.height
        assert len(canvas.elements) == 2

    def test_draw_translocation_with_template(self):
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        d1 = protein.Domain('PF0001', [(55, 61), (71, 73)])
        d2 = protein.Domain('PF0002', [(10, 20), (30, 34)])
        g1 = genomic.Gene(TEMPLATE_METADATA['1'], 150, 1000, strand=STRAND.POS, aliases=['HUGO2'])
        g2 = genomic.Gene(TEMPLATE_METADATA['X'], 5000, 7500, strand=STRAND.NEG, aliases=['HUGO3'])
        t1 = build_transcript(
            gene=g1,
            name='transcript1',
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1],
        )
        t2 = build_transcript(
            gene=g2,
            name='transcript2',
            cds_start=120,
            cds_end=700,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[],
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('2', 6520, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_seq='')
        ann = variant.Annotation(
            bpp, transcript1=t1, transcript2=t2, event_type=SVTYPE.ITRANS, protocol=PROTOCOL.GENOME
        )
        # genes 1
        ann.add_gene(genomic.Gene('1', 1500, 1950, strand=STRAND.POS, aliases=['HUGO5']))
        ann.add_gene(genomic.Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('1', 3700, 4400, strand=STRAND.NEG))
        # genes 2
        ann.add_gene(genomic.Gene('2', 1500, 1950, strand=STRAND.NEG))
        ann.add_gene(genomic.Gene('2', 5500, 9000, strand=STRAND.POS))
        ann.add_gene(genomic.Gene('2', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {
            '1': MockObject(seq=MockString('A')),
            '2': MockObject(seq=MockString('A')),
        }

        ft = variant.FusionTranscript.build(ann, reference_genome)
        ann.fusion = ft
        canvas, legend = draw_sv_summary_diagram(
            d, ann, draw_reference_templates=True, templates=TEMPLATE_METADATA
        )
        if OUTPUT_SVG:
            canvas.saveas('test_draw_translocation_with_template.svg')
        assert len(canvas.elements) == 8  # defs counts as element
        expected_height = (
            d.top_margin
            + d.bottom_margin
            + d.track_height * 2
            + d.padding
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.inner_margin
            + d.track_height
            + d.splice_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.padding
            + d.translation_track_height
            + d.padding * 2
            + d.domain_track_height * 2
            + d.inner_margin * 2
            + d.track_height
            + d.breakpoint_bottom_margin
            + d.breakpoint_top_margin
            + d.splice_height
            + d.template_track_height
        )
        assert pytest.approx(canvas.attribs['height']) == expected_height

    def test_draw_overlay(self):
        gene = genomic.Gene('12', 25357723, 25403870, strand=STRAND.NEG, name='KRAS')
        marker = BioInterval('12', 25403865, name='splice site mutation')
        t = build_transcript(
            cds_start=193,
            cds_end=759,
            exons=[
                (25403685, 25403865),
                (25398208, 25398329),
                (25380168, 25380346),
                (25378548, 25378707),
                (25357723, 25362845),
            ],
            gene=gene,
            domains=[],
        )
        build_transcript(
            cds_start=198,
            cds_end=425,
            exons=[(25403685, 25403870), (25398208, 25398329), (25362102, 25362845)],
            gene=gene,
            domains=[],
        )
        build_transcript(
            cds_start=65,
            cds_end=634,
            exons=[
                (25403685, 25403737),
                (25398208, 25398329),
                (25380168, 25380346),
                (25378548, 25378707),
                (25368371, 25368494),
                (25362365, 25362845),
            ],
            gene=gene,
            domains=[protein.Domain('domain1', [(1, 10)]), protein.Domain('domain1', [(4, 10)])],
            is_best_transcript=True,
        )
        build_transcript(
            cds_start=65,
            cds_end=634,
            exons=[(25403698, 25403863), (25398208, 25398329), (25386753, 25388160)],
            gene=gene,
            domains=[],
        )
        d = DiagramSettings(domain_name_regex_filter=r'.*')
        for i, t in enumerate(gene.transcripts):
            t.name = 'transcript {}'.format(i + 1)
        scatterx = [x + 100 for x in range(gene.start, gene.end + 1, 400)]
        scattery = [random.uniform(-0.2, 0.2) for x in scatterx]
        s = ScatterPlot(list(zip(scatterx, scattery)), 'cna', ymin=-1, ymax=1, yticks=[-1, 0, 1])

        d.gene_min_buffer = 0
        canvas = draw_multi_transcript_overlay(d, gene, vmarkers=[marker], plots=[s, s])
        assert len(canvas.elements) == 2  # defs counts as element
        if OUTPUT_SVG:
            canvas.saveas('test_draw_overlay.svg')

    def test_single_bp_ins_exon(self):
        transcript = fusion.FusionTranscript()
        transcript.position = Interval(401258, 408265)
        transcript.exons = [
            genomic.Exon(401258, 401461, transcript=transcript),
            genomic.Exon(404799, 405254, intact_end_splice=False, transcript=transcript),
            genomic.Exon(
                405255,
                405255,
                intact_start_splice=False,
                intact_end_splice=False,
                transcript=transcript,
            ),
            genomic.Exon(405256, 408265, intact_start_splice=False, transcript=transcript),
        ]

        cfg = DiagramSettings(width=1500)
        canvas = Drawing(size=(cfg.width, 1000))
        drawing_width = cfg.width - cfg.label_left_margin - cfg.left_margin - cfg.right_margin
        canvas.add(draw_ustranscript(cfg, canvas, transcript, target_width=drawing_width))
        if OUTPUT_SVG:
            canvas.saveas('test_single_bp_ins_exon.svg')

    def test_single_bp_dup(self):
        transcript = fusion.FusionTranscript()
        transcript.position = Interval(1, 500)
        transcript.exons = [
            genomic.Exon(1, 7, transcript=transcript, intact_end_splice=False),
            genomic.Exon(
                8, 8, transcript=transcript, intact_start_splice=False, intact_end_splice=False
            ),
            genomic.Exon(9, 100, transcript=transcript, intact_start_splice=False),
            genomic.Exon(200, 500, transcript=transcript),
        ]
        cfg = DiagramSettings(width=1500)
        canvas = Drawing(size=(cfg.width, 1000))
        drawing_width = cfg.width - cfg.label_left_margin - cfg.left_margin - cfg.right_margin
        canvas.add(draw_ustranscript(cfg, canvas, transcript, target_width=drawing_width))
        if OUTPUT_SVG:
            canvas.saveas('test_single_bp_dup.svg')

    def test_two_exon_transcript(self):
        transcript = fusion.FusionTranscript()
        transcript.position = Interval(1, 555)
        transcript.exons = [
            genomic.Exon(55820038, 55820969, transcript=transcript),
            genomic.Exon(55820971, 55820976, transcript=transcript),
        ]
        transcript.exon_mapping[Interval(55820971, 55820976)] = MockObject(
            transcript=MockObject(exon_number=lambda *x: 2)
        )
        cfg = DiagramSettings(width=1500)
        canvas = Drawing(size=(cfg.width, 1000))
        drawing_width = cfg.width - cfg.label_left_margin - cfg.left_margin - cfg.right_margin
        canvas.add(draw_ustranscript(cfg, canvas, transcript, target_width=drawing_width))
        if OUTPUT_SVG:
            canvas.saveas('test_two_exon_transcript.svg')
