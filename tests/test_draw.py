import unittest
from structural_variant.draw import Diagram, HEX_BLACK, HEX_WHITE
from structural_variant.annotate.base import BioInterval
from structural_variant.annotate.genomic import Gene, Exon, IntergenicRegion, Template
from structural_variant.annotate.protein import Domain
from structural_variant.annotate.variant import Annotation, FusionTranscript
from svgwrite import Drawing
from structural_variant.constants import STRAND, ORIENT, SVTYPE, GIESMA_STAIN
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.interval import Interval
from tests import MockSeq, MockString, build_transcript
from copy import copy


class TestDraw(unittest.TestCase):
    def setUp(self):
        self.canvas = Drawing(height=100, width=1000)
        self.template_1 = Template('1', 1, 135006517)
        self.template_2 = Template('2', 1, 135006517)
        bands = []
        bands.append(BioInterval(None, 1, 2800001, 'p15.5'))
        bands.append(BioInterval(None, 2800001, 10700000, 'p15.4'))
        bands.append(BioInterval(None, 10700001, 12700000, 'p15.3'))
        bands.append(BioInterval(None, 12700001, 16200000, 'p15.2'))
        bands.append(BioInterval(None, 16200001, 21700000, 'p15.1'))
        bands.append(BioInterval(None, 21700001, 26100000, 'p14.3'))
        bands.append(BioInterval(None, 26100001, 27200000, 'p14.2'))
        bands.append(BioInterval(None, 27200001, 31000000, 'p14.1'))
        bands.append(BioInterval(None, 31000001, 36400000, 'p13'))
        bands.append(BioInterval(None, 36400001, 43500000, 'p12'))
        bands.append(BioInterval(None, 43500001, 48800000, 'p11.2'))
        bands.append(BioInterval(None, 48800001, 51600000, 'p11.12'))
        bands.append(BioInterval(None, 51600001, 53700000, 'p11.11', data={'giesma_stain': GIESMA_STAIN.ACEN}))
        bands.append(BioInterval(None, 53700001, 55700000, 'q11', data={'giesma_stain': GIESMA_STAIN.ACEN}))
        bands.append(BioInterval(None, 55700001, 59900000, 'q12.1'))
        bands.append(BioInterval(None, 59900001, 61700000, 'q12.2'))
        bands.append(BioInterval(None, 61700001, 63400000, 'q12.3'))
        bands.append(BioInterval(None, 63400001, 65900000, 'q13.1'))
        bands.append(BioInterval(None, 65900001, 68400000, 'q13.2'))
        bands.append(BioInterval(None, 68400001, 70400000, 'q13.3'))
        bands.append(BioInterval(None, 70400001, 75200000, 'q13.4'))
        bands.append(BioInterval(None, 75200001, 77100000, 'q13.5'))
        bands.append(BioInterval(None, 77100001, 85600000, 'q14.1'))
        bands.append(BioInterval(None, 85600001, 88300000, 'q14.2'))
        bands.append(BioInterval(None, 88300001, 92800000, 'q14.3'))
        bands.append(BioInterval(None, 92800001, 97200000, 'q21'))
        bands.append(BioInterval(None, 97200001, 102100000, 'q22.1'))
        bands.append(BioInterval(None, 102100001, 102900000, 'q22.2'))
        bands.append(BioInterval(None, 102900001, 110400000, 'q22.3'))
        bands.append(BioInterval(None, 110400001, 112500000, 'q23.1'))
        bands.append(BioInterval(None, 112500001, 114500000, 'q23.2'))
        bands.append(BioInterval(None, 114500001, 121200000, 'q23.3'))
        bands.append(BioInterval(None, 121200001, 123900000, 'q24.1'))
        bands.append(BioInterval(None, 123900001, 127800000, 'q24.2'))
        bands.append(BioInterval(None, 127800001, 130800000, 'q24.3'))
        bands.append(BioInterval(None, 130800001, 135006517, 'q25'))
        for b in bands:
            b.reference_object = self.template_1
            self.template_1.bands.append(b)
            b = copy(b)
            b.reference_object = self.template_2
            self.template_2.bands.append(b)


    def test__generate_interval_mapping(self):
        x = Interval(150, 1000)
        y = Interval(1500, 1950)
        z = Interval(5000, 7500)
        genic_length = 3803
        intergenic_length = 6197
        genic_intervals = 3
        intergenic_intervals = 5
        min_inter_width = 10
        min_width = 20

        temp = Diagram._generate_interval_mapping(
            [x, y, z], target_width=1000, ratio=5, min_width=min_width, start=1, end=10000, min_inter_width=min_inter_width)
        self.assertEqual(7, len(temp.keys()))
        expt = []
        st = 1
        end = st + min_inter_width + 4 - 1
        expt.append((Interval(1, 149), Interval(st, end)))
        st = end + 1
        end = st + min_width  + 168 - 1
        expt.append((Interval(150, 1000), Interval(st, end)))
        st = end + 1
        end  = st + min_inter_width + 12 - 1
        expt.append((Interval(1001, 1499), Interval(st, end)))
        st = end + 1
        end = st + min_width + 89 - 1
        expt.append((Interval(1500, 1950), Interval(st, end)))
        st = end + 1
        end = st + min_inter_width + 74 - 1
        expt.append((Interval(1951, 4999), Interval(st, end)))
        st = end + 1
        end = st + min_width + 493 - 1
        expt.append((Interval(5000, 7500), Interval(st, end)))
        st = end + 1
        end = st + min_inter_width + 61 - 1
        expt.append((Interval(7501, 10000), Interval(st, 1000)))
        for e, a in zip(expt, sorted(temp.items())):
            self.assertEqual(e, a)

    def test__generate_gene_mapping(self):
        d = Diagram()
        a = Gene('1', 1000, 2000)
        b = Gene('1', 5000, 7000)
        c = Gene('1', 1500, 2500)
        genes = [a, b, c]
        """return self._generate_interval_mapping(
            target_width,
            genes,
            self.GENE_INTERGENIC_RATIO,
            self.MIN_WIDTH + self.GENE_ARROW_WIDTH,
            buffer=self.GENE_MIN_BUFFER
        )
        with self.assertRaises(AttributeError):
            m = d._generate_gene_mapping(100, genes)
        m = d._generate_gene_mapping(500, genes)
        u = Interval.union(*m.values())
        self.assertLessEqual(1, u.start)
        self.assertGreaterEqual(500, u.end)"""

    def test__generate_gene_mapping_err(self):
        #  _generate_interval_mapping [IntergenicRegion(11:77361962_77361962+)] 1181.39453125 5 30 None 77356962 77366962)
        ir = IntergenicRegion('11', 5000, 5000, STRAND.POS)
        tgt_width = 1000
        d = Diagram()
        d.GENE_MIN_BUFFER = 10
        # (self, canvas, gene, width, height, fill, label='', REFERENCE_GENOME=None)
        d.draw_genes(self.canvas, [ir], tgt_width, [])

        # _generate_interval_mapping ['Interval(29684391, 29684391)', 'Interval(29663998, 29696515)'] 1181.39453125 5 60 None 29662998 29697515
        # def _generate_interval_mapping(cls, input_intervals, target_width, ratio, min_width, buffer_length=None, start=None, end=None, min_inter_width=None)
        itvls = [Interval(29684391, 29684391), Interval(29663998, 29696515)]
        mapping = Diagram._generate_interval_mapping(
            itvls, 1181.39, 5, 60, None, 29662998, 29697515)



    def test__split_intervals_into_tracks(self):
        # ----======---------
        # ------======--------
        # -----===============
        t = Diagram._split_intervals_into_tracks(
            [(1, 3), (3, 7), (2, 2), (4, 5), (3, 10)]
        )
        self.assertEqual(3, len(t))
        self.assertEqual([(1, 3), (4, 5)], t[0])
        self.assertEqual([(2, 2), (3, 7)], t[1])
        self.assertEqual([(3, 10)], t[2])

    def test_draw_genes(self):

        x = Gene('1', 1000, 2000, strand=STRAND.POS)
        y = Gene('1', 5000, 7000, strand=STRAND.NEG)
        z = Gene('1', 1500, 2500, strand=STRAND.POS)

        d = Diagram()
        breakpoints = [
            Breakpoint('1', 1100, 1200, orient=ORIENT.RIGHT)
        ]
        g = d.draw_genes(
            self.canvas, [x, y, z], 500, breakpoints, {x: d.GENE1_COLOR, y: d.GENE2_COLOR_SELECTED, z: d.GENE2_COLOR})

        # test the class structure
        self.assertEqual(5, len(g.elements))
        self.assertEqual('scaffold', g.elements[0].attribs.get('class', ''))
        for i in range(1, 4):
            self.assertEqual('gene', g.elements[i].attribs.get('class', ''))
        self.assertEqual('breakpoint', g.elements[4].attribs.get('class', ''))
        self.assertEqual(
            d.TRACK_HEIGHT * 2 + d.PADDING + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN,
            g.height
        )
        self.canvas.add(g)
        self.assertEqual(len(g.labels), 4)
        self.assertEqual(x, g.labels['G1'])
        self.assertEqual(z, g.labels['G2'])
        self.assertEqual(y, g.labels['G3'])
        self.assertEqual(breakpoints[0], g.labels['B1'])

    def test_draw_ustranscript(self):
        d = Diagram()
        # domains = [Domain()]
        d1 = Domain('first', [(55, 61), (71, 73)])
        d2 = Domain('second', [(10, 20), (30, 34)])

        t = build_transcript(
            gene=None,
            cds_start=50,
            cds_end=249,
            exons=[(1, 99), (200, 299), (400, 499)],
            strand=STRAND.NEG,
            domains=[d2, d1]
        )
        b = Breakpoint('1', 350, 410, orient=ORIENT.LEFT)
        g = d.draw_ustranscript(
            self.canvas, t, 500,
            colors={t.exons[1]: '#FFFF00'},
            breakpoints=[b]
        )
        self.canvas.add(g)
        self.assertEqual(2, len(self.canvas.elements))
        self.assertEqual(4, len(g.elements))
        for el, cls in zip(g.elements, ['splicing', 'exon_track', 'protein', 'breakpoint']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        for el, cls in zip(g.elements[1].elements, ['scaffold', 'exon', 'exon', 'exon']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        for el, cls in zip(g.elements[2].elements, ['translation', 'domain', 'domain']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        self.assertEqual(
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT
            + 2 * d.PADDING + d.DOMAIN_TRACK_HEIGHT * 2
            + d.TRANSLATION_TRACK_HEIGHT + d.PADDING
            + d.BREAKPOINT_TOP_MARGIN
            + d.BREAKPOINT_BOTTOM_MARGIN,
            g.height)
        self.assertEqual(d1, g.labels['D1'])
        self.assertEqual(d2, g.labels['D2'])

    def test_dynamic_label_color(self):
        self.assertEqual(HEX_WHITE, Diagram.dynamic_label_color(HEX_BLACK))
        self.assertEqual(HEX_BLACK, Diagram.dynamic_label_color(HEX_WHITE))

    def test_draw_legend(self):
        d = Diagram()
        swatches = [
            ('#000000', 'black'),
            ('#FF0000', 'red'),
            ('#0000FF', 'blue'),
            ('#00FF00', 'green'),
            ('#FFFF00', 'yellow')
        ]
        g = d.draw_legend(self.canvas, swatches)
        self.canvas.add(g)

        self.assertEqual('legend', g.attribs.get('class', ''))
        self.assertEqual(
            d.LEGEND_SWATCH_SIZE * len(swatches) + d.PADDING * (len(swatches) - 1 + 2),
            g.height
        )
        self.assertEqual(6, len(g.elements))
        self.assertEqual(
            6 * d.LEGEND_FONT_SIZE * d.FONT_WIDTH_HEIGHT_RATIO + d.PADDING * 3 + d.LEGEND_SWATCH_SIZE,
            g.width
        )

    def test_draw_layout_single_transcript(self):
        d = Diagram()
        d1 = Domain('first', [(55, 61), (71, 73)])
        d2 = Domain('second', [(10, 20), (30, 34)])
        g1 = Gene('1', 150, 1000, strand=STRAND.POS)
        t = build_transcript(g1, [(200, 299), (400, 499), (700, 899)], 50, 249, [d2, d1])
        b1 = Breakpoint('1', 350, orient=ORIENT.RIGHT)
        b2 = Breakpoint('1', 600, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_sequence='')
        ann = Annotation(bpp, transcript1=t, transcript2=t, event_type=SVTYPE.DUP)
        ann.add_gene(Gene('1', 1500, 1950, strand=STRAND.POS))

        reference_genome = {'1': MockSeq(MockString('A'))}
        ft = FusionTranscript.build(ann, reference_genome)

        canvas = d.draw(ann, ft)
        self.assertEqual(4, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_layout_single_genomic(self):
        d = Diagram()
        d1 = Domain('first', [(55, 61), (71, 73)])
        d2 = Domain('second', [(10, 20), (30, 34)])
        g1 = Gene('1', 150, 1000, strand=STRAND.POS)
        g2 = Gene('1', 5000, 7500, strand=STRAND.POS)
        t1 = build_transcript(
            gene=g1,
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1]
        )
        t2 = build_transcript(
            gene=g2,
            cds_start=20,
            cds_end=500,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[]
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('1', 6500, orient=ORIENT.RIGHT)
        bpp = BreakpointPair(b1, b2, opposing_strands=False, untemplated_sequence='')
        ann = Annotation(bpp, transcript1=t1, transcript2=t2)
        ann.add_gene(Gene('1', 1500, 1950, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {'1': MockSeq(MockString('A'))}

        ft = FusionTranscript.build(ann, reference_genome)
        self.assertEqual(t1.exons[0], ft.exon_mapping[ft.exons[0]])
        self.assertEqual(t2.exons[2], ft.exon_mapping[ft.exons[1]])
        self.assertEqual(t2.exons[3], ft.exon_mapping[ft.exons[2]])

        canvas = d.draw(ann, ft)
        self.assertEqual(5, len(canvas.elements))  # defs counts as element

        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_area_plot(self):
        d = Diagram()
        mapping = {Interval(1, 100): Interval(1, 100)}
        data = [
            (1, 10),
            (2, 10),
            (10, 20),
            (40, 50),
            (67, 0),
            (75, 100),
            (85, 150),
            (95, 120)
        ]
        g = d.draw_area_plot(self.canvas, data, 100, '#FF0000')
        self.canvas.add(g)
        self.assertEqual(2, len(self.canvas.elements))

    def test_draw_layout_translocation(self):
        d = Diagram()
        d1 = Domain('first', [(55, 61), (71, 73)])
        d2 = Domain('second', [(10, 20), (30, 34)])
        g1 = Gene('1', 150, 1000, strand=STRAND.POS)
        g2 = Gene('2', 5000, 7500, strand=STRAND.NEG)
        t1 = build_transcript(
            gene=g1,
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1]
        )
        t2 = build_transcript(
            gene=g2,
            cds_start=120,
            cds_end=700,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[]
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('2', 6520, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_sequence='')
        ann = Annotation(bpp, transcript1=t1, transcript2=t2)
        # genes 1
        ann.add_gene(Gene('1', 1500, 1950, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3700, 4400, strand=STRAND.NEG))
        # genes 2
        ann.add_gene(Gene('2', 1500, 1950, strand=STRAND.NEG))
        ann.add_gene(Gene('2', 5500, 9000, strand=STRAND.POS))
        ann.add_gene(Gene('2', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {'1': MockSeq(MockString('A')), '2': MockSeq(MockString('A'))}

        ft = FusionTranscript.build(ann, reference_genome)

        canvas = d.draw(ann, ft)
        self.assertEqual(6, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT * 2 + d.PADDING  + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_template(self):
        # def draw_template(self, canvas, template, target_width, height, labels=None, colors=None):
        d = Diagram()
        canvas = Drawing(size=(1000, 50))
        t = Template(
            '1', 1, 100000,
            bands=[
                BioInterval(None, 1, 8000, 'p1'),
                BioInterval(None, 10000, 15000, 'p2')
            ])
        g = d.draw_template(canvas, t, 1000, 50)
        canvas.add(g)
        canvas.attribs['height'] = g.height
        canvas = Drawing(size=(1000, 50))

        g = d.draw_template(canvas, self.template_1, 1000, 50)
        self.assertEqual(d.BREAKPOINT_TOP_MARGIN + d.BREAKPOINT_BOTTOM_MARGIN + d.TEMPLATE_TRACK_HEIGHT, g.height)
        canvas.add(g)
        canvas.attribs['height'] = g.height
        self.assertEqual(2, len(canvas.elements))


    def test_draw_translocation_with_template(self):
        d = Diagram()
        d1 = Domain('first', [(55, 61), (71, 73)])
        d2 = Domain('second', [(10, 20), (30, 34)])
        g1 = Gene(self.template_1, 150, 1000, strand=STRAND.POS)
        g2 = Gene(self.template_2, 5000, 7500, strand=STRAND.NEG)
        templates = {self.template_1.name: self.template_1, self.template_2.name: self.template_2}
        t1 = build_transcript(
            gene=g1,
            cds_start=50,
            cds_end=249,
            exons=[(200, 299), (400, 499), (700, 899)],
            domains=[d2, d1]
        )
        t2 = build_transcript(
            gene=g2,
            cds_start=120,
            cds_end=700,
            exons=[(5100, 5299), (5800, 6199), (6500, 6549), (6700, 6799)],
            domains=[]
        )
        b1 = Breakpoint('1', 350, orient=ORIENT.LEFT)
        b2 = Breakpoint('2', 6520, orient=ORIENT.LEFT)
        bpp = BreakpointPair(b1, b2, opposing_strands=True, untemplated_sequence='')
        ann = Annotation(bpp, transcript1=t1, transcript2=t2)
        # genes 1
        ann.add_gene(Gene('1', 1500, 1950, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3700, 4400, strand=STRAND.NEG))
        # genes 2
        ann.add_gene(Gene('2', 1500, 1950, strand=STRAND.NEG))
        ann.add_gene(Gene('2', 5500, 9000, strand=STRAND.POS))
        ann.add_gene(Gene('2', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {'1': MockSeq(MockString('A')), '2': MockSeq(MockString('A'))}

        ft = FusionTranscript.build(ann, reference_genome)

        canvas = d.draw(ann, ft, draw_template=True, templates=templates)
        canvas.saveas('test_figure.svg')
        self.assertEqual(8, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT * 2 + d.PADDING  + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.TEMPLATE_TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_overlay(self):
        gene = Gene('12', 25357723, 25403870, strand=STRAND.NEG, name='KRAS')
        marker = BioInterval('12', 25403865, name='splice site mutation')
        t = build_transcript(
            cds_start=193, cds_end=759,
            exons=[
                Exon(25403685, 25403865),
                Exon(25398208, 25398329),
                Exon(25380168, 25380346),
                Exon(25378548, 25378707),
                Exon(25357723, 25362845)],
            gene=gene, domains=[])
        build_transcript(
            cds_start=198, cds_end=425,
            exons=[Exon(25403685, 25403870), Exon(25398208, 25398329), Exon(25362102, 25362845)],
            gene=gene, domains=[])
        build_transcript(
            cds_start=65, cds_end=634,
            exons=[
                Exon(25403685, 25403737),
                Exon(25398208, 25398329),
                Exon(25380168, 25380346),
                Exon(25378548, 25378707),
                Exon(25368371, 25368494),
                Exon(25362365, 25362845)],
            gene=gene, domains=[],
            is_best_transcript=True)
        build_transcript(
            cds_start=65, cds_end=634,
            exons=[Exon(25403698, 25403863), Exon(25398208, 25398329), Exon(25386753, 25388160)],
            gene=gene, domains=[])
        d = Diagram()
        d.GENE_MIN_BUFFER = 0
        canvas = d.draw_ustranscripts_overlay(gene, markers=[marker])
        self.assertEqual(2, len(canvas.elements))  # defs counts as element
        canvas.saveas('test_overlay_figure.svg')
        raise unittest.SkipTest('TODO. add height calculation assert')
