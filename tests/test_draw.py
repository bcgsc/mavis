import unittest
from structural_variant.illustrate.settings import DiagramSettings
from structural_variant.illustrate.scatter import ScatterPlot
from structural_variant.illustrate.draw import *
from structural_variant.illustrate.draw import _draw_template, _draw_ustranscript, _draw_genes
from structural_variant.annotate import *
from svgwrite import Drawing
from structural_variant.constants import STRAND, ORIENT, SVTYPE
from structural_variant.breakpoint import Breakpoint, BreakpointPair
from structural_variant.interval import Interval
from tests import MockSeq, MockString, build_transcript, TEMPLATE_METADATA_FILE
import random

TEMPLATE_METADATA = None


def setUpModule():
    global TEMPLATE_METADATA
    TEMPLATE_METADATA = load_templates(TEMPLATE_METADATA_FILE)


class TestDraw(unittest.TestCase):
    def setUp(self):
        self.canvas = Drawing(height=100, width=1000)

    def test_generate_interval_mapping(self):
        x = Interval(150, 1000)
        y = Interval(1500, 1950)
        z = Interval(5000, 7500)
        genic_length = 3803
        intergenic_length = 6197
        genic_intervals = 3
        intergenic_intervals = 5
        min_inter_width = 10
        min_width = 20

        temp = generate_interval_mapping(
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
        actual = sorted(temp.items())
        for e, a in zip(expt, actual):
            self.assertEqual(e[0], a[0])
        self.assertAlmostEqual(1, actual[0][1].start)
        self.assertAlmostEqual(1000, actual[-1][1].end)

    def test_generate_interval_mapping_outside_range_error(self):
        temp = [
            Interval(48556470, 48556646),
            Interval(48573290, 48573665),
            Interval(48575056, 48575078)
        ]
        mapping = generate_interval_mapping(
            input_intervals=temp,
            target_width=431.39453125,
            ratio=20,
            min_width=14,
            buffer_length=None,
            end=None,
            start=None,
            min_inter_width=10
        )
        st = min([x.start for x in temp])
        end = min([x.end for x in temp])
        Interval.convert_pos(mapping, st)
        Interval.convert_pos(mapping, end)

    def test_generate_gene_mapping(self):
        d = DiagramSettings()
        a = Gene('1', 1000, 2000)
        b = Gene('1', 5000, 7000)
        c = Gene('1', 1500, 2500)
        genes = [a, b, c]
        """return self.generate_interval_mapping(
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

    def test_generate_gene_mapping_err(self):
        #  _generate_interval_mapping [IntergenicRegion(11:77361962_77361962+)] 1181.39453125 5 30 None 77356962 77366962)
        ir = IntergenicRegion('11', 5000, 5000, STRAND.POS)
        tgt_width = 1000
        d = DiagramSettings()
        d.GENE_MIN_BUFFER = 10
        # (self, canvas, gene, width, height, fill, label='', REFERENCE_GENOME=None)
        _draw_genes(d, self.canvas, [ir], tgt_width, [])

        # _generate_interval_mapping ['Interval(29684391, 29684391)', 'Interval(29663998, 29696515)'] 1181.39453125 5 60 None 29662998 29697515
        # def generate_interval_mapping(cls, input_intervals, target_width, ratio, min_width, buffer_length=None, start=None, end=None, min_inter_width=None)
        itvls = [Interval(29684391, 29684391), Interval(29663998, 29696515)]
        mapping = generate_interval_mapping(itvls, 1181.39, 5, 60, None, 29662998, 29697515)

    def test_split_intervals_into_tracks(self):
        # ----======---------
        # ------======--------
        # -----===============
        t = split_intervals_into_tracks(
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

        d = DiagramSettings()
        breakpoints = [
            Breakpoint('1', 1100, 1200, orient=ORIENT.RIGHT)
        ]
        g = _draw_genes(
            d, self.canvas, [x, y, z], 500, breakpoints,
            {x: d.GENE1_COLOR, y: d.GENE2_COLOR_SELECTED, z: d.GENE2_COLOR})

        # test the class structure
        self.assertEqual(6, len(g.elements))
        self.assertEqual('scaffold', g.elements[0].attribs.get('class', ''))
        for i in range(1, 4):
            self.assertEqual('gene', g.elements[i].attribs.get('class', ''))
        self.assertEqual('mask', g.elements[4].attribs.get('class', ''))
        self.assertEqual('breakpoint', g.elements[5].attribs.get('class', ''))
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
        d = DiagramSettings()
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
        g = _draw_ustranscript(
            d, self.canvas, t, 500,
            colors={t.exons[1]: '#FFFF00'},
            breakpoints=[b]
        )
        self.canvas.add(g)
        self.canvas.saveas('test_draw_ustranscript.svg')
        self.assertEqual(2, len(self.canvas.elements))
        self.assertEqual(3, len(g.elements))
        for el, cls in zip(g.elements[0].elements, ['splicing', 'exon_track', 'protein']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        for el, cls in zip(g.elements[0].elements[1].elements, ['scaffold', 'exon', 'exon', 'exon']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        for el, cls in zip(g.elements[0].elements[2].elements, ['translation', 'domain', 'domain']):
            self.assertEqual(cls, el.attribs.get('class', ''))

        self.assertEqual(
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT
            + 2 * d.PADDING + d.DOMAIN_TRACK_HEIGHT * 2
            + d.TRANSLATION_TRACK_HEIGHT + d.PADDING
            + d.BREAKPOINT_TOP_MARGIN
            + d.BREAKPOINT_BOTTOM_MARGIN,
            g.height)
        self.assertEqual(d1.name, g.labels['D1'])
        self.assertEqual(d2.name, g.labels['D2'])

    def test_dynamic_label_color(self):
        self.assertEqual(HEX_WHITE, dynamic_label_color(HEX_BLACK))
        self.assertEqual(HEX_BLACK, dynamic_label_color(HEX_WHITE))

    def test_draw_legend(self):
        d = DiagramSettings()
        swatches = [
            ('#000000', 'black'),
            ('#FF0000', 'red'),
            ('#0000FF', 'blue'),
            ('#00FF00', 'green'),
            ('#FFFF00', 'yellow')
        ]
        g = draw_legend(d, self.canvas, swatches)
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
        d = DiagramSettings()
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

        canvas, legend = draw(d, ann, ft)
        self.assertEqual(4, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + d.SPLICE_HEIGHT
        canvas.saveas('test_draw_layout_single_transcript.svg')
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_layout_single_genomic(self):
        d = DiagramSettings()
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
        self.assertEqual(t1.exons[0], ft.exon_mapping[ft.exons[0].position])
        self.assertEqual(t2.exons[2], ft.exon_mapping[ft.exons[1].position])
        self.assertEqual(t2.exons[3], ft.exon_mapping[ft.exons[2].position])

        canvas, legend = draw(d, ann, ft)
        self.assertEqual(5, len(canvas.elements))  # defs counts as element

        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT
        self.assertEqual(expected_height, canvas.attribs['height'])
        canvas.saveas('test_draw_layout_single_genomic.svg')

    def test_draw_layout_translocation(self):
        d = DiagramSettings()
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

        canvas, legend = draw(d, ann, ft)
        self.assertEqual(6, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT * 2 + d.PADDING + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT
        self.assertEqual(expected_height, canvas.attribs['height'])

    def test_draw_template(self):
        # def _draw_template(self, canvas, template, target_width, height, labels=None, colors=None):
        d = DiagramSettings()
        canvas = Drawing(size=(1000, 50))
        t = Template(
            '1', 1, 100000,
            bands=[
                BioInterval(None, 1, 8000, 'p1'),
                BioInterval(None, 10000, 15000, 'p2')
            ])
        g = _draw_template(d, canvas, t, 1000)
        canvas.add(g)
        canvas.attribs['height'] = g.height
        canvas = Drawing(size=(1000, 50))

        g = _draw_template(d, canvas, TEMPLATE_METADATA['1'], 1000)
        self.assertEqual(d.BREAKPOINT_TOP_MARGIN + d.BREAKPOINT_BOTTOM_MARGIN + d.TEMPLATE_TRACK_HEIGHT, g.height)
        canvas.add(g)
        canvas.attribs['height'] = g.height
        self.assertEqual(2, len(canvas.elements))

    def test_draw_translocation_with_template(self):
        d = DiagramSettings()
        d1 = Domain('PF0001', [(55, 61), (71, 73)])
        d2 = Domain('PF0002', [(10, 20), (30, 34)])
        g1 = Gene(TEMPLATE_METADATA['1'], 150, 1000, strand=STRAND.POS, aliases=['HUGO2'])
        g2 = Gene(TEMPLATE_METADATA['X'], 5000, 7500, strand=STRAND.NEG, aliases=['HUGO3'])
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
        ann.add_gene(Gene('1', 1500, 1950, strand=STRAND.POS, aliases=['HUGO5']))
        ann.add_gene(Gene('1', 3000, 3980, strand=STRAND.POS))
        ann.add_gene(Gene('1', 3700, 4400, strand=STRAND.NEG))
        # genes 2
        ann.add_gene(Gene('2', 1500, 1950, strand=STRAND.NEG))
        ann.add_gene(Gene('2', 5500, 9000, strand=STRAND.POS))
        ann.add_gene(Gene('2', 3700, 4400, strand=STRAND.NEG))

        reference_genome = {'1': MockSeq(MockString('A')), '2': MockSeq(MockString('A'))}

        ft = FusionTranscript.build(ann, reference_genome)

        canvas, legend = draw(d, ann, ft, show_template=True, templates=TEMPLATE_METADATA)
        canvas.saveas('test_draw_translocation_with_template.svg')
        self.assertEqual(8, len(canvas.elements))  # defs counts as element
        expected_height = d.TOP_MARGIN + d.BOTTOM_MARGIN + \
            d.TRACK_HEIGHT * 2 + d.PADDING  + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.SPLICE_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + \
            d.PADDING + d.TRANSLATION_TRACK_HEIGHT + \
            d.PADDING * 2 + d.DOMAIN_TRACK_HEIGHT * 2 + \
            d.INNER_MARGIN + \
            d.TRACK_HEIGHT + d.BREAKPOINT_BOTTOM_MARGIN + d.BREAKPOINT_TOP_MARGIN + d.SPLICE_HEIGHT + \
            d.TEMPLATE_TRACK_HEIGHT
        self.assertAlmostEqual(expected_height, canvas.attribs['height'])

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
            gene=gene, domains=[Domain('domain1', [(1, 10)]), Domain('domain1', [(4, 10)])],
            is_best_transcript=True)
        build_transcript(
            cds_start=65, cds_end=634,
            exons=[Exon(25403698, 25403863), Exon(25398208, 25398329), Exon(25386753, 25388160)],
            gene=gene, domains=[])
        d = DiagramSettings()
        for i, t in enumerate(gene.transcripts):
            t.name = 'transcript {}'.format(i + 1)
        scatterx = [Interval(x, x + 200) for x in range(gene.start, gene.end + 1, 400)]
        scattery = [Interval(random.uniform(-0.2, 0.2)) for x in scatterx]
        s = ScatterPlot(
            list(zip(scatterx, scattery)),
            'cna',
            ymin=-1,
            ymax=1,
            #hmarkers=[-1, 1],
            yticks=[-1, 0, 1]
        )

        d.GENE_MIN_BUFFER = 0
        canvas = draw_multi_transcript_overlay(d, gene, vmarkers=[marker], plots=[s, s])
        self.assertEqual(2, len(canvas.elements))  # defs counts as element
        canvas.saveas('test_draw_overlay.svg')
