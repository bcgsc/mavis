import unittest
import re
from structural_variant.draw import Diagram
from structural_variant.annotate import Gene, Transcript
from svgwrite import Drawing
from structural_variant.constants import STRAND
from structural_variant.interval import Interval
from structural_variant.breakpoint import Breakpoint


class TestDraw(unittest.TestCase):
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

    def test_draw_gene_subdiagram(self):
        canvas = Drawing(height=100, width=1000)
        d = Diagram(
            PADDING=5,
            GENE1_COLOR='#325556',
        )
        genes = [
            Gene('1', 1000, 2000, strand=STRAND.POS),
            Gene('1', 5000, 7000, strand=STRAND.NEG),
            Gene('1', 1500, 2500, strand=STRAND.POS)
        ]
        breakpoints = [
            Breakpoint('1', 1100, 1200)
        ]
        g = d.draw_gene_subdiagram(canvas, 500, genes, breakpoints)
        canvas.add(g)
        #canvas.saveas('test.svg')

    def test_draw_transcript(self):
        canvas = Drawing(height=100, width=1000)
        d = Diagram(
            PADDING=5,
        )
        t = Transcript(gene=None, cds_start=50, cds_end=249, exons=[(1, 99), (200, 299), (400, 499)], strand=STRAND.POS)

        g = d.draw_transcript(canvas, 500, t, exon_color='#325556', utr_color='#FFFF00')
        print(g)
        canvas.add(g)
        print(canvas.tostring())
        canvas.saveas('test.svg')
        self.assertFalse(True)
