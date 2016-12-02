import unittest
from structural_variant.draw import Diagram
from structural_variant.annotate import Gene

class TestDraw(unittest.TestCase):
    def test__draw_gene_track(self):
        # don't actually need to use genes here
        d = Diagram()
        genes = [Gene('1', 1000, 2000), Gene('1', 5000, 7000), Gene('1', 1500, 2500)]
        d._generate_gene_mapping(100, genes)
