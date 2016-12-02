import unittest
from structural_variant.constants import *


class TestConstants(unittest.TestCase):

    def test_strand_compare(self):
        self.assertTrue(STRAND.compare(STRAND.NS, STRAND.POS))
        self.assertTrue(STRAND.compare(STRAND.NS, STRAND.NEG))
        self.assertTrue(STRAND.compare(STRAND.POS, STRAND.POS))
        self.assertTrue(STRAND.compare(STRAND.NEG, STRAND.NEG))
        self.assertFalse(STRAND.compare(STRAND.POS, STRAND.NEG))
        self.assertFalse(STRAND.compare(STRAND.NEG, STRAND.POS))

    def test_orient_compare(self):
        self.assertTrue(ORIENT.compare(ORIENT.NS, ORIENT.RIGHT))
        self.assertTrue(ORIENT.compare(ORIENT.NS, ORIENT.LEFT))
        self.assertTrue(ORIENT.compare(ORIENT.RIGHT, ORIENT.RIGHT))
        self.assertTrue(ORIENT.compare(ORIENT.LEFT, ORIENT.LEFT))
        self.assertFalse(ORIENT.compare(ORIENT.RIGHT, ORIENT.LEFT))
        self.assertFalse(ORIENT.compare(ORIENT.LEFT, ORIENT.RIGHT))
