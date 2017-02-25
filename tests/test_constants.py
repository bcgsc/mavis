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

    def test_reverse_complement(self):
        self.assertEqual('ATCG', reverse_complement('CGAT'))
        self.assertEqual('', reverse_complement(''))

    def test_translate(self):
        s = 'ATG' 'AAT' 'TCT' 'GGA' 'TGA'
        t = translate(s, 0)
        self.assertEqual('MNSG*', t)  # ATG AAT TCT GGA TGA
        t = translate(s, 1)
        self.assertEqual('*ILD', t)  # A TGA ATT CTG GAT GA
        t = translate(s, 2)
        self.assertEqual('EFWM', t)  # AT GAA TTC TGG ATG A

    def test_sort_columns(self):
        temp = ['NEW', 'NEW2', COLUMNS.break1_seq, COLUMNS.break2_seq, COLUMNS.break1_chromosome]
        self.assertEqual(
            [COLUMNS.break1_chromosome, COLUMNS.break1_seq, COLUMNS.break2_seq, 'NEW', 'NEW2'],
            sort_columns(temp))

    def test_column_matches_column_name(self):
        self.assertEqual(COLUMNS.library, COLUMNS.library)
        s = set([COLUMNS.library, COLUMNS.library])
        self.assertEqual(1, len(s))
