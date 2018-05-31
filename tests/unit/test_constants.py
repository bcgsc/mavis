import unittest
from mavis.constants import COLUMNS, MavisNamespace, ORIENT, reverse_complement, sort_columns, STRAND, translate


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
        seq = 'ATG' 'AAT' 'TCT' 'GGA' 'TGA'
        translated_seq = translate(seq, 0)
        self.assertEqual('MNSG*', translated_seq)  # ATG AAT TCT GGA TGA
        translated_seq = translate(seq, 1)
        self.assertEqual('*ILD', translated_seq)  # A TGA ATT CTG GAT GA
        translated_seq = translate(seq, 2)
        self.assertEqual('EFWM', translated_seq)  # AT GAA TTC TGG ATG A

    def test_sort_columns(self):
        temp = ['NEW', 'NEW2', COLUMNS.break1_seq, COLUMNS.break2_seq, COLUMNS.break1_chromosome]
        self.assertEqual(
            [COLUMNS.break1_chromosome, COLUMNS.break1_seq, COLUMNS.break2_seq, 'NEW', 'NEW2'],
            sort_columns(temp))

    def test_column_matches_column_name(self):
        self.assertEqual(COLUMNS.library, COLUMNS.library)
        s = set([COLUMNS.library, COLUMNS.library])
        self.assertEqual(1, len(s))


class TestMavisNamespace(unittest.TestCase):

    def setUp(self):
        self.namespace = MavisNamespace(a=1, b=2, c=3)

    def test_get_item(self):
        self.assertEqual(1, self.namespace['a'])
        self.assertEqual(1, self.namespace.a)
        self.assertEqual(1, self.namespace.get('a', None))

    def test_to_dict(self):
        self.assertEqual({'a': 1, 'b': 2, 'c': 3}, self.namespace.to_dict())

    def test_get_with_default(self):
        self.assertEqual(4, self.namespace.get('d', 4))

    def test_get_without_default_errors(self):
        self.assertEqual(None, self.namespace.get('d', None))

    def test_error_on_undefined(self):
        with self.assertRaises(KeyError):
            self.namespace.define('a')

    def test_infered_typing(self):
        self.assertEqual(int, self.namespace.type('a'))

    def test_keys(self):
        self.assertEqual(['a', 'b', 'c'], self.namespace.keys())

    def test_add(self):
        self.namespace.add('d', 4, defn='this is the letter d', cast_type=float)
        self.assertEqual(float, self.namespace.type('d'))
        self.assertEqual('this is the letter d', self.namespace.define('d'))
        self.assertEqual(4, self.namespace.d)

    def test_add_infer_type(self):
        self.namespace.add('d', 4, defn='this is the letter d')
        self.assertEqual(int, self.namespace.type('d'))
        self.assertEqual('this is the letter d', self.namespace.define('d'))
        self.assertEqual(4, self.namespace.d)

    def test_error_on_enforce_bad_value(self):
        with self.assertRaises(KeyError):
            self.namespace.enforce(5)

    def test_reverse(self):
        self.assertEqual('a', self.namespace.reverse(1))

    def test_reverse_nonunique_error(self):
        self.namespace['d'] = 1
        with self.assertRaises(KeyError):
            self.namespace.reverse(1)

    def test_reverse_bad_value_error(self):
        with self.assertRaises(KeyError):
            self.namespace.reverse(5)

    def test_get_argument_error(self):
        with self.assertRaises(TypeError):
            self.namespace.get('a', 1, 1)
        with self.assertRaises(AttributeError):
            self.namespace.get('d')

    def test_iterating(self):
        for act, exp in zip(self.namespace, ['a', 'b', 'c']):
            self.assertEqual(exp, act)
