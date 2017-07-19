from mavis.util import ChrListString
import unittest


class TestChrListString(unittest.TestCase):
    def test_cast_from_list(self):
        c = ChrListString(['1', '2', '3'])
        self.assertTrue('1' in c)
        self.assertTrue('2' in c)
        self.assertTrue('3' in c)
        self.assertFalse('4' in c)

    def test_cast_from_string(self):
        c = ChrListString('1;2;3;4')
        self.assertEqual(['1', '2', '3', '4'], list(c))
        self.assertTrue('1' in c)
        self.assertFalse('x' in c)
