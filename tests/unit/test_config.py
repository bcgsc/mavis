import unittest
from argparse import ArgumentTypeError

from mavis.config import float_fraction


class TestFloatFraction(unittest.TestCase):

    def test_bad_string(self):
        with self.assertRaises(ArgumentTypeError):
            float_fraction('a')

    def test_float_too_big(self):
        with self.assertRaises(ArgumentTypeError):
            float_fraction('1.1')

    def test_float_negative_error(self):
        with self.assertRaises(ArgumentTypeError):
            float_fraction('-0.1')

    def test_zero_ok(self):
        self.assertEqual(0, float_fraction('0'))

    def test_one_ok(self):
        self.assertEqual(1, float_fraction('1'))
