import unittest
from argparse import ArgumentTypeError

from mavis.config import float_fraction, nameable_string


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


class TestNoReservedChars(unittest.TestCase):

    def test_semicolon_error(self):
        with self.assertRaises(TypeError):
            nameable_string('thing;thing')

    def test_comma_error(self):
        with self.assertRaises(TypeError):
            nameable_string('thing,thing')

    def test_underscore_error(self):
        with self.assertRaises(TypeError):
            nameable_string('thing_thing')

    def test_space_error(self):
        with self.assertRaises(TypeError):
            nameable_string(' ')

        with self.assertRaises(TypeError):
            nameable_string('thing thing')

    def test_ok(self):
        lib = 'libName'
        self.assertEqual('libName', nameable_string(lib))

    def test_number_start_error(self):
        with self.assertRaises(TypeError):
            nameable_string('1thing')

        with self.assertRaises(TypeError):
            nameable_string('1')

    def test_empty_error(self):
        with self.assertRaises(TypeError):
            nameable_string('')

    def test_none_error(self):
        with self.assertRaises(TypeError):
            nameable_string('none')

        with self.assertRaises(TypeError):
            nameable_string(None)
