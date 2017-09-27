from mavis.util import ChrListString, cast, get_env_variable, MavisNamespace, WeakMavisNamespace, ENV_VAR_PREFIX
import unittest
import os


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


class TestCast(unittest.TestCase):
    def test_float(self):
        self.assertEqual(type(1.0), type(cast('1', float)))
        self.assertNotEqual(type(1.0), type(cast('1', int)))

    def test_boolean(self):
        self.assertEqual(type(False), type(cast('f', bool)))
        self.assertEqual(type(False), type(cast('false', bool)))
        self.assertFalse(cast('f', bool))
        self.assertFalse(cast('false', bool))
        self.assertFalse(cast('0', bool))
        self.assertFalse(cast('F', bool))


class TestGetEnvVariable(unittest.TestCase):
    def setUp(self):
        if 'MAVIS_TEST_ENV' in os.environ:
            del os.environ['MAVIS_TEST_ENV']

    def test_not_set(self):
        self.assertEqual(1, get_env_variable('test_env', 1))

    def test_needs_casting(self):
        os.environ['MAVIS_TEST_ENV'] = '15'
        self.assertEqual(15, get_env_variable('test_env', 1))


class TestMavisNamespace(unittest.TestCase):
    def setUp(self):
        self.namespace = MavisNamespace()

    def test_item_getter(self):
        self.namespace.thing = 2
        self.assertEqual(2, self.namespace['thing'])
        self.assertEqual(2, self.namespace.thing)

    def test_items(self):
        self.namespace.thing = 2
        self.namespace.otherthing = 3
        self.assertEqual([('otherthing', 3), ('thing', 2)], list(sorted(self.namespace.items())))


class TestWeakMavisNamespace(unittest.TestCase):
    def setUp(self):
        self.namespace = WeakMavisNamespace(a=1, b=2, c=3)
        for v in ['a', 'b', 'c']:
            v = ENV_VAR_PREFIX + v.upper()
            if v in os.environ:
                del os.environ[v]

    def test_no_env_set(self):
        self.assertEqual(1, self.namespace.a)
        self.assertEqual(1, self.namespace['a'])

    def test_env_overrides_default(self):
        os.environ['MAVIS_A'] = '5'
        self.assertEqual(5, self.namespace.a)
        self.assertEqual(1, self.namespace.__dict__['a'])
        self.assertEqual(5, self.namespace['a'])

    def test_error_on_invalid_attr(self):
        with self.assertRaises(AttributeError):
            self.namespace.other
