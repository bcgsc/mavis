import os
import unittest

from mavis.constants import COLUMNS, ORIENT, STRAND
from mavis.error import NotSpecifiedError
from mavis.util import cast, ENV_VAR_PREFIX, get_env_variable, MavisNamespace, WeakMavisNamespace, read_bpp_from_input_file, get_connected_components

from .mock import Mock


class MockFileHandle(Mock):
    def __init__(self, lines):
        Mock.__init__(self, lines=lines)

    def readlines(self):
        return self.lines


class TestGetConnectedComponents(unittest.TestCase):

    def test_no_nodes(self):
        self.assertEqual([], get_connected_components({}))

    def test_no_connections(self):
        graph = {1: {}, 2: {}, 3: {}}
        components = get_connected_components(graph)
        self.assertEqual(3, len(components))

    def test_fully_connected(self):
        graph = {1: {2, 3, 1}, 2: {1, 2, 2}, 3: {3, 2}}
        components = get_connected_components(graph)
        self.assertEqual(1, len(components))
        self.assertEqual([{1, 2, 3}], components)

    def test_multiple_components(self):
        graph = {1: {2}, 2: {3}, 3: {4}, 6: {7, 8}}
        components = get_connected_components(graph)
        self.assertEqual(2, len(components))
        self.assertEqual({1, 2, 3, 4}, components[0])
        self.assertEqual({6, 7, 8}, components[1])


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
        print(self.namespace)
        self.namespace.thing = 2
        print(self.namespace)
        self.namespace.otherthing = 3
        print(self.namespace)
        self.assertEqual({'thing': 2, 'otherthing': 3}, self.namespace._members)
        self.assertEqual([('otherthing', 3), ('thing', 2)], list(sorted(self.namespace.items())))


class TestWeakMavisNamespace(unittest.TestCase):
    def setUp(self):
        self.namespace = WeakMavisNamespace(a=1, b=2, c=3)
        print(self.namespace._members)
        for v in ['a', 'b', 'c']:
            v = ENV_VAR_PREFIX + v.upper()
            if v in os.environ:
                del os.environ[v]

    def test_no_env_set(self):
        self.assertEqual(1, self.namespace.a)
        self.assertEqual(1, self.namespace['a'])

    def test_env_overrides_default(self):
        os.environ['MAVIS_A'] = '5'
        env_name = self.namespace.get_env_name('a')
        self.assertEqual('MAVIS_A', env_name)
        self.assertEqual('5', os.environ[env_name])
        self.assertTrue(self.namespace.is_env_overwritable('a'))
        self.assertEqual(5, self.namespace.a)
        self.assertEqual(1, self.namespace._members['a'])
        self.assertEqual(5, self.namespace['a'])

    def test_error_on_invalid_attr(self):
        with self.assertRaises(AttributeError):
            self.namespace.other

    def test_iterate_keys(self):
        self.assertEqual(['a', 'b', 'c'], list(self.namespace.keys()))
        self.assertEqual(['a', 'b', 'c'], [k for k in self.namespace])


class TestReadBreakpointPairsFromFile(unittest.TestCase):

    def build_filehandle(self, row):
        header = [c for c in row]
        line = [row[c] for c in header]
        lines = ['\t'.join(header), '\t'.join([str(v) for v in line])]
        return MockFileHandle(lines)

    def test_break1_strand_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: False
        })
        with self.assertRaises(NotSpecifiedError):
            bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=True)
            for b in bpps:
                print(b)

    def test_stranded_no_expand_error(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)

    def test_break2_strand_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: False
        })

        with self.assertRaises(NotSpecifiedError) as err:
            print(err)
            bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=False)

        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)

    def test_stranded_expand_strands_and_orient(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=False)
        self.assertEqual(2, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break2.strand)

    def test_expand_strands_and_orient(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

    def test_stranded_expand_strands_not_orient(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=True)
        self.assertEqual(2, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.POS, bpps[0].break2.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[1].break2.strand)

    def test_expand_orient_not_strand(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.NS,
            COLUMNS.break1_orientation: ORIENT.LEFT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

    def test_break1_orient_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.NS,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)

    def test_break2_orient_ns(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.NS,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.POS,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: False
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)
        raise unittest.SkipTest('TODO')

    def test_both_break_orient_ns(self):
        raise unittest.SkipTest('TODO')

    def test_base_case(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.RIGHT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NEG,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: True
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(ORIENT.RIGHT, bpps[0].break1.orient)
        self.assertEqual(True, bpps[0].opposing_strands)

    def test_unstranded_with_strand_calls(self):
        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.RIGHT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NEG,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: False,
            COLUMNS.opposing_strands: True
        })
        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, expand_strand=False, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.NS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NS, bpps[0].break2.strand)

        fh = self.build_filehandle({
            COLUMNS.break1_chromosome: '1',
            COLUMNS.break1_position_start: 1,
            COLUMNS.break1_position_end: 1,
            COLUMNS.break1_strand: STRAND.POS,
            COLUMNS.break1_orientation: ORIENT.RIGHT,
            COLUMNS.break2_chromosome: '1',
            COLUMNS.break2_position_start: 10,
            COLUMNS.break2_position_end: 10,
            COLUMNS.break2_strand: STRAND.NEG,
            COLUMNS.break2_orientation: ORIENT.RIGHT,
            COLUMNS.stranded: True,
            COLUMNS.opposing_strands: True
        })

        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=False)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[0].break2.strand)

        bpps = read_bpp_from_input_file(fh, expand_strand=True, expand_orient=True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(STRAND.POS, bpps[0].break1.strand)
        self.assertEqual(STRAND.NEG, bpps[0].break2.strand)
