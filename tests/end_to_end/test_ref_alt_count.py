import tempfile
import unittest
import os
import shutil

from mavis.annotate.file_io import load_reference_genome
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, SVTYPE
from tools.calculate_ref_alt_counts import RefAltCalculator

from ..util import get_data
from . import glob_exists


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if 'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT' != REFERENCE_GENOME['fake'].seq[0:50].upper():
        raise AssertionError('fake genome file does not have the expected contents')


def print_file_tree(dirname):
    for root, dirs, files in os.walk(dirname):
        level = root.replace(dirname, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))


class TestFullCalculator(unittest.TestCase):
    def setUp(self):
        # create the temp output directory to store file outputs
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

        self.calculator = RefAltCalculator([("TEST", get_data('mock_reads_for_events.sorted.bam'))], REFERENCE_GENOME,
                                           max_event_size=100, buffer=20)

    def test_calculate_all_counts(self):
        self.calculator.calculate_all_counts([get_data("mavis_summary_all_mock-A36971_mock-A47933.tab")],
                                             os.path.join(self.temp_output, "ref_alt_output.tab"))
        self.assertTrue(glob_exists(self.temp_output, "ref_alt_output.tab"))

    def tearDown(self):
        # remove the temp directory and outputs
        print_file_tree(self.temp_output)
        shutil.rmtree(self.temp_output)


class TestRefAltCalulator(unittest.TestCase):

    def setUp(self):
        self.calculator = RefAltCalculator([("TEST", get_data('mock_reads_for_events.sorted.bam'))], REFERENCE_GENOME,
                                           max_event_size=100, buffer=20)

    def test_calculate_count(self):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 5999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 6003, orient=ORIENT.RIGHT),
            opposing_strands=False, event_type=SVTYPE.DEL
        )
        bpp = self.calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        self.assertEqual(27, bpp.data["TEST_ref_count"])
        self.assertEqual(14, bpp.data["TEST_alt_count"])
        self.assertEqual(188, bpp.data['TEST_ignored_count'])

    def test_calculate_count2(self):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 9999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 10030, orient=ORIENT.RIGHT),
            opposing_strands=False, event_type=SVTYPE.DEL
        )
        bpp = self.calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        self.assertEqual(0, bpp.data["TEST_ref_count"])
        self.assertEqual(63, bpp.data["TEST_alt_count"])
        self.assertEqual(197, bpp.data['TEST_ignored_count'])

    def test_calculate_count3(self):
        ev1 = BreakpointPair(
            Breakpoint('reference1', 2002, orient=ORIENT.LEFT),
            Breakpoint('reference1', 2003, orient=ORIENT.RIGHT),
            opposing_strands=False, event_type=SVTYPE.INS, untemplated_seq='TT'
        )
        bpp = self.calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        self.assertEqual(0, bpp.data["TEST_ref_count"])
        self.assertEqual(23, bpp.data["TEST_alt_count"])
        self.assertEqual(145, bpp.data['TEST_ignored_count'])

    def test_calculate_count4(self):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 1999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False, event_type=SVTYPE.DEL
        )
        bpp = self.calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        self.assertEqual(0, bpp.data["TEST_ref_count"])
        self.assertEqual(50, bpp.data["TEST_alt_count"])
        self.assertEqual(191, bpp.data['TEST_ignored_count'])
