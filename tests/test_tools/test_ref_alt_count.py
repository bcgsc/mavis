import os
import shutil
import tempfile

import pytest
from mavis.annotate.file_io import load_reference_genome
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, SVTYPE
from tools.calculate_ref_alt_counts import RefAltCalculator

from ..util import get_data, glob_exists


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')


def print_file_tree(dirname):
    for root, dirs, files in os.walk(dirname):
        level = root.replace(dirname, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))


@pytest.fixture
def calculator():
    return RefAltCalculator(
        [("TEST", get_data('mock_reads_for_events.sorted.bam'))],
        REFERENCE_GENOME,
        max_event_size=100,
        buffer=20,
    )


@pytest.fixture
def temp_output():
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


class TestFullCalculator:
    def test_calculate_all_counts(self, calculator, temp_output):
        calculator.calculate_all_counts(
            [get_data("mavis_summary_all_mock-A36971_mock-A47933.tab")],
            os.path.join(temp_output, "ref_alt_output.tab"),
        )
        assert glob_exists(temp_output, "ref_alt_output.tab")


class TestRefAltCalulator:
    def test_calculate_count(self, calculator):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 5999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 6003, orient=ORIENT.RIGHT),
            opposing_strands=False,
            event_type=SVTYPE.DEL,
        )
        bpp = calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        assert bpp.data["TEST_ref_count"] == 27
        assert bpp.data["TEST_alt_count"] == 14
        assert bpp.data['TEST_ignored_count'] == 188

    def test_calculate_count2(self, calculator):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 9999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 10030, orient=ORIENT.RIGHT),
            opposing_strands=False,
            event_type=SVTYPE.DEL,
        )
        bpp = calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        assert bpp.data["TEST_ref_count"] == 0
        assert bpp.data["TEST_alt_count"] == 63
        assert bpp.data['TEST_ignored_count'] == 197

    def test_calculate_count3(self, calculator):
        ev1 = BreakpointPair(
            Breakpoint('reference1', 2002, orient=ORIENT.LEFT),
            Breakpoint('reference1', 2003, orient=ORIENT.RIGHT),
            opposing_strands=False,
            event_type=SVTYPE.INS,
            untemplated_seq='TT',
        )
        bpp = calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        assert bpp.data["TEST_ref_count"] == 0
        assert bpp.data["TEST_alt_count"] == 23
        assert bpp.data['TEST_ignored_count'] == 145

    def test_calculate_count4(self, calculator):
        ev1 = BreakpointPair(
            Breakpoint('reference11', 1999, orient=ORIENT.LEFT),
            Breakpoint('reference11', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False,
            event_type=SVTYPE.DEL,
        )
        bpp = calculator.calculate_ref_counts(ev1)
        print(bpp.data)
        assert bpp.data["TEST_ref_count"] == 0
        assert bpp.data["TEST_alt_count"] == 50
        assert bpp.data['TEST_ignored_count'] == 191
