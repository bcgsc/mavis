import unittest

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.cluster.cluster import merge_breakpoint_pairs, merge_integer_intervals
from mavis.constants import COLUMNS, PROTOCOL, SVTYPE
from mavis.interval import Interval
from mavis.util import read_bpp_from_input_file


from ..util import get_data


FULL_BASE_EVENTS = get_data('mock_sv_events.tsv')
REF_CHR = 'fake'


class TestFullClustering(unittest.TestCase):
    def test_mocked_events(self):
        # none of the 24 events in the mocked file should cluster together
        # if we change the mock file we may need to update this function
        bpps = []
        for bpp in sorted(read_bpp_from_input_file(FULL_BASE_EVENTS), key=lambda x: (x.break1.chr, x.break2.chr)):
            if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
                bpps.append(bpp)
                print(bpp)
        self.assertEqual(28, len(bpps))
        clusters = merge_breakpoint_pairs(bpps, 10, 10)

        for cluster, input_pairs in sorted(clusters.items(), key=lambda x: (x[1][0].break1.chr, x[1][0].break2.chr)):
            print(cluster)
            for ip in input_pairs:
                print('\t', ip)
            self.assertEqual(1, len(input_pairs))
        self.assertEqual(len(bpps), len(clusters))


class TestMergeBreakpointPairs(unittest.TestCase):

    def test_order_is_retained(self):
        # BPP(Breakpoint(1:1925143-1925155R), Breakpoint(1:1925144L), opposing=False)
        # >>  BPP(Breakpoint(1:1925144L), Breakpoint(1:1925144-1925158R), opposing=False)
        # >>  BPP(Breakpoint(1:1925143L), Breakpoint(1:1925143-1925158R), opposing=False)
        pairs = [
            BreakpointPair(
                Breakpoint('2', 1925144, 1925144, 'L'),
                Breakpoint('2', 1925144, 1925158, 'R'),
                event_type='deletion',
                opposing_strands=False),
            BreakpointPair(
                Breakpoint('2', 1925143, 1925143, 'L'),
                Breakpoint('2', 1925143, 1925158, 'R'),
                event_type='deletion',
                opposing_strands=False)
        ]
        mapping = merge_breakpoint_pairs(pairs, 100, 25)
        for merge, inputs in mapping.items():
            print(merge)
            print(inputs)
        self.assertEqual(1, len(mapping))
        merge = list(mapping)[0]
        self.assertEqual('L', merge.break1.orient)
        self.assertEqual('R', merge.break2.orient)

    def test_merging_identical_large_inputs(self):
        b1 = BreakpointPair(Breakpoint(11, 12856838, 12897006, 'L'), Breakpoint(11, 12856838, 12897006, 'R'), opposing_strands=False)
        b2 = BreakpointPair(Breakpoint(11, 12856838, 12897006, 'L'), Breakpoint(11, 12856838, 12897006, 'R'), opposing_strands=False)
        mapping = merge_breakpoint_pairs([b1, b2], 100, 25, verbose=True)
        self.assertEqual(1, len(mapping))
        merge = list(mapping)[0]
        self.assertEqual(2, len(mapping[merge]))
        self.assertEqual('L', merge.break1.orient)
        self.assertEqual('R', merge.break2.orient)
        self.assertEqual('11', merge.break1.chr)
        self.assertEqual('11', merge.break2.chr)
        self.assertEqual(12856838, merge.break1.start)
        self.assertEqual(12856838, merge.break2.start)
        self.assertEqual(12897006, merge.break1.end)
        self.assertEqual(12897006, merge.break2.end)

    def test_events_separate(self):
        bpps = [
            BreakpointPair(Breakpoint('4', 157002032, 157002046, orient='L'), Breakpoint('4', 157002343, 157002343, orient='R'), event_type=SVTYPE.DEL, untemplated_seq='', protocol='genome', tracking_id='manta-MantaDEL:55718:0:0:1:0:0'),
            BreakpointPair(Breakpoint('4', 156935061, orient='L'), Breakpoint('4', 156935245, orient='R'), event_type=SVTYPE.DEL, untemplated_seq='', protocol='genome', tracking_id='delly-DEL00011007'),
            BreakpointPair(Breakpoint('4', 157002046, orient='L'), Breakpoint('4', 157002358, orient='R'), event_type=SVTYPE.DEL, untemplated_seq='', protocol='genome', tracking_id='delly-DEL00011008')
        ]
        mapping = merge_breakpoint_pairs(bpps, 100, 25, verbose=True)
        self.assertEqual(2, len(mapping))


class TestMergeIntervals(unittest.TestCase):

    def test_merge_even_length(self):
        i1 = Interval(1001, 1002)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        self.assertEqual(i1, result)

    def test_merge_odd_length(self):
        i1 = Interval(1001, 1003)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        self.assertEqual(i1, result)

    def test_merge_large_length(self):
        i1 = Interval(1001, 5003)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        self.assertEqual(i1, result)

        i1 = Interval(12856838, 12897006)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        self.assertEqual(i1, result)


if __name__ == '__main__':
    unittest.main()
