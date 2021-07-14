import unittest

from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.cluster.cluster import merge_breakpoint_pairs, merge_integer_intervals
from mavis.constants import COLUMNS, PROTOCOL, SVTYPE
from mavis.interval import Interval
from mavis.util import read_bpp_from_input_file


from ..util import get_data


FULL_BASE_EVENTS = get_data('mock_sv_events.tsv')
CLUSTERED_EVENTS = get_data('clustering_input.tab')
REF_CHR = 'fake'


class TestFullClustering:
    def test_mocked_events(self):
        # none of the 24 events in the mocked file should cluster together
        # if we change the mock file we may need to update this function
        bpps = []
        for bpp in sorted(
            read_bpp_from_input_file(FULL_BASE_EVENTS), key=lambda x: (x.break1.chr, x.break2.chr)
        ):
            if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
                bpps.append(bpp)
                print(bpp)
        assert len(bpps) == 28
        clusters = merge_breakpoint_pairs(bpps, 10, 10)

        for cluster, input_pairs in sorted(
            clusters.items(), key=lambda x: (x[1][0].break1.chr, x[1][0].break2.chr)
        ):
            print(cluster)
            for ip in input_pairs:
                print('\t', ip)
            assert len(input_pairs) == 1
        assert len(clusters) == len(bpps)

    def test_clustering_events(self):
        # this file contains 2 events that should be clustered and produce a valid bpp
        bpps = []
        for bpp in sorted(
            read_bpp_from_input_file(CLUSTERED_EVENTS), key=lambda x: (x.break1.chr, x.break2.chr)
        ):
            if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
                bpps.append(bpp)
                print(bpp)
        assert len(bpps) == 2
        clusters = merge_breakpoint_pairs(bpps, 200, 25)

        assert len(clusters) == 1

        for cluster, input_pairs in sorted(
            clusters.items(), key=lambda x: (x[1][0].break1.chr, x[1][0].break2.chr)
        ):
            print(cluster)
            for ip in input_pairs:
                print('\t', ip)
            print(cluster.flatten())
            # BPP(Breakpoint(15:67333604L), Breakpoint(15:67333606R), opposing=False)
            assert cluster.break1.orient == 'L'
            assert cluster.break2.orient == 'R'
            assert cluster.break1.chr == '15'
            assert cluster.break2.chr == '15'
            assert cluster.break1.start == 67333604
            assert cluster.break2.start == 67333606
            assert cluster.break1.end == 67333604
            assert cluster.break2.end == 67333606


class TestMergeBreakpointPairs:
    def test_order_is_retained(self):
        # BPP(Breakpoint(1:1925143-1925155R), Breakpoint(1:1925144L), opposing=False)
        # >>  BPP(Breakpoint(1:1925144L), Breakpoint(1:1925144-1925158R), opposing=False)
        # >>  BPP(Breakpoint(1:1925143L), Breakpoint(1:1925143-1925158R), opposing=False)
        pairs = [
            BreakpointPair(
                Breakpoint('2', 1925144, 1925144, 'L'),
                Breakpoint('2', 1925144, 1925158, 'R'),
                event_type='deletion',
                opposing_strands=False,
            ),
            BreakpointPair(
                Breakpoint('2', 1925143, 1925143, 'L'),
                Breakpoint('2', 1925143, 1925158, 'R'),
                event_type='deletion',
                opposing_strands=False,
            ),
        ]
        mapping = merge_breakpoint_pairs(pairs, 100, 25)
        for merge, inputs in mapping.items():
            print(merge)
            print(inputs)
        assert len(mapping) == 1
        merge = list(mapping)[0]
        assert merge.break1.orient == 'L'
        assert merge.break2.orient == 'R'

    def test_merging_identical_large_inputs(self):
        b1 = BreakpointPair(
            Breakpoint(11, 12856838, 12897006, 'L'),
            Breakpoint(11, 12856838, 12897006, 'R'),
            opposing_strands=False,
        )
        b2 = BreakpointPair(
            Breakpoint(11, 12856838, 12897006, 'L'),
            Breakpoint(11, 12856838, 12897006, 'R'),
            opposing_strands=False,
        )
        mapping = merge_breakpoint_pairs([b1, b2], 100, 25, verbose=True)
        assert len(mapping) == 1
        merge = list(mapping)[0]
        assert len(mapping[merge]) == 2
        assert merge.break1.orient == 'L'
        assert merge.break2.orient == 'R'
        assert merge.break1.chr == '11'
        assert merge.break2.chr == '11'
        assert merge.break1.start == 12856838
        assert merge.break2.start == 12856840  # putative indel will be shifted
        assert merge.break1.end == 12897006
        assert merge.break2.end == 12897006

    def test_events_separate(self):
        bpps = [
            BreakpointPair(
                Breakpoint('4', 157002032, 157002046, orient='L'),
                Breakpoint('4', 157002343, 157002343, orient='R'),
                event_type=SVTYPE.DEL,
                untemplated_seq='',
                protocol='genome',
                tracking_id='manta-MantaDEL:55718:0:0:1:0:0',
            ),
            BreakpointPair(
                Breakpoint('4', 156935061, orient='L'),
                Breakpoint('4', 156935245, orient='R'),
                event_type=SVTYPE.DEL,
                untemplated_seq='',
                protocol='genome',
                tracking_id='delly-DEL00011007',
            ),
            BreakpointPair(
                Breakpoint('4', 157002046, orient='L'),
                Breakpoint('4', 157002358, orient='R'),
                event_type=SVTYPE.DEL,
                untemplated_seq='',
                protocol='genome',
                tracking_id='delly-DEL00011008',
            ),
        ]
        mapping = merge_breakpoint_pairs(bpps, 100, 25, verbose=True)
        assert len(mapping) == 2


class TestMergeIntervals:
    def test_merge_even_length(self):
        i1 = Interval(1001, 1002)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        assert result == i1

    def test_merge_odd_length(self):
        i1 = Interval(1001, 1003)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        assert result == i1

    def test_merge_large_length(self):
        i1 = Interval(1001, 5003)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        assert result == i1

        i1 = Interval(12856838, 12897006)
        result = merge_integer_intervals(i1, i1, weight_adjustment=25)
        assert result == i1


if __name__ == '__main__':
    unittest.main()
