from mavis.cluster.cluster import merge_breakpoint_pairs
from mavis.breakpoint import read_bpp_from_input_file, Breakpoint, BreakpointPair
from mavis.constants import PROTOCOL, COLUMNS
from . import FULL_BASE_EVENTS

import unittest

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

if __name__ == "__main__":
    unittest.main()
