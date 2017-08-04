from mavis.cluster.cluster import merge_breakpoint_pairs
from mavis.breakpoint import read_bpp_from_input_file
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


if __name__ == "__main__":
    unittest.main()
