from mavis.cluster.cluster import cluster_breakpoint_pairs
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
        for bpp in read_bpp_from_input_file(FULL_BASE_EVENTS):
            if bpp.data[COLUMNS.protocol] == PROTOCOL.GENOME:
                bpps.append(bpp)
        self.assertEqual(28, len(bpps))
        clusters = cluster_breakpoint_pairs(bpps, 10, 10)
        self.assertEqual(len(bpps), len(clusters))
        for cluster, input_pairs in clusters.items():
            self.assertEqual(1, len(input_pairs))


if __name__ == "__main__":
    unittest.main()
