import tempfile
import shutil
import unittest
import os
import subprocess
from mavis.constants import STRAND
from mavis.breakpoint import read_bpp_from_input_file


data_prefix = os.path.join(os.path.dirname(__file__), 'data')
temp_output = None
main_run_script = os.path.join(os.path.dirname(__file__), './../../bin/mavis_run.py')


def setUpModule():
    global temp_output
    # create the temp output directory to store file outputs
    temp_output = tempfile.mkdtemp()


class TestPairing(unittest.TestCase):
    def test_pairing(self):
        command = 'python {2} pairing -n {0}/pairing_annotations.tab -f {0}/pairing_sequences.fa -o {1} --annotations' \
            ' {0}/pairing_reference_annotations_file.tab'.format(data_prefix, temp_output, main_run_script)
        print(command)
        subprocess.check_output(command, shell=True)
        # make sure the output file exists
        output = os.path.join(temp_output, 'mavis_paired_A36971_A36973.tab')
        self.assertTrue(os.path.exists(output))
        # check that the expected pairings are present
        bpps = read_bpp_from_input_file(output, explicit_strand=True, expand_ns=False)
        self.assertEqual(6, len(bpps))
        for bpp in bpps:
            self.assertTrue(bpp.break1.strand != STRAND.NS)
            self.assertTrue(bpp.break2.strand != STRAND.NS)


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(temp_output)


if __name__ == "__main__":
    unittest.main()
