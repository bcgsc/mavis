import unittest
from unittest.mock import mock_open, patch
import configparser

from mavis.config import MavisConfig


STUB = """
[reference]
template_metadata = tests/data/cytoBand.txt
annotations = tests/data/mock_annotations.json
masking = tests/data/mock_masking.tab
reference_genome = tests/data/mock_reference_genome.fa
aligner_reference = tests/data/mock_reference_genome.2bit
dgv_annotation = tests/data/mock_dgv_annotation.txt

[mock-A36971]
read_length = 150
median_fragment_size = 400
stdev_fragment_size = 97
bam_file = tests/data/mock_reads_for_events.sorted.bam
protocol = genome
inputs = mock_converted
strand_specific = False
disease_status=diseased

[mock-A47933]
read_length = 75
median_fragment_size = 188
stdev_fragment_size = 50
bam_file = tests/data/mock_trans_reads_for_events.sorted.bam
protocol = transcriptome
inputs = mock_converted
strand_specific = True
disease_status=diseased

[convert]
assume_no_untemplated = True
# addfile twice to check this notation is ok (will collapse them anyway)
mock_converted = convert_tool_output
    tests/data/mock_sv_events.tsv
    tests/data/mock_sv_events.tsv
    mavis
    False
"""


class TestConfig(unittest.TestCase):

    def mock_config(self, content=""):
        with patch('configparser.ConfigParser.read', configparser.ConfigParser.read_string), patch('os.path.isfile') as isfile, patch('os.path.exists') as exists:
            isfile.return_value = True
            exists.return_value = True
            return MavisConfig.read(content)

    def test_error_in_schedule(self):
        with self.assertRaises(TypeError):
            content = STUB + '\n[schedule]\nmail_type=\n'
            print(content)
            self.mock_config(content)

    def test_ok(self):
        self.mock_config(STUB)
