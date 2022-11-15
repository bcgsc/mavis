import os
import pytest
from unittest.mock import patch
from tools.convert_dgv import main as convert_dgv_main
import sys


@pytest.mark.parametrize(
    "filename,expected_file",
    [
        ["dgv_test.tab", "dgv_test_expected.tab"],
    ],
)
def test_dgv_examples(tmp_path, filename, expected_file):
    data_dir = os.path.join(os.path.dirname(__file__), "data")

    output_path = str(tmp_path / "tmp_data.tab")
    args = [
        "python",
        "--input",
        os.path.join(data_dir, filename),
        "--output",
        output_path,
    ]

    with patch.object(convert_dgv_main, "main", create=True):

        with patch.object(sys, "argv", args):
            convert_dgv_main()

    with open(os.path.join(data_dir, expected_file), 'r') as fh:
        expected = fh.read().strip()

    with open(output_path, 'r') as fh:
        observed = fh.read().strip()

    assert expected == observed
