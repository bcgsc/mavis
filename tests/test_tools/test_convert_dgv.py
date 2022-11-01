import os
import pandas as pd
import pytest
from unittest.mock import Mock, patch
from tools.convert_dgv import main as convert_dgv_main
import sys


@pytest.mark.parametrize(
    'filename,expected_file',
    [
        ['dgv_test.tab', 'dgv_test_expected.tab'],
    ],
)

def test_dgv_examples(tmp_path, filename, expected_file):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')

    output_path = str(tmp_path / "tmp_data.tab")
    expected = pd.read_csv(os.path.join(data_dir, expected_file), sep='\t', header=0)
    args = [
        'python',
        '--input',
        os.path.join(data_dir, filename),
        '--output',
        output_path
    ]

    with patch.object(convert_dgv_main, 'main', Mock(), create = True):

        with patch.object(sys, 'argv', args) as m:
            convert_dgv_main()

    result = pd.read_csv(output_path, sep="\t")    
    assert result.shape[0] == len(expected)
