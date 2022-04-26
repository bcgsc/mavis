import os
import pandas as pd
import pytest
from tools.convert_dgv import (
    convert_to_dictionary,
)


@pytest.mark.parametrize(
    'filename,expected_file',
    [
        ['dgv_test.tab', 'dgv_test_expected.tab'],
    ],
)
def test_dgv_examples(filename, expected_file):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')

    expected = pd.read_csv(os.path.join(data_dir, expected_file), sep='\t', header=0)
    input_dict = convert_to_dictionary(os.path.join(data_dir, filename))
    assert len(input_dict) == expected.shape[0]
