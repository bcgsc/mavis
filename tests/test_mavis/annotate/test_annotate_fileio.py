import json

import pytest
from mavis.annotate.file_io import load_annotations


@pytest.mark.parametrize(
    'annotations,error_message_include',
    [
        [{'genes': []}, "schema['properties']['genes']"],
        [
            {'genes': [{'start': '1'}]},
            "schema['properties']['genes']['items']['properties']['start']",
        ],
    ],
)
def test_min_genes_error(annotations, error_message_include, tmp_path):
    filename = tmp_path / "annotations.json"
    filename.write_text(json.dumps(annotations))
    with pytest.raises(AssertionError) as exc:
        load_annotations(str(filename))
    assert error_message_include in str(exc.value)
