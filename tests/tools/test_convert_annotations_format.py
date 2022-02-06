import json
import os

import pytest

from tools.convert_annotations_format import (
    convert_gff2_to_mavis,
    convert_gff3_to_mavis,
    convert_mavis_json_2to3,
)

CONVERTERS = {
    'gff3': convert_gff3_to_mavis,
    'gtf': convert_gff2_to_mavis,
    'v2-json': convert_mavis_json_2to3,
}


def sort_elements(data):
    """
    Sort lists of exons, domains, genes, etc by position and name to facilitate comparison
    """
    if not isinstance(data, dict):
        if isinstance(data, list):
            items = [sort_elements(e) for e in data]

            if all(isinstance(elem, dict) for elem in data):
                return sorted(
                    items, key=lambda elem: (elem.get('start'), elem.get('end'), elem.get('name'))
                )
            return items
        else:
            return data

    for key, value in data.items():
        data[key] = sort_elements(value)
    return data


@pytest.mark.parametrize(
    'filename,expected_file,input_type',
    [
        ['K02718.1.gff3', 'K02718.1.gff3.json', 'gff3'],
        ['K02718.1.gtf', 'K02718.1.gtf.json', 'gtf'],
        ['Homo_sapiens.GRCh38.kras.gff3', 'Homo_sapiens.GRCh38.kras.gff3.json', 'gff3'],
        ['Homo_sapiens.GRCh38.kras.gtf', 'Homo_sapiens.GRCh38.kras.gtf.json', 'gtf'],
        ['example_genes.v2.json', 'example_genes.v3.json', 'v2-json'],
    ],
)
def test_gff_examples(filename, expected_file, input_type):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    input_file = os.path.join(data_dir, filename)
    with open(os.path.join(data_dir, expected_file), 'r') as fh:
        expected = json.load(fh)

    # order doesn't matter
    data = sort_elements(CONVERTERS[input_type](input_file))
    expected = sort_elements(expected)

    assert len(data['genes']) == len(expected['genes'])
    assert data == expected
