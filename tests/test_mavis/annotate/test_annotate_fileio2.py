from mavis.annotate.file_io import load_annotations

from ...util import get_data

JSON = get_data('annotations_subsample.json')


class TestAnnotationLoading:
    def test_load_json(self):
        result = load_annotations(JSON)
        assert len(result.keys()) == 12
