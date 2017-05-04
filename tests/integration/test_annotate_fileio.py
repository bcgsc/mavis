import unittest
import os
from mavis.annotate.file_io import load_annotations, convert_tab_to_json
from . import DATA_DIR, REFERENCE_ANNOTATIONS_FILE


class TestAnnotationLoading(unittest.TestCase):
    def setUp(self):
        self.tab = os.path.join(DATA_DIR, 'annotations_subsample.tab')
        self.json = os.path.join(DATA_DIR, 'annotations_subsample.json')

    def test_convert_tab_to_json(self):
        json = convert_tab_to_json(self.tab, print)
        self.assertEqual(31, len(json['genes']))

    def test_tab_equivalent_to_json(self):
        tab_result = load_annotations(self.tab, print)
        json_result = load_annotations(self.json, print)
        self.assertEqual(sorted(tab_result.keys()), sorted(json_result.keys()))

    def test_load_tab(self):
        result = load_annotations(self.tab, print)
        self.assertEqual(12, len(result.keys()))
        
        result = load_annotations(REFERENCE_ANNOTATIONS_FILE, print)
        self.assertEqual(1, len(result.keys()))


    def test_load_json(self):
        result = load_annotations(self.json, print)
        self.assertEqual(12, len(result.keys()))
