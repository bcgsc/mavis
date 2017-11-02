import os
import unittest

from mavis.annotate.file_io import convert_tab_to_json, load_annotations

from . import DATA_DIR, REFERENCE_ANNOTATIONS_FILE


class TestAnnotationLoading(unittest.TestCase):
    def setUp(self):
        self.tab = os.path.join(DATA_DIR, 'annotations_subsample.tab')
        self.json = os.path.join(DATA_DIR, 'annotations_subsample.json')

    def test_convert_tab_to_json(self):
        json = convert_tab_to_json(self.tab, print)
        self.assertEqual(32, len(json['genes']))

    def test_tab_equivalent_to_json(self):
        tab_result = load_annotations(self.tab, print)
        json_result = load_annotations(self.json, print)
        self.assertEqual(sorted(tab_result.keys()), sorted(json_result.keys()))

    def test_load_tab(self):
        result = load_annotations(self.tab, print)
        self.assertEqual(12, len(result.keys()))
        domains = []
        for gene in result['12']:
            for t in gene.spliced_transcripts:
                print(t)
                if t.unspliced_transcript.name == 'ENST00000550458':
                    tl = t.translations[0]
                    domains = tl.domains
                    break
            if domains:
                break
        for d in domains:
            print(d.name, d.regions)
        self.assertEqual(2, len(domains))
        result = load_annotations(REFERENCE_ANNOTATIONS_FILE, print)
        self.assertEqual(1, len(result.keys()))

    def test_load_json(self):
        result = load_annotations(self.json, print)
        self.assertEqual(12, len(result.keys()))
