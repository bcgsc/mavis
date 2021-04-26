from mavis.annotate.file_io import convert_tab_to_json, load_annotations

from ..util import get_data

TAB = get_data('annotations_subsample.tab')
JSON = get_data('annotations_subsample.json')


class TestAnnotationLoading:
    def test_convert_tab_to_json(self):
        json = convert_tab_to_json(TAB, warn=print)
        assert len(json['genes']) == 32

    def test_tab_equivalent_to_json(self):
        tab_result = load_annotations(TAB, warn=print)
        json_result = load_annotations(JSON, warn=print)
        assert sorted(json_result.keys()) == sorted(tab_result.keys())

    def test_load_tab(self):
        result = load_annotations(TAB, warn=print)
        assert len(result.keys()) == 12
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
        assert len(domains) == 2
        result = load_annotations(get_data('mock_reference_annotations.tsv'), warn=print)
        assert len(result.keys()) == 1

    def test_load_json(self):
        result = load_annotations(JSON, warn=print)
        assert len(result.keys()) == 12
