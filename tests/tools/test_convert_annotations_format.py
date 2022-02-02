import os

from tools.convert_annotations_format import convert_gff2_to_mavis, convert_gff3_to_mavis


def test_load_gff3():
    input = os.path.join(os.path.dirname(__file__), 'data', 'Homo_sapiens.GRCh38.105.chr.kras.gtf')
    data = convert_gff2_to_mavis(input, False)
    assert len(data['genes']) == 2
    assert sum([len(g['transcripts']) for g in data['genes']]) == 15
    exons = 0
    for gene in data['genes']:
        for transcript in gene['transcripts']:
            exons += len(transcript['exons'])
    assert exons == 62


def test_load_gtf():
    input = os.path.join(os.path.dirname(__file__), 'data', 'Homo_sapiens.GRCh38.105.kras.gff3')
    data = convert_gff3_to_mavis(input, False)
    assert len(data['genes']) == 4
    assert sum([len(g['transcripts']) for g in data['genes']]) == 15
