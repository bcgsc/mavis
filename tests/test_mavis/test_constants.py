from mavis.constants import COLUMNS, ORIENT, STRAND, reverse_complement, sort_columns, translate


class TestConstants:
    def test_strand_compare(self):
        assert STRAND.compare(STRAND.NS, STRAND.POS)
        assert STRAND.compare(STRAND.NS, STRAND.NEG)
        assert STRAND.compare(STRAND.POS, STRAND.POS)
        assert STRAND.compare(STRAND.NEG, STRAND.NEG)
        assert not STRAND.compare(STRAND.POS, STRAND.NEG)
        assert not STRAND.compare(STRAND.NEG, STRAND.POS)

    def test_orient_compare(self):
        assert ORIENT.compare(ORIENT.NS, ORIENT.RIGHT)
        assert ORIENT.compare(ORIENT.NS, ORIENT.LEFT)
        assert ORIENT.compare(ORIENT.RIGHT, ORIENT.RIGHT)
        assert ORIENT.compare(ORIENT.LEFT, ORIENT.LEFT)
        assert not ORIENT.compare(ORIENT.RIGHT, ORIENT.LEFT)
        assert not ORIENT.compare(ORIENT.LEFT, ORIENT.RIGHT)

    def test_reverse_complement(self):
        assert reverse_complement('CGAT') == 'ATCG'
        assert reverse_complement('') == ''

    def test_translate(self):
        seq = 'ATG' 'AAT' 'TCT' 'GGA' 'TGA'
        translated_seq = translate(seq, 0)
        assert translated_seq == 'MNSG*'  # ATG AAT TCT GGA TGA
        translated_seq = translate(seq, 1)
        assert translated_seq == '*ILD'  # A TGA ATT CTG GAT GA
        translated_seq = translate(seq, 2)
        assert translated_seq == 'EFWM'  # AT GAA TTC TGG ATG A

    def test_sort_columns(self):
        temp = ['NEW', 'NEW2', COLUMNS.break1_seq, COLUMNS.break2_seq, COLUMNS.break1_chromosome]
        assert sort_columns(temp) == [
            COLUMNS.break1_chromosome,
            COLUMNS.break1_seq,
            COLUMNS.break2_seq,
            'NEW',
            'NEW2',
        ]

    def test_column_matches_column_name(self):
        assert COLUMNS.library == COLUMNS.library
        s = set([COLUMNS.library, COLUMNS.library])
        assert len(s) == 1
