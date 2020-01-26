from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.convert import vcf

from mavis.constants import ORIENT


class TestGenerateBndAlt:
    def test_insertion(self):
        pass

    def test_duplication(self):
        pass

    def test_invertion(self):
        pass

    def test_deletion(self):
        pass

    def test_translocation(self):
        pass

    def test_inverted_translocation(self):
        pass


class TestParseBndAlt:
    def test_right(self):
        # '[4:190898243[AGGT'
        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('[4:190898243[A')
        assert '4' == chrom
        assert 190898243 == pos
        assert ORIENT.RIGHT == orient
        assert '' == seq
        assert 'A' == ref

    def test_right_untemp_seq(self):
        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('[5:190898243[AGGT')
        assert '5' == chrom
        assert 190898243 == pos
        assert ORIENT.RIGHT == orient
        assert 'AGG' == seq
        assert 'T' == ref

        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('CAGTNNNCA[5:190898243[')
        assert '5' == chrom
        assert 190898243 == pos
        assert ORIENT.RIGHT == orient
        assert 'AGTNNNCA' == seq
        assert 'C' == ref

        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('CTG[21:47575965[')
        assert '21' == chrom
        assert 47575965 == pos
        assert ORIENT.RIGHT == orient
        assert 'TG' == seq
        assert 'C' == ref

    def test_left(self):
        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('G]10:198982]')
        assert '10' == chrom
        assert 198982 == pos
        assert ORIENT.LEFT == orient
        assert '' == seq
        assert 'G' == ref

        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt(']10:198982]G')
        assert '10' == chrom
        assert 198982 == pos
        assert ORIENT.LEFT == orient
        assert '' == seq
        assert 'G' == ref

    def test_alternate_chrom(self):
        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('G]GL000.01:198982]')
        assert 'GL000.01' == chrom
        assert 198982 == pos
        assert ORIENT.LEFT == orient
        assert '' == seq
        assert 'G' == ref

    def test_left_untemp_seq(self):
        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt(']11:123456]AGTNNNCAT')
        assert '11' == chrom
        assert 123456 == pos
        assert ORIENT.LEFT == orient
        assert 'AGTNNNCA' == seq
        assert 'T' == ref

        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt(']8:1682443]TGC')
        assert '8' == chrom
        assert 1682443 == pos
        assert ORIENT.LEFT == orient
        assert 'TG' == seq
        assert 'C' == ref

        chrom, pos, orient, ref, seq = vcf.parse_bnd_alt('AAGTG]11:66289601]')
        assert '11' == chrom
        assert 66289601 == pos
        assert ORIENT.LEFT == orient
        assert 'AGTG' == seq
        assert 'A' == ref
