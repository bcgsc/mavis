from mavis.tools import _convert_tool_row, SUPPORTED_TOOL
import unittest
from mavis.constants import ORIENT, STRAND, SVTYPE
from .mock import Mock


class TestConvertToolRow(unittest.TestCase):

    def test_delly_insertion(self):
        row = Mock(
            CHROM='1', POS=247760043,
            INFO={'SVTYPE': 'INS', 'CT': 'NtoN', 'END': 247760044, 'CHR2': '1', 'CIEND': [-10, 10], 'CIPOS': [-10, 10]}
        )
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        self.assertEqual(1, len(bpp_list))
        self.assertEqual('1', bpp_list[0].break1.chr)
        self.assertEqual(247760043 - 10, bpp_list[0].break1.start)
        self.assertEqual(247760043 + 10, bpp_list[0].break1.end)
        self.assertEqual(ORIENT.LEFT, bpp_list[0].break1.orient)
        self.assertEqual(STRAND.NS, bpp_list[0].break1.strand)
        self.assertEqual(247760044 - 10, bpp_list[0].break2.start)
        self.assertEqual(247760044 + 10, bpp_list[0].break2.end)
        self.assertEqual(ORIENT.RIGHT, bpp_list[0].break2.orient)
        self.assertEqual(STRAND.NS, bpp_list[0].break2.strand)
        self.assertEqual('1', bpp_list[0].break2.chr)
        self.assertEqual(SVTYPE.INS, bpp_list[0].event_type)

    def test_delly_translocation(self):
        row = Mock(
            CHROM='7', POS=21673582,
            INFO={
                'SVTYPE': 'TRA',
                'CT': '5to5',
                'CIEND': [-700, 700],
                'CIPOS': [-700, 700],
                'CHR2': '2',
                'END': 58921502
            }
        )
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        self.assertEqual(1, len(bpp_list))
        row.INFO['CT'] = 'NtoN'
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        self.assertEqual(4, len(bpp_list))

    def test_ta_indel_insertion(self):
        raise unittest.SkipTest('TODO')

    def test_ta_indel_deletion(self):
        raise unittest.SkipTest('TODO')

    def test_ta_indel_duplication(self):
        raise unittest.SkipTest('TODO')

    def test_ta_translocation(self):
        raise unittest.SkipTest('TODO')

    def test_ta_stranded(self):
        raise unittest.SkipTest('TODO')

    def test_manta_deletion(self):
        raise unittest.SkipTest('TODO')
