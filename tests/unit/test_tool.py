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

    def test_ta_stranded_indel_insertion(self):
        row = {
            'chr': '1', 'chr_start': '10015', 'chr_end': '10015', 'ctg_strand': '-', 'type': 'ins'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('1', bpp.break1.chr)
        self.assertEqual('1', bpp.break2.chr)
        self.assertEqual(10015, bpp.break1.start)
        self.assertEqual(10016, bpp.break2.start)
        self.assertEqual(SVTYPE.INS, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(STRAND.POS, bpp.break1.strand)
        self.assertEqual(STRAND.POS, bpp.break2.strand)
        self.assertEqual(True, bpp.stranded)

    def test_ta_indel_deletion(self):
        raise unittest.SkipTest('TODO')

    def test_ta_indel_duplication(self):
        raise unittest.SkipTest('TODO')

    def test_ta_translocation(self):
        raise unittest.SkipTest('TODO')

    def test_ta_stranded_translocation(self):
        row = {
            'strands': '+,-',
            'rearrangement': 'translocation',
            'breakpoint': '17:16342728|17:39766281',
            'orientations': 'L,L',
            'type': 'sense_fusion',
            '_index': 5261
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        self.assertEqual(1, len(bpp_list))

    def test_manta_deletion(self):
        raise unittest.SkipTest('TODO')

    def test_defuse_inverted_translocation(self):
        row = {
            'gene_chromosome1': 'X',
            'gene_chromosome2': '3',
            'genomic_break_pos1': '153063989',
            'genomic_break_pos2': '50294136',
            'genomic_strand1': '+',
            'genomic_strand2': '-'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('X', bpp.break2.chr)
        self.assertEqual(50294136, bpp.break1.start)
        self.assertEqual(153063989, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_defuse_translocation(self):
        row = {
            'gene_chromosome1': 'X',
            'gene_chromosome2': '3',
            'genomic_break_pos1': '153063989',
            'genomic_break_pos2': '50294136',
            'genomic_strand1': '+',
            'genomic_strand2': '+'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('X', bpp.break2.chr)
        self.assertEqual(50294136, bpp.break1.start)
        self.assertEqual(153063989, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_defuse_indel(self):
        row = {
            'gene_chromosome1': '1',
            'gene_chromosome2': '1',
            'genomic_break_pos1': '151732089',
            'genomic_break_pos2': '1663681',
            'genomic_strand1': '-',
            'genomic_strand2': '+'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('1', bpp.break1.chr)
        self.assertEqual('1', bpp.break2.chr)
        self.assertEqual(1663681, bpp.break1.start)
        self.assertEqual(151732089, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_defuse_inversion(self):
        row = {
            'gene_chromosome1': '1',
            'gene_chromosome2': '1',
            'genomic_break_pos1': '235294748',
            'genomic_break_pos2': '144898348',
            'genomic_strand1': '+',
            'genomic_strand2': '+'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('1', bpp.break1.chr)
        self.assertEqual('1', bpp.break2.chr)
        self.assertEqual(144898348, bpp.break1.start)
        self.assertEqual(235294748, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_chimerascan_pos_pos(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '+',
            'strand3p': '+'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('3', bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row['end5p']), bpp.break1.start)
        self.assertEqual(int(row['start3p']), bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_chimerascan_pos_neg(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '+',
            'strand3p': '-'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('3', bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row['end5p']), bpp.break1.start)
        self.assertEqual(int(row['end3p']), bpp.break2.start)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_chimerascan_neg_pos(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '-',
            'strand3p': '+'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('3', bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row['start5p']), bpp.break1.start)
        self.assertEqual(int(row['start3p']), bpp.break2.start)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_chimerascan_neg_neg(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '-',
            'strand3p': '-'
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual('3', bpp.break1.chr)
        self.assertEqual('3', bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row['start5p']), bpp.break1.start)
        self.assertEqual(int(row['end3p']), bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)
