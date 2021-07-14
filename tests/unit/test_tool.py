import unittest

import pytest
from mavis.constants import COLUMNS, ORIENT, STRAND, SVTYPE
from mavis.tools import SUPPORTED_TOOL, _convert_tool_row, _parse_transabyss
from mavis.tools.vcf import convert_record as _parse_vcf_record
from mavis.tools.vcf import parse_bnd_alt as _parse_bnd_alt

from .mock import Mock


class TestDelly:
    def test_convert_insertion(self):
        row = Mock(
            chrom='1',
            pos=247760043,
            id='1DEL00000330',
            info={
                'SVTYPE': 'INS',
                'CT': 'NtoN',
                'CHR2': '1',
                'CIEND': [-10, 10],
                'CIPOS': [-10, 10],
            },
            stop=247760044,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break1.start == 247760043 - 10
        assert bpp.break1.end == 247760043 + 10
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.start == 247760044 - 10
        assert bpp.break2.end == 247760044 + 10
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.break2.strand == STRAND.NS
        assert bpp.break2.chr == '1'
        assert bpp.event_type == SVTYPE.INS
        assert bpp.untemplated_seq is None

        bpp_list = _convert_tool_row(
            _parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False, assume_no_untemplated=True
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.untemplated_seq is None
        assert bpp.untemplated_seq != ''

    def test_convert_convert_translocation(self):
        row = Mock(
            chrom='7',
            pos=21673582,
            id='TRA00016056',
            info={
                'SVTYPE': 'TRA',
                'CT': '5to5',
                'CIEND': [-700, 700],
                'CIPOS': [-700, 700],
                'CHR2': '2',
            },
            stop=58921502,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        assert len(bpp_list) == 1
        row.info['CT'] = 'NtoN'
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        assert len(bpp_list) == 4


class TestCnvNator:
    def test_convert_deletion(self):
        row = {'event_type': 'deletion', 'coordinates': '1:1-10000'}
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CNVNATOR, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 1
        assert bpp.break1.end == 1
        assert bpp.break2.start == 10000
        assert bpp.break2.start == 10000
        assert bpp.event_type == SVTYPE.DEL
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'

    def test_convert_duplication(self):
        row = {'event_type': 'duplication', 'coordinates': '1:1-10000'}
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CNVNATOR, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 1
        assert bpp.break1.end == 1
        assert bpp.break2.start == 10000
        assert bpp.break2.start == 10000
        assert bpp.event_type == SVTYPE.DUP
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'


class TestStarFusion:
    def test_convert_standard_event(self):
        row = {
            'FusionName': 'GAS6--RASA3',
            'LeftBreakpoint': 'chr13:114529969:-',
            'RightBreakpoint': 'chr13:114751269:-',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, True)

        assert len(bpp_list) == 2
        bpp = bpp_list[0]
        assert bpp.break1.chr == 'chr13'
        assert bpp.break2.chr == 'chr13'
        assert bpp.break1.start == 114529969
        assert bpp.break2.start == 114751269
        assert bpp.opposing_strands is False
        assert bpp.stranded is True

    def test_convert_translocation(self):
        row = {
            'FusionName': 'BCAS4--BCAS3',
            'LeftBreakpoint': 'chr20:49411710:+',
            'RightBreakpoint': 'chr17:59445688:+',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, True)

        assert len(bpp_list) == 2
        bpp = bpp_list[0]
        assert bpp.break1.chr == 'chr17'
        assert bpp.break2.chr == 'chr20'
        assert bpp.break1.start == 59445688
        assert bpp.break2.start == 49411710
        assert bpp.opposing_strands is False
        assert bpp.stranded is True

    def test_malformed(self):
        row = {'FusionName': 'BCAS4--BCAS3', 'LeftBreakpoint': '', 'RightBreakpoint': None}
        with pytest.raises(AssertionError):
            _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, False)


class TestTransAbyss:
    def test_convert_stranded_indel_insertion(self):
        row = {
            'chr': '1',
            'chr_start': '10015',
            'chr_end': '10015',
            'ctg_strand': '-',
            'type': 'ins',
            'alt': 'aat',
            'id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        assert len(bpp_list) == 2
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 10015
        assert bpp.break2.start == 10016
        assert bpp.event_type == SVTYPE.INS
        assert bpp.opposing_strands is False
        assert bpp.stranded is True
        assert bpp.untemplated_seq == 'AAT'

    def test_convert_indel_deletion(self):
        row = {
            'id': '1177',
            'type': 'del',
            'chr': 'X',
            'chr_start': '153523769',
            'chr_end': '153523790',
            'alt': 'na',
            'ctg_strand': '+',
            '_index': 9,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        print('after call')
        print(_convert_tool_row)
        for bpp in bpp_list:
            print(bpp)
        assert len(bpp_list) == 2
        bpp = bpp_list[0]
        assert bpp.untemplated_seq == ''

    def test_convert_indel_unstranded_insertion(self):
        row = {
            'id': '1',
            'type': 'ins',
            'chr': '1',
            'chr_start': '8877520',
            'chr_end': '8877520',
            'alt': 'tt',
            'ctg_strand': '+',
            '_index': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, False)
        print([str(b) for b in bpp_list])
        assert len(bpp_list) == 1

        bpp = bpp_list[0]
        assert bpp.event_type == SVTYPE.INS
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.strand == STRAND.NS
        assert bpp.stranded is False
        assert bpp.opposing_strands is False
        assert bpp.untemplated_seq == 'TT'

    def test_convert_indel_duplication(self):
        row = {
            'id': '1185',
            'type': 'dup',
            'chr': 'GL000211.1',
            'chr_start': '108677',
            'chr_end': '108683',
            'alt': 'aaaaaaa',
            'ctg_strand': '+',
            '_index': 15,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, False)
        print([str(b) for b in bpp_list])
        assert len(bpp_list) == 1

        bpp = bpp_list[0]
        assert bpp.event_type == SVTYPE.DUP
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.strand == STRAND.NS
        assert bpp.stranded is False
        assert bpp.opposing_strands is False
        assert bpp.untemplated_seq == ''

    def test_convert_translocation(self):
        raise unittest.SkipTest('TODO')

    def test_convert_stranded_translocation(self):
        row = {
            'strands': '+,-',
            'rearrangement': 'translocation',
            'breakpoint': '17:16342728|17:39766281',
            'orientations': 'L,L',
            'type': 'sense_fusion',
            '_index': 5261,
            'id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        assert len(bpp_list) == 2

    def test_parse_stranded_translocation(self):
        row = {
            'strands': '+,-',
            'rearrangement': 'translocation',
            'breakpoint': '17:16342728|17:39766281',
            'orientations': 'L,L',
            'type': 'sense_fusion',
            '_index': 5261,
            'id': 1,
        }
        std = _parse_transabyss(row)
        print(std)
        assert 'event_type' not in std


class TestManta:
    def test_convert_deletion(self):
        row = Mock(
            chrom='21',
            pos=9412306,
            id='MantaDEL:20644:0:2:0:0:0',
            info={'SVTYPE': 'DEL', 'CIPOS': [0, 4], 'CIEND': [0, 4]},
            stop=9412400,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.MANTA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '21'
        assert bpp.break1.start == 9412306
        assert bpp.break1.end == 9412310
        assert bpp.break2.start == 9412400
        assert bpp.break2.end == 9412404
        assert bpp.break2.chr == '21'
        print(bpp, bpp.data['tracking_id'])
        assert bpp.data['tracking_id'] == 'manta-MantaDEL:20644:0:2:0:0:0'

    def test_convert_duplication(self):
        row = Mock(
            chrom='1',
            pos=224646602,
            id='MantaDUP:TANDEM:22477:0:1:0:9:0',
            info={'SVTYPE': 'DUP', 'SVINSSEQ': 'CAAAACTTACTATAGCAGTTCTGTGAGCTGCTCTAGC'},
            stop=224800120,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.MANTA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.data['tracking_id'] == 'manta-MantaDUP:TANDEM:22477:0:1:0:9:0'

    def test_non_trans_bnd(self):
        row = Mock(
            chrom='chr1',
            pos=17051724,
            id='MantaBND:207:0:1:0:0:0:0',
            info=dict(
                SVTYPE='BND',
                MATEID='MantaBND:207:0:1:0:0:0:1',
                SVINSLEN=7,
                SVINSSEQ='GCCCCAT',
                BND_DEPTH=5,
                MATE_BND_DEPTH=4,
            ),
            ref='C',
            alts=['[chr1:234912188[GCCCCATC'],
        )
        vcf_list = _parse_vcf_record(row)
        bpp_list = _convert_tool_row(vcf_list[0], SUPPORTED_TOOL.MANTA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 17051724
        assert bpp.break2.start == 234912188
        assert bpp.break1.orient == 'R'
        assert bpp.break2.orient == 'R'
        assert bpp.data['tracking_id'] == 'manta-MantaBND:207:0:1:0:0:0:0'
        assert len(bpp_list) == 1

    def test_non_trans_bnd_from_mate(self):
        row = Mock(
            chrom='chr1',
            pos=234912188,
            id='MantaBND:207:0:1:0:0:0:1',
            info=dict(
                SVTYPE='BND',
                MATEID='MantaBND:207:0:1:0:0:0:0',
                SVINSLEN=7,
                SVINSSEQ='ATGGGGC',
                BND_DEPTH=5,
                MATE_BND_DEPTH=4,
            ),
            ref='A',
            alts=['[chr1:17051724[ATGGGGCA'],
        )
        vcf_list = _parse_vcf_record(row)
        bpp_list = _convert_tool_row(vcf_list[0], SUPPORTED_TOOL.MANTA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 17051724
        assert bpp.break2.start == 234912188
        assert bpp.break1.orient == 'R'
        assert bpp.break2.orient == 'R'
        assert bpp.data['tracking_id'] == 'manta-MantaBND:207:0:1:0:0:0:1'
        assert len(bpp_list) == 1


class TestDefuse:
    def test_convert_inverted_translocation(self):
        row = {
            'gene_chromosome1': 'X',
            'gene_chromosome2': '3',
            'genomic_break_pos1': '153063989',
            'genomic_break_pos2': '50294136',
            'genomic_strand1': '+',
            'genomic_strand2': '-',
            'cluster_id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == 'X'
        assert bpp.break1.start == 50294136
        assert bpp.break2.start == 153063989
        assert bpp.event_type is None
        assert bpp.opposing_strands is False
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.stranded is False
        assert bpp.data['tracking_id'] == 'defuse-1'

    def test_convert_translocation(self):
        row = {
            'gene_chromosome1': 'X',
            'gene_chromosome2': '3',
            'genomic_break_pos1': '153063989',
            'genomic_break_pos2': '50294136',
            'genomic_strand1': '+',
            'genomic_strand2': '+',
            'cluster_id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == 'X'
        assert bpp.break1.start == 50294136
        assert bpp.break2.start == 153063989
        assert bpp.event_type is None
        assert bpp.opposing_strands is True
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.stranded is False
        assert bpp.data['tracking_id'] == 'defuse-1'

    def test_convert_indel(self):
        row = {
            'gene_chromosome1': '1',
            'gene_chromosome2': '1',
            'genomic_break_pos1': '151732089',
            'genomic_break_pos2': '1663681',
            'genomic_strand1': '-',
            'genomic_strand2': '+',
            'cluster_id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 1663681
        assert bpp.break2.start == 151732089
        assert bpp.event_type is None
        assert bpp.opposing_strands is False
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.stranded is False
        assert bpp.data['tracking_id'] == 'defuse-1'

    def test_convert_inversion(self):
        row = {
            'gene_chromosome1': '1',
            'gene_chromosome2': '1',
            'genomic_break_pos1': '235294748',
            'genomic_break_pos2': '144898348',
            'genomic_strand1': '+',
            'genomic_strand2': '+',
            'cluster_id': 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 144898348
        assert bpp.break2.start == 235294748
        assert bpp.event_type is None
        assert bpp.opposing_strands is True
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.stranded is False
        assert bpp.data['tracking_id'] == 'defuse-1'


class TestChimerascan:
    def test_convert_pos_pos(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '+',
            'strand3p': '+',
            'chimera_cluster_id': 'CLUSTER30',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == '3'
        print(bpp)
        assert bpp.break1.start == int(row['end5p'])
        assert bpp.break2.start == int(row['start3p'])
        assert bpp.opposing_strands is False
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.stranded is False

    def test_convert_pos_neg(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '+',
            'strand3p': '-',
            'chimera_cluster_id': 'CLUSTER30',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == '3'
        print(bpp)
        assert bpp.break1.start == int(row['end5p'])
        assert bpp.break2.start == int(row['end3p'])
        assert bpp.opposing_strands is True
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.stranded is False

    def test_convert_neg_pos(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '-',
            'strand3p': '+',
            'chimera_cluster_id': 'CLUSTER30',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == '3'
        print(bpp)
        assert bpp.break1.start == int(row['start5p'])
        assert bpp.break2.start == int(row['start3p'])
        assert bpp.opposing_strands is True
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.stranded is False

    def test_convert_neg_neg(self):
        row = {
            'chrom5p': 'chr3',
            'start5p': '48599150',
            'end5p': '48601200',
            'chrom3p': 'chr3',
            'start3p': '49555116',
            'end3p': '49587666',
            'strand5p': '-',
            'strand3p': '-',
            'chimera_cluster_id': 'CLUSTER30',
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '3'
        assert bpp.break2.chr == '3'
        print(bpp)
        assert bpp.break1.start == int(row['start5p'])
        assert bpp.break2.start == int(row['end3p'])
        assert bpp.opposing_strands is False
        assert bpp.break1.orient == ORIENT.RIGHT
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.stranded is False


class TestPindel:
    def test_convert_deletion(self):
        row = Mock(chrom='21', pos=9412306, info={'SVTYPE': 'DEL'}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '21'
        assert bpp.break2.chr == '21'
        assert bpp.event_type == SVTYPE.DEL
        assert bpp.break1.start == row.pos
        assert bpp.break1.end == row.pos
        assert bpp.break2.start == row.stop
        assert bpp.break2.end == row.stop
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.break2.strand == STRAND.NS
        assert bpp.stranded is False
        assert bpp.opposing_strands is False

    def test_convert_insertion(self):
        row = Mock(chrom='21', pos=9412306, info={'SVTYPE': 'INS'}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.chr == '21'
        assert bpp.break2.chr == '21'
        assert bpp.event_type == SVTYPE.INS
        assert bpp.break1.start == row.pos
        assert bpp.break1.end == row.pos
        assert bpp.break2.start == row.stop
        assert bpp.break2.end == row.stop
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.orient == ORIENT.RIGHT
        assert bpp.break2.strand == STRAND.NS
        assert bpp.stranded is False
        assert bpp.opposing_strands is False

    def test_convert_inversion(self):
        row = Mock(chrom='21', pos=9412306, info={'SVTYPE': 'INV'}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        assert len(bpp_list) == 2
        bpp = sorted(bpp_list, key=lambda x: x.break1)[0]
        assert bpp.break1.chr == '21'
        assert bpp.break2.chr == '21'
        assert bpp.event_type == SVTYPE.INV
        assert bpp.break1.start == row.pos
        assert bpp.break1.end == row.pos
        assert bpp.break2.start == row.stop
        assert bpp.break2.end == row.stop
        assert bpp.break1.orient == ORIENT.LEFT
        assert bpp.break1.strand == STRAND.NS
        assert bpp.break2.orient == ORIENT.LEFT
        assert bpp.break2.strand == STRAND.NS
        assert bpp.stranded is False
        assert bpp.opposing_strands is True


class TestParseBndAlt:
    def test_right(self):
        # '[4:190898243[AGGT'
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('[4:190898243[A')
        assert chrom == '4'
        assert pos == 190898243
        assert orient1 == ORIENT.RIGHT
        assert orient2 == ORIENT.RIGHT
        assert seq == ''
        assert ref == 'A'

    def test_right_untemp_seq(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('[5:190898243[AGGT')
        assert chrom == '5'
        assert pos == 190898243
        assert orient1 == ORIENT.RIGHT
        assert orient2 == ORIENT.RIGHT
        assert seq == 'AGG'
        assert ref == 'T'

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('CAGTNNNCA[5:190898243[')
        assert chrom == '5'
        assert pos == 190898243
        assert orient1 == ORIENT.LEFT
        assert orient2 == ORIENT.RIGHT
        assert seq == 'AGTNNNCA'
        assert ref == 'C'

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('CTG[21:47575965[')
        assert chrom == '21'
        assert pos == 47575965
        assert orient1 == ORIENT.LEFT
        assert orient2 == ORIENT.RIGHT
        assert seq == 'TG'
        assert ref == 'C'

    def test_left(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('G]10:198982]')
        assert chrom == '10'
        assert pos == 198982
        assert orient1 == ORIENT.LEFT
        assert orient2 == ORIENT.LEFT
        assert seq == ''
        assert ref == 'G'

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt(']10:198982]G')
        assert chrom == '10'
        assert pos == 198982
        assert orient2 == ORIENT.LEFT
        assert seq == ''
        assert ref == 'G'

    def test_alternate_chrom(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('G]GL000.01:198982]')
        assert chrom == 'GL000.01'
        assert pos == 198982
        assert orient2 == ORIENT.LEFT
        assert seq == ''
        assert ref == 'G'

    def test_left_untemp_seq(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt(']11:123456]AGTNNNCAT')
        assert chrom == '11'
        assert pos == 123456
        assert orient2 == ORIENT.LEFT
        assert seq == 'AGTNNNCA'
        assert ref == 'T'

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt(']8:1682443]TGC')
        assert chrom == '8'
        assert pos == 1682443
        assert orient2 == ORIENT.LEFT
        assert seq == 'TG'
        assert ref == 'C'

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt('AAGTG]11:66289601]')
        assert chrom == '11'
        assert pos == 66289601
        assert orient2 == ORIENT.LEFT
        assert seq == 'AGTG'
        assert ref == 'A'


class TestBreakDancer:
    def test_itx(self):
        row = {
            'Chr1': '1',
            'Pos1': '10001',
            'Orientation1': '83+126-',
            'Chr2': '1',
            'Pos2': '10546',
            'Orientation2': '83+126-',
            'Type': 'ITX',
            'Size': '-352',
            'Score': '99',
            'num_Reads': '43',
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        assert len(bpps) == 1
        assert bpps[0].event_type == SVTYPE.DUP
        assert bpps[0].break1.start == 10001
        assert bpps[0].break1.end == 10001
        assert bpps[0].break1.orient == ORIENT.RIGHT
        assert bpps[0].break2.start == 10546
        assert bpps[0].break2.end == 10546
        assert bpps[0].break2.orient == ORIENT.LEFT
        assert bpps[0].opposing_strands is False

    def test_deletion(self):
        row = {
            'Chr1': '1',
            'Pos1': '869445',
            'Orientation1': '89+21-',
            'Chr2': '1',
            'Pos2': '870225',
            'Orientation2': '5+93-',
            'Type': 'DEL',
            'Size': '892',
            'Score': '99',
            'num_Reads': '67',
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        assert len(bpps) == 1
        assert bpps[0].event_type == SVTYPE.DEL
        assert bpps[0].break1.start == 869445
        assert bpps[0].break1.end == 869445
        assert bpps[0].break1.orient == ORIENT.LEFT
        assert bpps[0].break2.start == 870225
        assert bpps[0].break2.end == 870225
        assert bpps[0].break2.orient == ORIENT.RIGHT
        assert bpps[0].opposing_strands is False

    def test_inversion(self):
        row = {
            'Chr1': '1',
            'Pos1': '13143396',
            'Orientation1': '18+4-',
            'Chr2': '1',
            'Pos2': '13218683',
            'Orientation2': '10+10-',
            'Type': 'INV',
            'Size': '74618',
            'Score': '31',
            'num_Reads': '2',
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        assert len(bpps) == 2
        assert bpps[0].event_type == SVTYPE.INV
        assert bpps[0].break1.start == 13143396
        assert bpps[0].break1.end == 13143396
        assert bpps[0].break1.orient == ORIENT.LEFT
        assert bpps[0].break2.start == 13218683
        assert bpps[0].break2.end == 13218683
        assert bpps[0].break2.orient == ORIENT.LEFT
        assert bpps[0].opposing_strands is True

        assert bpps[1].event_type == SVTYPE.INV
        assert bpps[1].break1.start == 13143396
        assert bpps[1].break1.end == 13143396
        assert bpps[1].break1.orient == ORIENT.RIGHT
        assert bpps[1].break2.start == 13218683
        assert bpps[1].break2.end == 13218683
        assert bpps[1].break2.orient == ORIENT.RIGHT
        assert bpps[1].opposing_strands is True

    def test_insertion(self):
        row = {
            'Chr1': '1',
            'Pos1': '20216146',
            'Orientation1': '23+26-',
            'Chr2': '1',
            'Pos2': '20218060',
            'Orientation2': '23+26-',
            'Type': 'INS',
            'Size': '-421',
            'Score': '99',
            'num_Reads': '3',
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        assert len(bpps) == 1
        assert bpps[0].event_type == SVTYPE.INS
        assert bpps[0].break1.start == 20216146
        assert bpps[0].break1.end == 20216146
        assert bpps[0].break1.orient == ORIENT.LEFT
        assert bpps[0].break2.start == 20218060
        assert bpps[0].break2.end == 20218060
        assert bpps[0].break2.orient == ORIENT.RIGHT
        assert bpps[0].opposing_strands is False


class TestStrelka:
    def testInsertion(self):
        event = Mock(
            chrom='1', pos=724986, id=None, info={}, ref='G', stop=724986, alts=('GGAATT',)
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.STRELKA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 724986
        assert bpp.break1.end == 724986
        assert bpp.break2.start == 724986
        assert bpp.break2.end == 724986
        assert bpp.event_type == SVTYPE.INS

    def testDeletion(self):
        event = Mock(
            chrom='1',
            pos=1265353,
            id=None,
            info={},
            ref='GCGTGTGCCATGCA',
            stop=1265366,
            alts=('G',),
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.STRELKA, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 1265353
        assert bpp.break1.end == 1265353
        assert bpp.break2.start == 1265366
        assert bpp.break2.end == 1265366
        assert bpp.event_type == SVTYPE.DEL

    def testMalformated(self):
        event = Mock(
            chrom='1',
            pos=53678660,
            id=None,
            info={'SVTYPE': 'BND'},
            ref='C',
            alts=('CTTTTAAATGTAACATGACATAATATATTTCCTAAATAATTTAAAATAATC.',),
            stop=53678660,
        )
        with pytest.raises(NotImplementedError):
            _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.STRELKA, False)


class TestMutect(unittest.TestCase):
    def testInsertion(self):
        event = Mock(
            chrom='1', pos=724986, id=None, info={}, ref='G', stop=724986, alts=('GGAATT',)
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.MUTECT, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 724986
        assert bpp.break1.end == 724986
        assert bpp.break2.start == 724986
        assert bpp.break2.end == 724986
        assert bpp.event_type == SVTYPE.INS

    def testDeletion(self):
        event = Mock(
            chrom='1',
            pos=1265353,
            id=None,
            info={},
            ref='GCGTGTGCCATGCA',
            stop=1265366,
            alts=('G',),
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.MUTECT, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 1265353
        assert bpp.break1.end == 1265353
        assert bpp.break2.start == 1265366
        assert bpp.break2.end == 1265366
        assert bpp.event_type == SVTYPE.DEL

    def testMalformated(self):
        event = Mock(
            chrom='1',
            pos=53678660,
            id=None,
            info={'SVTYPE': 'BND'},
            ref='C',
            alts=('CTTTTAAATGTAACATGACATAATATATTTCCTAAATAATTTAAAATAATC.',),
            stop=53678660,
        )
        with self.assertRaises(NotImplementedError):
            _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.MUTECT, False)


class TestVCF(unittest.TestCase):
    def setUp(self):
        self.tra = Mock(
            chrom='2',
            pos=21673582,
            id=None,
            info={'SVTYPE': 'TRA', 'CT': '5to5', 'CHR2': '3'},
            stop=58921502,
            alts=[],
        )


@pytest.fixture
def vcf_translocation():
    return Mock(
        chrom='2',
        pos=21673582,
        id=None,
        info={'SVTYPE': 'TRA', 'CT': '5to5', 'CHR2': '3'},
        stop=58921502,
        alts=[],
    )


class TestVCF:
    def test_no_ci(self, vcf_translocation):
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 21673582
        assert bpp.break1.end == 21673582
        assert bpp.break2.start == 58921502
        assert bpp.break2.end == 58921502

    def test_ci(self, vcf_translocation):
        vcf_translocation.info.update({'CIEND': [-700, 700], 'CIPOS': [-700, 700]})
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        print(bpp)
        assert bpp.break1.start == 21673582 - 700
        assert bpp.break1.end == 21673582 + 700
        assert bpp.break2.start == 58921502 - 700
        assert bpp.break2.end == 58921502 + 700

    def test_precise_flag_ignores_ci(self, vcf_translocation):
        vcf_translocation.info.update({'CIEND': [-700, 700], 'CIPOS': [-700, 700], 'PRECISE': True})
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 21673582
        assert bpp.break1.end == 21673582
        assert bpp.break2.start == 58921502
        assert bpp.break2.end == 58921502

    def test_no_id(self, vcf_translocation):
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.data[COLUMNS.tracking_id]

    def test_N_id(self, vcf_translocation):
        vcf_translocation.id = 'N'
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.data[COLUMNS.tracking_id]
        assert bpp.data[COLUMNS.tracking_id] != 'N'

    def test_id_given(self, vcf_translocation):
        vcf_translocation.id = 'thing-1'
        bpp_list = _convert_tool_row(
            _parse_vcf_record(vcf_translocation)[0], SUPPORTED_TOOL.VCF, False
        )
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.data[COLUMNS.tracking_id] == 'vcf-thing-1'
