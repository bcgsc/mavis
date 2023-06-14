import unittest

from mavis.constants import COLUMNS, ORIENT, STRAND, SVTYPE
from mavis.tools import (
    _convert_tool_row,
    SUPPORTED_TOOL,
    _parse_transabyss,
)
from mavis.tools.vcf import parse_bnd_alt as _parse_bnd_alt, convert_record as _parse_vcf_record

from .mock import Mock


class TestDelly(unittest.TestCase):
    def test_convert_insertion(self):
        row = Mock(
            chrom="1",
            pos=247760043,
            id="1DEL00000330",
            info={
                "SVTYPE": "INS",
                "CT": "NtoN",
                "CHR2": "1",
                "CIEND": [-10, 10],
                "CIPOS": [-10, 10],
            },
            stop=247760044,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual(247760043 - 10, bpp.break1.start)
        self.assertEqual(247760043 + 10, bpp.break1.end)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(247760044 - 10, bpp.break2.start)
        self.assertEqual(247760044 + 10, bpp.break2.end)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(SVTYPE.INS, bpp.event_type)
        self.assertEqual(None, bpp.untemplated_seq)

        bpp_list = _convert_tool_row(
            _parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False, assume_no_untemplated=True
        )
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(None, bpp.untemplated_seq)
        self.assertNotEqual("", bpp.untemplated_seq)

    def test_convert_convert_translocation(self):
        row = Mock(
            chrom="7",
            pos=21673582,
            id="TRA00016056",
            info={
                "SVTYPE": "TRA",
                "CT": "5to5",
                "CIEND": [-700, 700],
                "CIPOS": [-700, 700],
                "CHR2": "2",
            },
            stop=58921502,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        self.assertEqual(1, len(bpp_list))
        row.info["CT"] = "NtoN"
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.DELLY, False)
        for b in bpp_list:
            print(b)
        self.assertEqual(4, len(bpp_list))


class TestCnvNator(unittest.TestCase):
    def test_convert_deletion(self):
        row = {"event_type": "deletion", "coordinates": "1:1-10000"}
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CNVNATOR, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(1, bpp.break1.start)
        self.assertEqual(1, bpp.break1.end)
        self.assertEqual(10000, bpp.break2.start)
        self.assertEqual(10000, bpp.break2.start)
        self.assertEqual(SVTYPE.DEL, bpp.event_type)
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)

    def test_convert_duplication(self):
        row = {"event_type": "duplication", "coordinates": "1:1-10000"}
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CNVNATOR, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(1, bpp.break1.start)
        self.assertEqual(1, bpp.break1.end)
        self.assertEqual(10000, bpp.break2.start)
        self.assertEqual(10000, bpp.break2.start)
        self.assertEqual(SVTYPE.DUP, bpp.event_type)
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)


class TestArriba(unittest.TestCase):
    def test_convert_standard_event(self):
        row = {
            "breakpoint1": "13:114529969",
            "breakpoint2": "13:114751269",
            "type": "inversion",
            "strand1(gene/fusion)": "+/+",
            "strand2(gene/fusion)": "-/-",
            "direction1": "downstream",
            "direction2": "downstream",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.ARRIBA, True)

        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("13", bpp.break1.chr)
        self.assertEqual("13", bpp.break2.chr)
        self.assertEqual(114529969, bpp.break1.start)
        self.assertEqual(114751269, bpp.break2.start)
        self.assertEqual(SVTYPE.INV, bpp.event_type)
        self.assertEqual("L", bpp.break1.orient)
        self.assertEqual("L", bpp.break2.orient)
        self.assertEqual(True, bpp.opposing_strands)

    def test_convert_translocation(self):
        row = {
            "breakpoint1": "17:69313092",
            "breakpoint2": "20:58272875",
            "type": "translocation",
            "strand1(gene/fusion)": "-/-",
            "strand2(gene/fusion)": "-/-",
            "direction1": "upstream",
            "direction2": "downstream",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.ARRIBA, True)

        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr17", bpp.break1.chr)
        self.assertEqual("chr20", bpp.break2.chr)
        self.assertEqual(69313092, bpp.break1.start)
        self.assertEqual(58272875, bpp.break2.start)
        self.assertEqual(SVTYPE.TRANS, bpp.event_type)
        self.assertEqual("R", bpp.break1.orient)
        self.assertEqual("L", bpp.break2.orient)
        self.assertEqual(False, bpp.opposing_strands)

    def test_convert_translocation(self):
        row = {
            "breakpoint1": "20:57265705",
            "breakpoint2": "20:47786405",
            "type": "inversion/5'-5'",
            "strand1(gene/fusion)": "-/-",
            "strand2(gene/fusion)": "-/+",
            "direction1": "upstream",
            "direction2": "upstream",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.ARRIBA, True)

        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr20", bpp.break1.chr)
        self.assertEqual("chr20", bpp.break2.chr)
        self.assertEqual(57265705, bpp.break1.start)
        self.assertEqual(47786405, bpp.break2.start)
        self.assertEqual(SVTYPE.INV, bpp.event_type)
        self.assertEqual("R", bpp.break1.orient)
        self.assertEqual("R", bpp.break2.orient)
        self.assertEqual(True, bpp.opposing_strands)

    def test_convert_translocation(self):
        row = {
            "breakpoint1": "14:102877322",
            "breakpoint2": "14:102994672",
            "type": "deletion/read-through/5'-5'",
            "strand1(gene/fusion)": "+/+",
            "strand2(gene/fusion)": "-/+",
            "direction1": "downstream",
            "direction2": "upstream",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.ARRIBA, True)

        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr14", bpp.break1.chr)
        self.assertEqual("chr14", bpp.break2.chr)
        self.assertEqual(102877322, bpp.break1.start)
        self.assertEqual(102994672, bpp.break2.start)
        self.assertEqual(SVTYPE.DEL, bpp.event_type)
        self.assertEqual("L", bpp.break1.orient)
        self.assertEqual("R", bpp.break2.orient)
        self.assertEqual(False, bpp.opposing_strands)

    def test_malformed(self):
        row = {
            "breakpoint1": "",
            "breakpoint2": None,
            "type": "translocation",
            "strand1(gene/fusion)": "-/-",
            "strand2(gene/fusion)": "-/-",
            "direction1": "upstream",
            "direction2": "downstream",
        }
        with self.assertRaises(AssertionError):
            _convert_tool_row(row, SUPPORTED_TOOL.ARRIBA, False)


class TestStarFusion(unittest.TestCase):
    def test_convert_standard_event(self):
        row = {
            "FusionName": "GAS6--RASA3",
            "LeftBreakpoint": "chr13:114529969:-",
            "RightBreakpoint": "chr13:114751269:-",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, True)

        self.assertEqual(2, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr13", bpp.break1.chr)
        self.assertEqual("chr13", bpp.break2.chr)
        self.assertEqual(114529969, bpp.break1.start)
        self.assertEqual(114751269, bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(True, bpp.stranded)

    def test_convert_translocation(self):
        row = {
            "FusionName": "BCAS4--BCAS3",
            "LeftBreakpoint": "chr20:49411710:+",
            "RightBreakpoint": "chr17:59445688:+",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, True)

        self.assertEqual(2, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr17", bpp.break1.chr)
        self.assertEqual("chr20", bpp.break2.chr)
        self.assertEqual(59445688, bpp.break1.start)
        self.assertEqual(49411710, bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(True, bpp.stranded)

    def test_malformed(self):
        row = {"FusionName": "BCAS4--BCAS3", "LeftBreakpoint": "", "RightBreakpoint": None}
        with self.assertRaises(AssertionError):
            _convert_tool_row(row, SUPPORTED_TOOL.STARFUSION, False)


class TestTransAbyss(unittest.TestCase):
    def test_convert_stranded_indel_insertion(self):
        row = {
            "chr": "1",
            "chr_start": "10015",
            "chr_end": "10015",
            "ctg_strand": "-",
            "type": "ins",
            "alt": "aat",
            "id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        self.assertEqual(2, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(10015, bpp.break1.start)
        self.assertEqual(10016, bpp.break2.start)
        self.assertEqual(SVTYPE.INS, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(True, bpp.stranded)
        self.assertEqual("AAT", bpp.untemplated_seq)

    def test_convert_indel_deletion(self):
        row = {
            "id": "1177",
            "type": "del",
            "chr": "X",
            "chr_start": "153523769",
            "chr_end": "153523790",
            "alt": "na",
            "ctg_strand": "+",
            "_index": 9,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        print("after call")
        print(_convert_tool_row)
        for bpp in bpp_list:
            print(bpp)
        self.assertEqual(2, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("", bpp.untemplated_seq)

    def test_convert_indel_unstranded_insertion(self):
        row = {
            "id": "1",
            "type": "ins",
            "chr": "1",
            "chr_start": "8877520",
            "chr_end": "8877520",
            "alt": "tt",
            "ctg_strand": "+",
            "_index": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, False)
        print([str(b) for b in bpp_list])
        self.assertEqual(1, len(bpp_list))

        bpp = bpp_list[0]
        self.assertEqual(SVTYPE.INS, bpp.event_type)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual("TT", bpp.untemplated_seq)

    def test_convert_indel_duplication(self):
        row = {
            "id": "1185",
            "type": "dup",
            "chr": "GL000211.1",
            "chr_start": "108677",
            "chr_end": "108683",
            "alt": "aaaaaaa",
            "ctg_strand": "+",
            "_index": 15,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, False)
        print([str(b) for b in bpp_list])
        self.assertEqual(1, len(bpp_list))

        bpp = bpp_list[0]
        self.assertEqual(SVTYPE.DUP, bpp.event_type)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual("", bpp.untemplated_seq)

    def test_convert_translocation(self):
        raise unittest.SkipTest("TODO")

    def test_convert_stranded_translocation(self):
        row = {
            "strands": "+,-",
            "rearrangement": "translocation",
            "breakpoint": "17:16342728|17:39766281",
            "orientations": "L,L",
            "type": "sense_fusion",
            "_index": 5261,
            "id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.TA, True)
        self.assertEqual(2, len(bpp_list))

    def test_parse_stranded_translocation(self):
        row = {
            "strands": "+,-",
            "rearrangement": "translocation",
            "breakpoint": "17:16342728|17:39766281",
            "orientations": "L,L",
            "type": "sense_fusion",
            "_index": 5261,
            "id": 1,
        }
        std = _parse_transabyss(row)
        print(std)
        self.assertTrue("event_type" not in std)


class TestManta(unittest.TestCase):
    def test_convert_deletion(self):
        row = Mock(
            chrom="21",
            pos=9412306,
            id="MantaDEL:20644:0:2:0:0:0",
            info={"SVTYPE": "DEL", "CIPOS": [0, 4], "CIEND": [0, 4]},
            stop=9412400,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("21", bpp.break1.chr)
        self.assertEqual(9412306, bpp.break1.start)
        self.assertEqual(9412310, bpp.break1.end)
        self.assertEqual(9412400, bpp.break2.start)
        self.assertEqual(9412404, bpp.break2.end)
        self.assertEqual("21", bpp.break2.chr)
        print(bpp, bpp.tracking_id)
        self.assertEqual("manta-MantaDEL:20644:0:2:0:0:0", bpp.tracking_id)

    def test_convert_duplication(self):
        row = Mock(
            chrom="1",
            pos=224646602,
            id="MantaDUP:TANDEM:22477:0:1:0:9:0",
            info={"SVTYPE": "DUP", "SVINSSEQ": "CAAAACTTACTATAGCAGTTCTGTGAGCTGCTCTAGC"},
            stop=224800120,
            alts=[],
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual("manta-MantaDUP:TANDEM:22477:0:1:0:9:0", bpp.tracking_id)

    def test_non_trans_bnd(self):
        row = Mock(
            chrom="chr1",
            pos=17051724,
            id="MantaBND:207:0:1:0:0:0:0",
            info=dict(
                SVTYPE="BND",
                MATEID="MantaBND:207:0:1:0:0:0:1",
                SVINSLEN=7,
                SVINSSEQ="GCCCCAT",
                BND_DEPTH=5,
                MATE_BND_DEPTH=4,
            ),
            ref="C",
            alts=["[chr1:234912188[GCCCCATC"],
        )
        vcf_list = _parse_vcf_record(row)
        bpp_list = _convert_tool_row(vcf_list[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(17051724, bpp.break1.start)
        self.assertEqual(234912188, bpp.break2.start)
        self.assertEqual("R", bpp.break1.orient)
        self.assertEqual("R", bpp.break2.orient)
        self.assertEqual("manta-MantaBND:207:0:1:0:0:0:0", bpp.tracking_id)
        self.assertEqual(1, len(bpp_list))

    def test_non_trans_bnd_from_mate(self):
        row = Mock(
            chrom="chr1",
            pos=234912188,
            id="MantaBND:207:0:1:0:0:0:1",
            info=dict(
                SVTYPE="BND",
                MATEID="MantaBND:207:0:1:0:0:0:0",
                SVINSLEN=7,
                SVINSSEQ="ATGGGGC",
                BND_DEPTH=5,
                MATE_BND_DEPTH=4,
            ),
            ref="A",
            alts=["[chr1:17051724[ATGGGGCA"],
        )
        vcf_list = _parse_vcf_record(row)
        bpp_list = _convert_tool_row(vcf_list[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(17051724, bpp.break1.start)
        self.assertEqual(234912188, bpp.break2.start)
        self.assertEqual("R", bpp.break1.orient)
        self.assertEqual("R", bpp.break2.orient)
        self.assertEqual("manta-MantaBND:207:0:1:0:0:0:1", bpp.tracking_id)
        self.assertEqual(1, len(bpp_list))


class TestDefuse(unittest.TestCase):
    def test_convert_inverted_translocation(self):
        row = {
            "gene_chromosome1": "X",
            "gene_chromosome2": "3",
            "genomic_break_pos1": "153063989",
            "genomic_break_pos2": "50294136",
            "genomic_strand1": "+",
            "genomic_strand2": "-",
            "cluster_id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("X", bpp.break2.chr)
        self.assertEqual(50294136, bpp.break1.start)
        self.assertEqual(153063989, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual("defuse-1", bpp.tracking_id)

    def test_convert_translocation(self):
        row = {
            "gene_chromosome1": "X",
            "gene_chromosome2": "3",
            "genomic_break_pos1": "153063989",
            "genomic_break_pos2": "50294136",
            "genomic_strand1": "+",
            "genomic_strand2": "+",
            "cluster_id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("X", bpp.break2.chr)
        self.assertEqual(50294136, bpp.break1.start)
        self.assertEqual(153063989, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual("defuse-1", bpp.tracking_id)

    def test_convert_indel(self):
        row = {
            "gene_chromosome1": "1",
            "gene_chromosome2": "1",
            "genomic_break_pos1": "151732089",
            "genomic_break_pos2": "1663681",
            "genomic_strand1": "-",
            "genomic_strand2": "+",
            "cluster_id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(1663681, bpp.break1.start)
        self.assertEqual(151732089, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual("defuse-1", bpp.tracking_id)

    def test_convert_inversion(self):
        row = {
            "gene_chromosome1": "1",
            "gene_chromosome2": "1",
            "genomic_break_pos1": "235294748",
            "genomic_break_pos2": "144898348",
            "genomic_strand1": "+",
            "genomic_strand2": "+",
            "cluster_id": 1,
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.DEFUSE, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("1", bpp.break1.chr)
        self.assertEqual("1", bpp.break2.chr)
        self.assertEqual(144898348, bpp.break1.start)
        self.assertEqual(235294748, bpp.break2.start)
        self.assertEqual(None, bpp.event_type)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual("defuse-1", bpp.tracking_id)


class TestChimerascan(unittest.TestCase):
    def test_convert_pos_pos(self):
        row = {
            "chrom5p": "chr3",
            "start5p": "48599150",
            "end5p": "48601200",
            "chrom3p": "chr3",
            "start3p": "49555116",
            "end3p": "49587666",
            "strand5p": "+",
            "strand3p": "+",
            "chimera_cluster_id": "CLUSTER30",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("3", bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row["end5p"]), bpp.break1.start)
        self.assertEqual(int(row["start3p"]), bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_convert_pos_neg(self):
        row = {
            "chrom5p": "chr3",
            "start5p": "48599150",
            "end5p": "48601200",
            "chrom3p": "chr3",
            "start3p": "49555116",
            "end3p": "49587666",
            "strand5p": "+",
            "strand3p": "-",
            "chimera_cluster_id": "CLUSTER30",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("3", bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row["end5p"]), bpp.break1.start)
        self.assertEqual(int(row["end3p"]), bpp.break2.start)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_convert_neg_pos(self):
        row = {
            "chrom5p": "chr3",
            "start5p": "48599150",
            "end5p": "48601200",
            "chrom3p": "chr3",
            "start3p": "49555116",
            "end3p": "49587666",
            "strand5p": "-",
            "strand3p": "+",
            "chimera_cluster_id": "CLUSTER30",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("3", bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row["start5p"]), bpp.break1.start)
        self.assertEqual(int(row["start3p"]), bpp.break2.start)
        self.assertEqual(True, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)

    def test_convert_neg_neg(self):
        row = {
            "chrom5p": "chr3",
            "start5p": "48599150",
            "end5p": "48601200",
            "chrom3p": "chr3",
            "start3p": "49555116",
            "end3p": "49587666",
            "strand5p": "-",
            "strand3p": "-",
            "chimera_cluster_id": "CLUSTER30",
        }
        bpp_list = _convert_tool_row(row, SUPPORTED_TOOL.CHIMERASCAN, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("3", bpp.break1.chr)
        self.assertEqual("3", bpp.break2.chr)
        print(bpp)
        self.assertEqual(int(row["start5p"]), bpp.break1.start)
        self.assertEqual(int(row["end3p"]), bpp.break2.start)
        self.assertEqual(False, bpp.opposing_strands)
        self.assertEqual(ORIENT.RIGHT, bpp.break1.orient)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(False, bpp.stranded)


class TestPindel(unittest.TestCase):
    def test_convert_deletion(self):
        row = Mock(chrom="21", pos=9412306, info={"SVTYPE": "DEL"}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("21", bpp.break1.chr)
        self.assertEqual("21", bpp.break2.chr)
        self.assertEqual(SVTYPE.DEL, bpp.event_type)
        self.assertEqual(row.pos, bpp.break1.start)
        self.assertEqual(row.pos, bpp.break1.end)
        self.assertEqual(row.stop, bpp.break2.start)
        self.assertEqual(row.stop, bpp.break2.end)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual(False, bpp.opposing_strands)

    def test_convert_insertion(self):
        row = Mock(chrom="21", pos=9412306, info={"SVTYPE": "INS"}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("21", bpp.break1.chr)
        self.assertEqual("21", bpp.break2.chr)
        self.assertEqual(SVTYPE.INS, bpp.event_type)
        self.assertEqual(row.pos, bpp.break1.start)
        self.assertEqual(row.pos, bpp.break1.end)
        self.assertEqual(row.stop, bpp.break2.start)
        self.assertEqual(row.stop, bpp.break2.end)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(ORIENT.RIGHT, bpp.break2.orient)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual(False, bpp.opposing_strands)

    def test_convert_inversion(self):
        row = Mock(chrom="21", pos=9412306, info={"SVTYPE": "INV"}, stop=9412400, id=None, alts=[])
        bpp_list = _convert_tool_row(_parse_vcf_record(row)[0], SUPPORTED_TOOL.PINDEL, False)
        self.assertEqual(2, len(bpp_list))
        bpp = sorted(bpp_list, key=lambda x: x.break1)[0]
        self.assertEqual("21", bpp.break1.chr)
        self.assertEqual("21", bpp.break2.chr)
        self.assertEqual(SVTYPE.INV, bpp.event_type)
        self.assertEqual(row.pos, bpp.break1.start)
        self.assertEqual(row.pos, bpp.break1.end)
        self.assertEqual(row.stop, bpp.break2.start)
        self.assertEqual(row.stop, bpp.break2.end)
        self.assertEqual(ORIENT.LEFT, bpp.break1.orient)
        self.assertEqual(STRAND.NS, bpp.break1.strand)
        self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        self.assertEqual(STRAND.NS, bpp.break2.strand)
        self.assertEqual(False, bpp.stranded)
        self.assertEqual(True, bpp.opposing_strands)


class TestParseBndAlt(unittest.TestCase):
    def test_right(self):
        # '[4:190898243[AGGT'
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("[4:190898243[A")
        self.assertEqual("4", chrom)
        self.assertEqual(190898243, pos)
        self.assertEqual(ORIENT.RIGHT, orient1)
        self.assertEqual(ORIENT.RIGHT, orient2)
        self.assertEqual("", seq)
        self.assertEqual("A", ref)

    def test_right_untemp_seq(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("[5:190898243[AGGT")
        self.assertEqual("5", chrom)
        self.assertEqual(190898243, pos)
        self.assertEqual(ORIENT.RIGHT, orient1)
        self.assertEqual(ORIENT.RIGHT, orient2)
        self.assertEqual("AGG", seq)
        self.assertEqual("T", ref)

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("CAGTNNNCA[5:190898243[")
        self.assertEqual("5", chrom)
        self.assertEqual(190898243, pos)
        self.assertEqual(ORIENT.LEFT, orient1)
        self.assertEqual(ORIENT.RIGHT, orient2)
        self.assertEqual("AGTNNNCA", seq)
        self.assertEqual("C", ref)

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("CTG[21:47575965[")
        self.assertEqual("21", chrom)
        self.assertEqual(47575965, pos)
        self.assertEqual(ORIENT.LEFT, orient1)
        self.assertEqual(ORIENT.RIGHT, orient2)
        self.assertEqual("TG", seq)
        self.assertEqual("C", ref)

    def test_left(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("G]10:198982]")
        self.assertEqual("10", chrom)
        self.assertEqual(198982, pos)
        self.assertEqual(ORIENT.LEFT, orient1)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("", seq)
        self.assertEqual("G", ref)

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("]10:198982]G")
        self.assertEqual("10", chrom)
        self.assertEqual(198982, pos)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("", seq)
        self.assertEqual("G", ref)

    def test_alternate_chrom(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("G]GL000.01:198982]")
        self.assertEqual("GL000.01", chrom)
        self.assertEqual(198982, pos)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("", seq)
        self.assertEqual("G", ref)

    def test_left_untemp_seq(self):
        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("]11:123456]AGTNNNCAT")
        self.assertEqual("11", chrom)
        self.assertEqual(123456, pos)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("AGTNNNCA", seq)
        self.assertEqual("T", ref)

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("]8:1682443]TGC")
        self.assertEqual("8", chrom)
        self.assertEqual(1682443, pos)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("TG", seq)
        self.assertEqual("C", ref)

        chrom, pos, orient1, orient2, ref, seq = _parse_bnd_alt("AAGTG]11:66289601]")
        self.assertEqual("11", chrom)
        self.assertEqual(66289601, pos)
        self.assertEqual(ORIENT.LEFT, orient2)
        self.assertEqual("AGTG", seq)
        self.assertEqual("A", ref)


class TestBreakDancer(unittest.TestCase):
    def test_itx(self):
        row = {
            "Chr1": "1",
            "Pos1": "10001",
            "Orientation1": "83+126-",
            "Chr2": "1",
            "Pos2": "10546",
            "Orientation2": "83+126-",
            "Type": "ITX",
            "Size": "-352",
            "Score": "99",
            "num_Reads": "43",
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(SVTYPE.DUP, bpps[0].event_type)
        self.assertEqual(10001, bpps[0].break1.start)
        self.assertEqual(10001, bpps[0].break1.end)
        self.assertEqual(ORIENT.RIGHT, bpps[0].break1.orient)
        self.assertEqual(10546, bpps[0].break2.start)
        self.assertEqual(10546, bpps[0].break2.end)
        self.assertEqual(ORIENT.LEFT, bpps[0].break2.orient)
        self.assertEqual(False, bpps[0].opposing_strands)

    def test_deletion(self):
        row = {
            "Chr1": "1",
            "Pos1": "869445",
            "Orientation1": "89+21-",
            "Chr2": "1",
            "Pos2": "870225",
            "Orientation2": "5+93-",
            "Type": "DEL",
            "Size": "892",
            "Score": "99",
            "num_Reads": "67",
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(SVTYPE.DEL, bpps[0].event_type)
        self.assertEqual(869445, bpps[0].break1.start)
        self.assertEqual(869445, bpps[0].break1.end)
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)
        self.assertEqual(870225, bpps[0].break2.start)
        self.assertEqual(870225, bpps[0].break2.end)
        self.assertEqual(ORIENT.RIGHT, bpps[0].break2.orient)
        self.assertEqual(False, bpps[0].opposing_strands)

    def test_inversion(self):
        row = {
            "Chr1": "1",
            "Pos1": "13143396",
            "Orientation1": "18+4-",
            "Chr2": "1",
            "Pos2": "13218683",
            "Orientation2": "10+10-",
            "Type": "INV",
            "Size": "74618",
            "Score": "31",
            "num_Reads": "2",
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        self.assertEqual(2, len(bpps))
        self.assertEqual(SVTYPE.INV, bpps[0].event_type)
        self.assertEqual(13143396, bpps[0].break1.start)
        self.assertEqual(13143396, bpps[0].break1.end)
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)
        self.assertEqual(13218683, bpps[0].break2.start)
        self.assertEqual(13218683, bpps[0].break2.end)
        self.assertEqual(ORIENT.LEFT, bpps[0].break2.orient)
        self.assertEqual(True, bpps[0].opposing_strands)

        self.assertEqual(SVTYPE.INV, bpps[1].event_type)
        self.assertEqual(13143396, bpps[1].break1.start)
        self.assertEqual(13143396, bpps[1].break1.end)
        self.assertEqual(ORIENT.RIGHT, bpps[1].break1.orient)
        self.assertEqual(13218683, bpps[1].break2.start)
        self.assertEqual(13218683, bpps[1].break2.end)
        self.assertEqual(ORIENT.RIGHT, bpps[1].break2.orient)
        self.assertEqual(True, bpps[1].opposing_strands)

    def test_insertion(self):
        row = {
            "Chr1": "1",
            "Pos1": "20216146",
            "Orientation1": "23+26-",
            "Chr2": "1",
            "Pos2": "20218060",
            "Orientation2": "23+26-",
            "Type": "INS",
            "Size": "-421",
            "Score": "99",
            "num_Reads": "3",
        }
        bpps = _convert_tool_row(row, SUPPORTED_TOOL.BREAKDANCER, False, True)
        self.assertEqual(1, len(bpps))
        self.assertEqual(SVTYPE.INS, bpps[0].event_type)
        self.assertEqual(20216146, bpps[0].break1.start)
        self.assertEqual(20216146, bpps[0].break1.end)
        self.assertEqual(ORIENT.LEFT, bpps[0].break1.orient)
        self.assertEqual(20218060, bpps[0].break2.start)
        self.assertEqual(20218060, bpps[0].break2.end)
        self.assertEqual(ORIENT.RIGHT, bpps[0].break2.orient)
        self.assertEqual(False, bpps[0].opposing_strands)


class TestStrelka(unittest.TestCase):
    def testInsertion(self):
        event = Mock(
            chrom="1", pos=724986, id=None, info={}, ref="G", stop=724986, alts=("GGAATT",)
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.STRELKA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(724986, bpp.break1.start)
        self.assertEqual(724986, bpp.break1.end)
        self.assertEqual(724986, bpp.break2.start)
        self.assertEqual(724986, bpp.break2.end)
        self.assertEqual(SVTYPE.INS, bpp.event_type)

    def testDeletion(self):
        event = Mock(
            chrom="1",
            pos=1265353,
            id=None,
            info={},
            ref="GCGTGTGCCATGCA",
            stop=1265366,
            alts=("G",),
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.STRELKA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(1265353, bpp.break1.start)
        self.assertEqual(1265353, bpp.break1.end)
        self.assertEqual(1265366, bpp.break2.start)
        self.assertEqual(1265366, bpp.break2.end)
        self.assertEqual(SVTYPE.DEL, bpp.event_type)


class TestMutect(unittest.TestCase):
    def testInsertion(self):
        event = Mock(
            chrom="1", pos=724986, id=None, info={}, ref="G", stop=724986, alts=("GGAATT",)
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
            chrom="1",
            pos=1265353,
            id=None,
            info={},
            ref="GCGTGTGCCATGCA",
            stop=1265366,
            alts=("G",),
        )
        bpp_list = _convert_tool_row(_parse_vcf_record(event)[0], SUPPORTED_TOOL.MUTECT, False)
        assert len(bpp_list) == 1
        bpp = bpp_list[0]
        assert bpp.break1.start == 1265353
        assert bpp.break1.end == 1265353
        assert bpp.break2.start == 1265366
        assert bpp.break2.end == 1265366
        assert bpp.event_type == SVTYPE.DEL


class TestVCF(unittest.TestCase):
    def setUp(self):
        self.tra = Mock(
            chrom="2",
            pos=21673582,
            id=None,
            info={"SVTYPE": "TRA", "CT": "5to5", "CHR2": "3"},
            stop=58921502,
            alts=[],
        )

    def test_no_ci(self):
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(21673582, bpp.break1.start)
        self.assertEqual(21673582, bpp.break1.end)
        self.assertEqual(58921502, bpp.break2.start)
        self.assertEqual(58921502, bpp.break2.end)

    def test_ci(self):
        self.tra.info.update({"CIEND": [-700, 700], "CIPOS": [-700, 700]})
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        print(bpp)
        self.assertEqual(21673582 - 700, bpp.break1.start)
        self.assertEqual(21673582 + 700, bpp.break1.end)
        self.assertEqual(58921502 - 700, bpp.break2.start)
        self.assertEqual(58921502 + 700, bpp.break2.end)

    def test_precise_flag_ignores_ci(self):
        self.tra.info.update({"CIEND": [-700, 700], "CIPOS": [-700, 700], "PRECISE": True})
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual(21673582, bpp.break1.start)
        self.assertEqual(21673582, bpp.break1.end)
        self.assertEqual(58921502, bpp.break2.start)
        self.assertEqual(58921502, bpp.break2.end)

    def test_no_id(self):
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertTrue(bpp.data[COLUMNS.tracking_id])

    def test_N_id(self):
        self.tra.id = "N"
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertTrue(bpp.data[COLUMNS.tracking_id])
        self.assertNotEqual("N", bpp.data[COLUMNS.tracking_id])

    def test_id_given(self):
        self.tra.id = "thing-1"
        bpp_list = _convert_tool_row(_parse_vcf_record(self.tra)[0], SUPPORTED_TOOL.VCF, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("vcf-thing-1", bpp.data[COLUMNS.tracking_id])


class TestImprecise(unittest.TestCase):
    def setUp(self):
        self.event = Mock(
            chrom="chr5",
            pos=1853407,
            id=None,
            info={"SVTYPE": "INS", "CHR2": "chr5"},
            stop=1853408,
            alts=["AGG"],
            ref="A",
        )

    def test_convert_intrachromosomal_imprecise_breakend1(self):
        self.event.id = "Mock"
        self.event.info.update({"CIPOS": [-30, 30], "CIEND": [-65, 65]})

        bpp_list = _convert_tool_row(_parse_vcf_record(self.event)[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr5", bpp.break1.chr)
        self.assertEqual(1853377, bpp.break1.start)
        self.assertEqual(1853437, bpp.break1.end)
        self.assertEqual(1853377, bpp.break2.start)
        self.assertEqual(1853473, bpp.break2.end)
        self.assertEqual("chr5", bpp.break2.chr)

    def test_convert_intrachromosomal_imprecise_breakend2(self):
        self.event.id = "Mock"
        self.event.info.update({"CIPOS": [-30, 99999], "CIEND": [-10, 65]})

        bpp_list = _convert_tool_row(_parse_vcf_record(self.event)[0], SUPPORTED_TOOL.MANTA, False)
        self.assertEqual(1, len(bpp_list))
        bpp = bpp_list[0]
        self.assertEqual("chr5", bpp.break1.chr)
        self.assertEqual(1853377, bpp.break1.start)
        self.assertEqual(1853473, bpp.break1.end)
        self.assertEqual(1853398, bpp.break2.start)
        self.assertEqual(1853473, bpp.break2.end)
        self.assertEqual("chr5", bpp.break2.chr)

    def test_error_on_convert_intrachromosomal_imprecise_breakend(self):
        self.event.id = "Mock"
        self.event.info.update({"CIPOS": [-30, 99999], "CIEND": [70, 65]})
        with self.assertRaises(ValueError):
            _convert_tool_row(_parse_vcf_record(self.event)[0], SUPPORTED_TOOL.MANTA, False)

    def test_error_on_convert_intrachromosomal_imprecise_breakend2(self):
        # breakpoint_1_start > breakpoint_1_end
        variant_cilen4 = Mock(
            id="Sniffle.INS",
            pos=11184,
            chrom="chr2",
            alts=["AGG"],
            ref="N",
            stop=11183,
            info=dict(
                CHR2="chr2",
                IMPRECISE=True,
                SVTYPE="INS",
                END=11183,
            ),
        )
        with self.assertRaises(ValueError):
            bpp_list = _convert_tool_row(
                _parse_vcf_record(variant_cilen4)[0], SUPPORTED_TOOL.MANTA, False
            )
            bpp = bpp_list[0]
            import pdb

            pdb.set_trace()
            print(bpp)
