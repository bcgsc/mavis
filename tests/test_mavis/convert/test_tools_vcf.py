from mavis.convert import SUPPORTED_TOOL, _convert_tool_row
from mavis.convert.vcf import VcfInfoType, VcfRecordType, convert_record, pandas_vcf

from ...util import get_data


def test_read_vcf():
    header, df = pandas_vcf(get_data('sniffles.vcf'))
    assert len(header) == 231
    assert df.shape[0] == 106


def test_convert_telomeric_region():
    variant_imprecise = VcfRecordType(
        id='mock-BND-imprecise',
        pos=0,
        chrom='chr14_KI270722v1_random',
        alts=['N[chr17_GL000205v2_random:0['],
        ref='N',
        info=VcfInfoType(
            IMPRECISE=True,
            SVMETHOD="Snifflesv1.0.11",
            SVTYPE="BND",
            SUPTYPE="SR",
            SVLEN="0",
            STRANDS="+-",
            RE="5",
            REF_strand="0,0",
            AF="1",
        ),
    )
    variant_precise = VcfRecordType(
        id='mock-BND-precise',
        pos=0,
        chrom='chr14_KI270722v1_random',
        alts=[']chrUn_GL000216v2:142821]N'],
        ref='N',
        info=VcfInfoType(
            IMPRECISE=False,
            SVMETHOD="Snifflesv1.0.11",
            SVTYPE="BND",
            SUPTYPE="SR",
            SVLEN="0",
            STRANDS="+-",
            RE="5",
            REF_strand="0,0",
            AF="1",
        ),
    )
    imprecise_records = convert_record(variant_imprecise)
    assert len(imprecise_records) == 1
    imprecise_records = imprecise_records[0]
    assert imprecise_records.get('break1_position_end') == 1

    precise_records = convert_record(variant_precise)
    assert len(precise_records) == 1
    precise_records = precise_records[0]
    assert precise_records.get('break1_position_end') == 1

    assert precise_records.get('break1_chromosome') == 'chr14_KI270722v1_random'
    assert imprecise_records.get('break1_chromosome') == 'chr14_KI270722v1_random'


def test_convert_intrachromosomal_imprecise_breakend():
    # breakpoint_2_start < breakpoint_1_start
    variant_imprecise = VcfRecordType(
        id='vcf-cuteSV.INS',
        pos=1853407,
        chrom='chr5',
        alts=['AGG'],
        ref='A',
        info=VcfInfoType(
            CHR2="chr5",
            IMPRECISE=True,
            SVMETHOD="cuteSV-1.0.12",
            SVTYPE="INS",
            CIPOS=(-30, 30),
            CILEN=(-65, 65),
        ),
    )
    variant_ins_imprecise = convert_record(variant_imprecise)
    assert len(variant_ins_imprecise) == 1
    variant_ins_imprecise = variant_ins_imprecise[0]
    assert variant_ins_imprecise.get('break2_position_start') == 1853377
    assert variant_ins_imprecise.get('break2_position_end') == 1853472

    # breakpoint_1_end > breakpoint_2_end
    variant_imprecise2 = VcfRecordType(
        id='vcf-cuteSV.INS',
        pos=1853407,
        chrom='chr5',
        alts=['AGG'],
        ref='A',
        info=VcfInfoType(
            CHR2="chr5",
            IMPRECISE=True,
            SVMETHOD="cuteSV-1.0.12",
            SVTYPE="INS",
            SUPTYPE="None",
            STRANDS="None",
            CIPOS=(-30, 9999),
            CILEN=(-65, 65),
        ),
    )
    variant_ins_imprecise_2 = convert_record(variant_imprecise2)
    assert len(variant_ins_imprecise_2) == 1
    variant_ins_imprecise_2 = variant_ins_imprecise_2[0]
    assert variant_ins_imprecise_2.get('break1_position_start') == 1853377
    assert variant_ins_imprecise_2.get('break1_position_end') == 1853472
    assert variant_ins_imprecise_2.get('break2_position_end') == 1853472
    assert variant_ins_imprecise_2.get('break2_position_start') == 1853377

    # breakpoint_2_start > breakpoint_2_end
    variant_imprecise3 = VcfRecordType(
        id='mock-INS-imprecise',
        pos=1853407,
        chrom='chr5',
        alts=['AGG'],
        ref='A',
        info=VcfInfoType(
            CHR2="chr5",
            IMPRECISE=True,
            SVMETHOD="Snifflesv1.0.11",
            SVTYPE="INS",
            SUPTYPE="None",
            STRANDS="None",
            CIPOS=(-30, 9999),
            CILEN=(70, 65),
        ),
    )
    variant_ins_imprecise_3 = convert_record(variant_imprecise3)
    assert len(variant_ins_imprecise_3) == 1
    variant_ins_imprecise_3 = variant_ins_imprecise_3[0]
    assert variant_ins_imprecise_3.get('break2_position_end') == 1853472
    assert variant_ins_imprecise_3.get('break2_position_start') == 1853472

    # breakpoint_1_start > breakpoint_1_end
    variant_cilen4 = VcfRecordType(
        id='Sniffle.INS',
        pos=11184,
        chrom='chr2',
        alts=['AGG'],
        ref='N',
        info=VcfInfoType(
            CHR2="chr2",
            IMPRECISE=True,
            SVTYPE="INS",
            END=11183,
        ),
    )
    variant_ins_imprecise_4 = convert_record(variant_cilen4)
    assert len(variant_ins_imprecise_4) == 1
    variant_ins_imprecise_4 = variant_ins_imprecise_4[0]
    assert variant_ins_imprecise_4.get('break2_position_end') == 11183
    assert variant_ins_imprecise_4.get('break1_position_end') == 11183
