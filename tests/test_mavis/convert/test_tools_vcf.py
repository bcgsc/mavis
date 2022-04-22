import pytest
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


TEST_POS = 1853407


@pytest.mark.parametrize(
    'pos,break1_ci,break2_ci,break1,break2',
    [
        [
            TEST_POS,
            (-30, 30),
            (-65, 65),
            (TEST_POS - 30, TEST_POS + 30),
            (TEST_POS - 30, TEST_POS + 65),
        ],
        [
            TEST_POS,
            (-30, 99999),
            (-10, 65),
            (TEST_POS - 30, TEST_POS + 65),
            (TEST_POS - 10, TEST_POS + 65),
        ],
        [
            TEST_POS,
            (-30, 99999),
            (70, 65),
            (TEST_POS - 30, TEST_POS + 65),
            (TEST_POS + 65, TEST_POS + 65),
        ],
    ],
    ids=[
        'breakpoint_2_start < breakpoint_1_start',
        'breakpoint_1_end > breakpoint_2_end',
        'breakpoint_2_start > breakpoint_2_end',
    ],
)
def test_convert_intrachromosomal_imprecise_breakend(pos, break1_ci, break2_ci, break1, break2):
    variant_vcf = VcfRecordType(
        id='vcf-cuteSV.INS',
        pos=pos,
        chrom='chr5',
        alts=['AGG'],
        ref='A',
        info=VcfInfoType(
            CHR2="chr5",
            IMPRECISE=True,
            SVMETHOD="cuteSV-1.0.12",
            SVTYPE="INS",
            CIPOS=break1_ci,
            CILEN=break2_ci,
        ),
    )
    result = convert_record(variant_vcf)
    assert len(result) == 1
    variant = result[0]
    assert variant.get('break1_position_start') == break1[0]
    assert variant.get('break1_position_end') == break1[1]
    assert variant.get('break2_position_start') == break2[0]
    assert variant.get('break2_position_end') == break2[1]


def test_convert_intrachromosomal_imprecise_breakend_no_ci():
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
