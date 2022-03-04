from mavis.convert import SUPPORTED_TOOL, _convert_tool_row
from mavis.convert.vcf import VcfInfoType, VcfRecordType, convert_record, pandas_vcf

from ...util import get_data


def test_read_vcf():
    header, df = pandas_vcf(get_data('sniffles.vcf'))
    assert len(header) == 231
    assert df.shape[0] == 106


def test_convert_record():
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


def test_convert_record_cuteSV():
    variant_cilen = VcfRecordType(
        id='vcf-cuteSV.INS',
        pos=1853407,
        chrom='chr5',
        alts=[
            'AGGATCTATGTGGCTGTTGCAGGGTGACCCGAGGTCACGAGAGGCAAGGTCAGAGGACGATGTGAGGGCTGCAGGGTGACCCGAGGTCACGTAGGGCAAGGTCAGAGGACGATGTGGCGGTTGCAGGGAGACCCAGGTCACGCAGGCAAGGTCAGAGGACGATGTGAGGGAGTTGCAGGGTGACCCGAGGTCACGTAGGGCAAGGTCAGAGGACGATGTGGCGGTTGCAGGGTGACCCGAGGTCA'
        ],
        ref='A',
        info=VcfInfoType(
            CHR2="chr5",
            IMPRECISE=True,
            SVMETHOD="cuteSV-1.0.12",
            SVTYPE="INS",
            SUPTYPE="None",
            STRANDS="None",
            CIPOS=(-30, 30),
            CILEN=(-65, 65),
        ),
    )
    variant_ins_cutesv = convert_record(variant_cilen)
    assert len(variant_ins_cutesv) == 1
    variant_ins_cutesv = variant_ins_cutesv[0]
    assert variant_ins_cutesv.get('break2_position_start') == 1853342
    assert variant_ins_cutesv.get('break2_position_end') == 1853472
