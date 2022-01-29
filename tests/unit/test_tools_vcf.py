from mavis.tools import SUPPORTED_TOOL, _convert_tool_row
from mavis.tools.vcf import VcfInfoType, VcfRecordType, convert_record, pandas_vcf

from ..util import get_data


def test_read_vcf():
    header, df = pandas_vcf(get_data('sniffles.vcf'))
    assert len(header) == 231
    assert df.shape[0] == 106


def test_convert_record():
    variant = VcfRecordType(
        1,
        0,
        'chr14_KI270722v1_random',
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
    records = convert_record(variant)
    assert len(records) == 1
    record = records[0]
    assert record.get('break1_position_start') == 1
    assert record.get('break1_position_end') == 1
    assert record.get('break2_position_start') == 1
    assert record.get('break2_position_end') == 1
    assert record.get('break2_chromosome') == 'chr17_GL000205v2_random'
