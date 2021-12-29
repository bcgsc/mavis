from mavis.tools.vcf import pandas_vcf

from ..util import get_data


def test_read_vcf():
    header, df = pandas_vcf(get_data('delly_events.vcf'))
    assert len(header) == 63
    assert df.shape[0] == 31
