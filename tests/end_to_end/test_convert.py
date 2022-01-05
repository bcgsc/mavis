import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

from mavis.constants import ORIENT, SUBCOMMAND, SVTYPE
from mavis.main import main
from mavis.tools import SUPPORTED_TOOL
from mavis.util import read_bpp_from_input_file

from ..util import get_data, glob_exists

TEMP_OUTPUT = None


def setUpModule():
    global TEMP_OUTPUT
    # create the temp output directory to store file outputs
    TEMP_OUTPUT = tempfile.mkdtemp()
    print('output dir', TEMP_OUTPUT)


class TestConvert:
    def run_main(self, inputfile, file_type, strand_specific=False):
        outputfile = os.path.join(TEMP_OUTPUT, file_type + '.tab')
        args = [
            'mavis',
            SUBCOMMAND.CONVERT,
            '-o',
            outputfile,
            '-n',
            inputfile,
            '--file_type',
            file_type,
            '--strand_specific',
            strand_specific,
        ]
        with patch.object(sys, 'argv', args):
            main()
            print('output', outputfile)
            assert glob_exists(outputfile, n=1)
        result = {}
        for pair in read_bpp_from_input_file(outputfile):
            result.setdefault(pair.data['tracking_id'], []).append(pair)
        return result

    def test_chimerascan(self):
        self.run_main(get_data('chimerascan_output.bedpe'), SUPPORTED_TOOL.CHIMERASCAN, False)

    def test_defuse(self):
        self.run_main(get_data('defuse_output.tsv'), SUPPORTED_TOOL.DEFUSE, False)

    def test_delly(self):
        result = self.run_main(get_data('delly_events.vcf'), SUPPORTED_TOOL.DELLY, False)
        # test the contents were converted successfully
        assert len(result['delly-DUP00000424']) == 1
        bpp = result['delly-DUP00000424'][0]
        print(bpp.data)
        print(bpp)
        assert bpp.event_type == SVTYPE.DUP
        assert bpp.break1.chr == '1'
        assert bpp.break2.chr == '1'
        assert bpp.break1.start == 224646569
        assert bpp.break1.end == 224646569
        assert bpp.break2.start == 224800120
        assert bpp.break2.end == 224800120
        assert len(result['delly-TRA00020624']) == 1
        bpp = result['delly-TRA00020624'][0]
        assert bpp.event_type == SVTYPE.TRANS
        assert bpp.break1.chr == '10'
        assert bpp.break2.chr == '19'
        assert bpp.break1.start == 7059510 - 670
        assert bpp.break1.end == 7059510 + 670
        assert bpp.break2.start == 17396810 - 670
        assert bpp.break2.end == 17396810 + 670
        assert 31 == len(result)

    def test_manta(self):
        result = self.run_main(get_data('manta_events.vcf'), SUPPORTED_TOOL.MANTA, False)
        # ensure weird bnd type is converted correctly
        bnd_id = 'manta-MantaBND:173633:0:1:0:0:0:0'
        assert len(result[bnd_id]) == 1
        bpp = result[bnd_id][0]
        assert bpp.event_type == SVTYPE.TRANS
        assert bpp.break1.chr == '10'
        assert bpp.break2.chr == '19'
        assert bpp.break1.start == 7059511 - 0
        assert bpp.break1.end == 7059511 + 1
        assert bpp.break2.start == 17396810
        assert bpp.break2.end == 17396810
        assert bpp.break2.orient == ORIENT.LEFT
        somatic_event = result['manta-MantaDEL:20644:0:2:0:0:0'][0]
        assert somatic_event.data.get('SOMATIC', False) is True

    def test_pindel(self):
        self.run_main(get_data('pindel_events.vcf'), SUPPORTED_TOOL.PINDEL, False)

    def test_transabyss(self):
        self.run_main(get_data('transabyss_indels_output.tab'), SUPPORTED_TOOL.TA, False)
        self.run_main(get_data('transabyss_events.tab'), SUPPORTED_TOOL.TA, False)

    def test_vcf(self):
        results = self.run_main(get_data('clinvar_short_test.vcf'), SUPPORTED_TOOL.VCF, False)
        print(results.keys())
        record = results['vcf-460818'][0]
        print(record, record.data)
        assert record.data['CLNSIG'] == 'Pathogenic'

    def test_sniffle(self):
        results = self.run_main(get_data('sniffles.vcf'), SUPPORTED_TOOL.VCF, False)
        print(results.keys())
        record = results['vcf-35777'][0]
        print(record, record.data)
        assert record.data['event_type'] == 'translocation'
    
    def test_cuteSV(self):
        results = self.run_main(get_data('cuteSV.vcf'), SUPPORTED_TOOL.VCF, False)
        print(results.keys())
        record = results['vcf-cuteSV.BND.0'][0]
        print(record, record.data)
        assert record.data['event_type'] == 'inverted translocation'
    
    def test_breakseq2(self):
        self.run_main(get_data('breakseq.vcf'), SUPPORTED_TOOL.BREAKSEQ, False)

    def test_cnvnator(self):
        self.run_main(get_data('cnvnator.tab'), SUPPORTED_TOOL.CNVNATOR, False)

    def test_breakdancer(self):
        self.run_main(get_data('breakdancer_output.txt'), SUPPORTED_TOOL.BREAKDANCER, False)

    def test_strelka(self):
        self.run_main(get_data('strelka.vcf'), SUPPORTED_TOOL.STRELKA, False)


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(TEMP_OUTPUT)


if __name__ == '__main__':
    unittest.main()
