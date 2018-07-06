import glob
import os
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch


from mavis.constants import ORIENT, SUBCOMMAND, SVTYPE
from mavis.main import main
from mavis.tools import SUPPORTED_TOOL
from mavis.util import unique_exists, read_bpp_from_input_file

from ..util import get_data


TEMP_OUTPUT = None


def setUpModule():
    global TEMP_OUTPUT
    # create the temp output directory to store file outputs
    TEMP_OUTPUT = tempfile.mkdtemp()
    print('output dir', TEMP_OUTPUT)


class TestConvert(unittest.TestCase):

    def run_main(self, inputfile, file_type, strand_specific=False):
        outputfile = os.path.join(TEMP_OUTPUT, file_type + '.tab')
        args = [
            'mavis', SUBCOMMAND.CONVERT,
            '-o', outputfile,
            '-n', inputfile,
            '--file_type', file_type,
            '--strand_specific', strand_specific
        ]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())
            print('output', outputfile)
            self.assertTrue(unique_exists(outputfile))
        result = {}
        for pair in read_bpp_from_input_file(outputfile):
            result.setdefault(pair.tracking_id, []).append(pair)
        return result

    def test_chimerascan(self):
        self.run_main(get_data('chimerascan_output.bedpe'), SUPPORTED_TOOL.CHIMERASCAN, False)

    def test_defuse(self):
        self.run_main(get_data('defuse_output.tsv'), SUPPORTED_TOOL.DEFUSE, False)

    def test_delly(self):
        result = self.run_main(get_data('delly_events.vcf'), SUPPORTED_TOOL.DELLY, False)
        # test the contents were converted successfully
        self.assertEqual(1, len(result['delly-DUP00000424']))
        bpp = result['delly-DUP00000424'][0]
        self.assertEqual(SVTYPE.DUP, bpp.event_type)
        self.assertEqual('1', bpp.break1.chr)
        self.assertEqual('1', bpp.break2.chr)
        self.assertEqual(224646569, bpp.break1.start)
        self.assertEqual(224646569, bpp.break1.end)
        self.assertEqual(224800120, bpp.break2.start)
        self.assertEqual(224800120, bpp.break2.end)
        self.assertEqual(1, len(result['delly-TRA00020624']))
        bpp = result['delly-TRA00020624'][0]
        self.assertEqual(SVTYPE.TRANS, bpp.event_type)
        self.assertEqual('10', bpp.break1.chr)
        self.assertEqual('19', bpp.break2.chr)
        self.assertEqual(7059510 - 670, bpp.break1.start)
        self.assertEqual(7059510 + 670, bpp.break1.end)
        self.assertEqual(17396810 - 670, bpp.break2.start)
        self.assertEqual(17396810 + 670, bpp.break2.end)
        self.assertEqual(len(result), 31)

    def test_manta(self):
        result = self.run_main(get_data('manta_events.vcf'), SUPPORTED_TOOL.MANTA, False)
        # ensure weird bnd type is converted correctly
        bnd_id = 'manta-MantaBND:173633:0:1:0:0:0:0'
        self.assertEqual(2, len(result[bnd_id]))
        bpp1, bpp2 = result[bnd_id]
        if bpp1.event_type == SVTYPE.ITRANS:
            bpp2, bpp1 = bpp1, bpp2
        self.assertEqual(SVTYPE.TRANS, bpp1.event_type)
        self.assertEqual(SVTYPE.ITRANS, bpp2.event_type)
        for bpp in [bpp1, bpp2]:
            self.assertEqual('10', bpp.break1.chr)
            self.assertEqual('19', bpp.break2.chr)
            self.assertEqual(7059511 - 0, bpp.break1.start)
            self.assertEqual(7059511 + 1, bpp.break1.end)
            self.assertEqual(17396810, bpp.break2.start)
            self.assertEqual(17396810, bpp.break2.end)
            self.assertEqual(ORIENT.LEFT, bpp.break2.orient)
        somatic_event = result['manta-MantaDEL:20644:0:2:0:0:0'][0]
        self.assertEqual('True', somatic_event.data.get('SOMATIC', False))

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
        self.assertEqual('Pathogenic', record.data['CLNSIG'])

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
