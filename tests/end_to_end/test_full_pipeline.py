import glob
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


from mavis.constants import PIPELINE_STEP
from mavis.main import main
from mavis.util import unique_exists
from mock import patch


data_prefix = os.path.join(os.path.dirname(__file__), 'data')
temp_output = None
mavis_home = os.path.join(os.path.dirname(__file__), './../..')
main_run_script = os.path.join(mavis_home, 'bin/mavis_run.py')
config = os.path.join(data_prefix, 'pipeline_config.cfg')
mock_genome = 'mock-A36971'
mock_trans = 'mock-A47933'


def setUpModule():
    global temp_output
    # create the temp output directory to store file outputs
    temp_output = tempfile.mkdtemp()
    print('output dir', temp_output)


def glob_exists(*pos, strict=True):
    globexpr = os.path.join(*pos)
    l = glob.glob(globexpr)
    if strict and len(l) == 1:
        return l[0]
    elif not strict and len(l) > 0:
        return l
    else:
        print(globexpr)
        print(l)
        return False


def tail(output, count=50):
    output = output.decode('UTF8').split('\n')
    for l in output[max(len(output) - count, 0):]:
        print(l)


def convert_qsub_to_args(filename, sub=None):
    with open(filename, 'r') as fh:
        lines = [l.strip() for l in fh.readlines() if not l.startswith('#') and l]
        lines = ' '.join(lines)
        if sub:
            for original, replacement in sub:
                lines.replace(original, replacement)
        lex = shlex.split(lines)
        return [a.strip() for a in lex]



@unittest.skipIf(
    not int(os.environ.get('RUN_FULL', 1)),
    'slower tests will not be run unless the environment variable RUN_FULL is given')

class TestFullPipeline(unittest.TestCase):

    def test_full_pipeline(self):
        args = ['mavis', PIPELINE_STEP.PIPELINE, config, '-o', temp_output]
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())
        
        # check that the subdirectories were built
        for lib in[mock_genome + '_*', mock_trans + '_*']:
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER, 'batch-*-1.tab'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER, 'uninformative_clusters.txt'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER, 'clusters.bed'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER, 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.CLUSTER, '*.COMPLETE'))

            # run validation
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.VALIDATE))
            qsub = unique_exists(os.path.join(temp_output, lib, PIPELINE_STEP.VALIDATE, 'qsub.sh'))
            args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])  # read the arguments from the file
            print(args)
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())

            for suffix in [
                'contigs.bam',
                'contigs.blat_out.pslx',
                'contigs.fa',
                'contigs.sorted.bam',
                'contigs.sorted.bam.bai',
                'evidence.bed',
                'igv.batch',
                'raw_evidence.bam',
                'raw_evidence.sorted.bam',
                'raw_evidence.sorted.bam.bai',
                'validation-failed.tab',
                'validation-passed.tab'
            ]:
                self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.VALIDATE + '/*-1', '*.' + suffix))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.VALIDATE + '/*-1', '*.COMPLETE'))

            # run annotation
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE))
            qsub = unique_exists(os.path.join(temp_output, lib, PIPELINE_STEP.ANNOTATE, 'qsub.sh'))
            args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
            with patch.object(sys, 'argv', args):
                self.assertEqual(0, main())
            # check the generated files
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', 'annotations.tab'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', 'annotations.fusion-cdna.fa'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', 'drawings'))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', 'drawings', '*svg', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', 'drawings', '*json', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, PIPELINE_STEP.ANNOTATE + '/*-1', '*.COMPLETE'))
        # now run the pairing
        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.PAIR))
        qsub = unique_exists(os.path.join(temp_output, PIPELINE_STEP.PAIR, 'qsub.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.PAIR, 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.PAIR, '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.SUMMARY))
        qsub = unique_exists(os.path.join(temp_output, PIPELINE_STEP.SUMMARY, 'qsub.sh'))
        args = convert_qsub_to_args(qsub, sub=[('SGE_TASK_ID', '1')])
        with patch.object(sys, 'argv', args):
            self.assertEqual(0, main())

        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.SUMMARY, 'mavis_summary*.tab'))
        self.assertTrue(glob_exists(temp_output, PIPELINE_STEP.SUMMARY, '*.COMPLETE'))

        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.CHECKER, '-o', temp_output]):
            self.assertEqual(0, main())


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(temp_output)
    pass


if __name__ == '__main__':
    unittest.main()
