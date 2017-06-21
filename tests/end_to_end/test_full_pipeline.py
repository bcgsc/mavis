import tempfile
import shutil
import unittest
import os
import subprocess
import glob
from mavis.constants import PROTOCOL


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
    l = len(glob.glob(globexpr))
    if strict and l == 1:
        return True
    elif not strict and l > 0:
        return True
    else:
        print(globexpr)
        print(l)
        return False


def tail(output, count=10):
    for l in output.decode('UTF8').split('\n')[-1 * count:]:
        print(l)


@unittest.skipIf(
    not int(os.environ.get('RUN_FULL', 1)),
    'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestFullPipeline(unittest.TestCase):

    def test_mocked(self):
        command = 'python {} pipeline {} -o {}'.format(main_run_script, config, temp_output)
        print(command)
        output = subprocess.check_output(command, shell=True)
        tail(output)

        # check that the subdirectories were built
        for lib in[mock_genome + '_' + PROTOCOL.GENOME, mock_trans + '_' + PROTOCOL.TRANS]:
            self.assertTrue(glob_exists(temp_output, lib, 'clustering'))
            self.assertTrue(glob_exists(temp_output, lib, 'clustering', 'batch-*-1.tab'))
            self.assertTrue(glob_exists(temp_output, lib, 'clustering', 'uninformative_clusters.txt'))
            self.assertTrue(glob_exists(temp_output, lib, 'clustering', 'clusters.bed'))
            self.assertTrue(glob_exists(temp_output, lib, 'clustering', 'cluster_assignment.tab'))
            self.assertTrue(glob_exists(temp_output, lib, 'clustering', '*.COMPLETE'))

            # run validation
            self.assertTrue(glob_exists(temp_output, lib, 'validation'))
            qsub = os.path.join(temp_output, lib, 'validation', 'qsub.sh')
            self.assertTrue(glob_exists(qsub))
            command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
            print(command)
            output = subprocess.check_output(command, shell=True)
            tail(output)

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
                self.assertTrue(glob_exists(temp_output, lib, 'validation', '*.' + suffix))
            self.assertTrue(glob_exists(temp_output, lib, 'validation', '*.COMPLETE'))

            # run annotation
            self.assertTrue(glob_exists(temp_output, lib, 'annotation'))
            qsub = os.path.join(temp_output, lib, 'annotation', 'qsub.sh')
            self.assertTrue(glob_exists(qsub))
            command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
            print(command)
            output = subprocess.check_output(command, shell=True)
            tail(output)
            # check the generated files
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', 'annotations.tab'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', 'annotations.fusion-cdna.fa'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', 'drawings'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', 'drawings', '*svg', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', 'drawings', '*json', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*', '*.COMPLETE'))
        # now run the pairing
        self.assertTrue(glob_exists(temp_output, 'pairing'))
        qsub = os.path.join(temp_output, 'pairing', 'qsub.sh')
        self.assertTrue(glob_exists(qsub))
        command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
        print(command)
        output = subprocess.check_output(command, shell=True)
        tail(output)

        self.assertTrue(glob_exists(temp_output, 'pairing', 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(temp_output, 'pairing', '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(temp_output, 'summary'))
        qsub = os.path.join(temp_output, 'summary', 'qsub.sh')
        self.assertTrue(glob_exists(qsub))
        command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
        print(command)
        output = subprocess.check_output(command, shell=True)
        tail(output)

        self.assertTrue(glob_exists(temp_output, 'summary', 'mavis_summary*.tab'))
        self.assertTrue(glob_exists(temp_output, 'summary', '*.COMPLETE'))


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(temp_output)
    pass

if __name__ == "__main__":
    unittest.main()
