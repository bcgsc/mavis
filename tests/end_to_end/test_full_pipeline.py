import tempfile
import re
import shutil
import unittest
import os
import subprocess
import glob


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


@unittest.skipIf(
    not int(os.environ.get('RUN_FULL', 1)),
    'slower tests will not be run unless the environment variable RUN_FULL is given')
class TestFullPipeline(unittest.TestCase):

    def test_mocked(self):
        command = 'mavis pipeline {} -o {}'.format(config, temp_output)
        print(command)
        try:
            output = subprocess.check_output(command, shell=True, env=os.environ)
        except subprocess.CalledProcessError as err:
            print('failed command', command)
            print(err.output.decode('UTF-8'))
            raise err

        # check that the subdirectories were built
        for lib in[mock_genome + '_*', mock_trans + '_*']:
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
            try:
                output = subprocess.check_output(command, shell=True, env=os.environ)
            except subprocess.CalledProcessError as err:
                print('failed command', command)
                print(err.output.decode('UTF-8'))
                raise err

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
                self.assertTrue(glob_exists(temp_output, lib, 'validation/*-1', '*.' + suffix))
            self.assertTrue(glob_exists(temp_output, lib, 'validation/*-1', '*.COMPLETE'))

            # run annotation
            self.assertTrue(glob_exists(temp_output, lib, 'annotation'))
            qsub = os.path.join(temp_output, lib, 'annotation', 'qsub.sh')
            self.assertTrue(glob_exists(qsub))
            command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
            try:
                output = subprocess.check_output(command, shell=True, env=os.environ)
            except subprocess.CalledProcessError as err:
                print('failed command', command)
                print(err.output.decode('UTF-8'))
                raise err
            # check the generated files
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', 'annotations.tab'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', 'annotations.fusion-cdna.fa'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', 'drawings'))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', 'drawings', '*svg', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', 'drawings', '*json', strict=False))
            self.assertTrue(glob_exists(temp_output, lib, 'annotation/*-1', '*.COMPLETE'))
        # now run the pairing
        self.assertTrue(glob_exists(temp_output, 'pairing'))
        qsub = os.path.join(temp_output, 'pairing', 'qsub.sh')
        self.assertTrue(glob_exists(qsub))
        command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
        try:
            output = subprocess.check_output(command, shell=True, env=os.environ)
        except subprocess.CalledProcessError as err:
            print('failed command', command)
            print(err.output.decode('UTF-8'))
            raise err

        self.assertTrue(glob_exists(temp_output, 'pairing', 'mavis_paired*.tab'))
        self.assertTrue(glob_exists(temp_output, 'pairing', '*.COMPLETE'))

        # now run the summary
        self.assertTrue(glob_exists(temp_output, 'summary'))
        qsub = os.path.join(temp_output, 'summary', 'qsub.sh')
        self.assertTrue(glob_exists(qsub))
        command = 'export SGE_TASK_ID=1; bash {}'.format(qsub)
        try:
            output = subprocess.check_output(command, shell=True, env=os.environ)
        except subprocess.CalledProcessError as err:
            print('failed command', command)
            print(err.output.decode('UTF-8'))
            raise err

        self.assertTrue(glob_exists(temp_output, 'summary', 'mavis_summary*.tab'))
        self.assertTrue(glob_exists(temp_output, 'summary', '*.COMPLETE'))


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(temp_output)
    pass

if __name__ == "__main__":
    unittest.main()
