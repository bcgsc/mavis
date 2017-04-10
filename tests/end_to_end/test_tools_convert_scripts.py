import tempfile
import shutil
import unittest
import os
import subprocess


script_prefix = os.path.join(os.path.dirname(__file__), './../../tools')
output_checker = os.path.join(script_prefix, 'test_input_file.py')
data_prefix = os.path.join(os.path.dirname(__file__), 'data')
temp_output = None


def setUpModule():
    global temp_output
    # create the temp output directory to store file outputs
    temp_output = tempfile.mkdtemp()


class TestConvertScripts(unittest.TestCase):
    def test_chimerascan(self):
        raise unittest.SkipTest('TODO')

    def test_delly(self):
        input_file = os.path.join(data_prefix, 'delly_events.vcf')
        script = os.path.join(script_prefix, 'convert_delly.py')
        output = os.path.join(temp_output, 'delly_out.tab')
        subprocess.check_output('python {} -f {} -o {}'.format(script, input_file, output), shell=True)
        subprocess.check_output(['python', output_checker, output]) 

    def test_transabyss(self):
        input_file = os.path.join(data_prefix, 'transabyss_events.tab')
        script = os.path.join(script_prefix, 'convert_transabyss.py')
        output = os.path.join(temp_output, 'transabyss_out.tab')
        subprocess.check_output(
            'python {} -n {} -o {} -l testlibrary -p genome'.format(script, input_file, output), shell=True)
        subprocess.check_output(['python', output_checker, output]) 

    def test_manta(self):
        input_file = os.path.join(data_prefix, 'manta_events.vcf')
        script = os.path.join(script_prefix, 'convert_manta.py')
        output = os.path.join(temp_output, 'manta_out.tab')
        subprocess.check_output('python {} -i {} -o {} -l testlibrary'.format(script, input_file, output), shell=True)
        subprocess.check_output(['python', output_checker, output])

    def test_pindel(self):
        raise unittest.SkipTest('TODO')
    
    def test_defuse(self):
        raise unittest.SkipTest('TODO')

    def test_breakdancer(self):
        raise unittest.SkipTest('TODO')


def tearDownModule():
    # remove the temp directory and outputs
    shutil.rmtree(temp_output)

if __name__ == "__main__":
    unittest.main()
