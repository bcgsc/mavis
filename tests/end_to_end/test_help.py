import unittest
import os
import subprocess


class TestHelpMenu(unittest.TestCase):

    def test_main(self):
        subprocess.check_output(['mavis', '-h'], env=os.environ)

    def test_pipeline(self):
        subprocess.check_output(['mavis', 'pipeline', '-h'], env=os.environ)

    def test_config(self):
        subprocess.check_output(['mavis', 'config', '-h'], env=os.environ)

    def test_cluster(self):
        subprocess.check_output(['mavis', 'cluster', '-h'], env=os.environ)

    def test_validate(self):
        subprocess.check_output(['mavis', 'validate', '-h'], env=os.environ)

    def test_annotate(self):
        subprocess.check_output(['mavis', 'annotate', '-h'], env=os.environ)

    def test_pairing(self):
        subprocess.check_output(['mavis', 'pairing', '-h'], env=os.environ)

    def test_summary(self):
        subprocess.check_output(['mavis', 'summary', '-h'], env=os.environ)
