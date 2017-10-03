import os
import subprocess
import sys
import unittest


from mavis.constants import PIPELINE_STEP
from mavis.main import main
from mock import patch


class TestHelpMenu(unittest.TestCase):

    def test_main(self):
        with patch.object(sys, 'argv', ['mavis', '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_pipeline(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.PIPELINE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_config(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.CONFIG, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_cluster(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.CLUSTER, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_validate(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.VALIDATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_annotate(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.ANNOTATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_pairing(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.PAIR, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_summary(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.SUMMARY, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_bad_option(self):
        with patch.object(sys, 'argv', ['mavis', PIPELINE_STEP.PIPELINE, '--blargh']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertNotEqual(0, err.code)
            else:
                self.assertNotEqual(0, returncode)
