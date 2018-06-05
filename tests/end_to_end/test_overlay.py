import os
import shutil
import subprocess
import sys
import tempfile
import unittest

from unittest.mock import patch

from mavis.constants import SUBCOMMAND
from mavis.main import main

from . import glob_exists
from ..util import get_data


ANNOTATIONS = get_data('annotations_subsample.json')
BAM = get_data('mock_reads_for_events.sorted.bam')


class TestOverlayOptions(unittest.TestCase):
    def setUp(self):
        # create the temp output directory to store file outputs
        self.temp_output = tempfile.mkdtemp()
        print('output dir', self.temp_output)

    def test_basic(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output
        ]):
            try:
                print(sys.argv)
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_marker(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--marker', 'm', '49364900'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_marker_range(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--marker', 'm', '49364900', '49365900'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_marker_not_enough_args(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--marker', 'm'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertNotEqual(0, err.code)
            else:
                self.assertNotEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, '--annotations', ANNOTATIONS, 'GAGE4', '--output', self.temp_output,
            '--marker'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertNotEqual(0, err.code)
            else:
                self.assertNotEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_marker_not_int(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--marker', 'm', 'k'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertNotEqual(0, err.code)
            else:
                self.assertNotEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_read_depth_plot(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--read_depth_plot', 'axis', BAM
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_read_depth_plot_binned(self):
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--read_depth_plot', 'axis', BAM, '0.5'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def test_read_depth_plot_not_binned_but_stranded(self):
        # no ymax
        with patch.object(sys, 'argv', [
            'mavis', SUBCOMMAND.OVERLAY, 'GAGE4', '--annotations', ANNOTATIONS, '--output', self.temp_output,
            '--read_depth_plot', 'axis', BAM, '1', 'none', 'True'
        ]):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)
                self.assertTrue(glob_exists(os.path.join(self.temp_output, '*GAGE4*.svg')))

    def tearDown(self):
        # remove the temp directory and outputs
        shutil.rmtree(self.temp_output)
