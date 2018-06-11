import os
import subprocess
import sys
import unittest
from unittest.mock import patch


from mavis.constants import SUBCOMMAND
from mavis.main import main


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
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SETUP, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_config(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CONFIG, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_cluster(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CLUSTER, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_validate(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.VALIDATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_annotate(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.ANNOTATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_pairing(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.PAIR, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_summary(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SUMMARY, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_convert(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CONVERT, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_overlay(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.OVERLAY, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertEqual(0, err.code)
            else:
                self.assertEqual(0, returncode)

    def test_bad_option(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SETUP, '--blargh']):
            try:
                returncode = main()
            except SystemExit as err:
                self.assertNotEqual(0, err.code)
            else:
                self.assertNotEqual(0, returncode)
