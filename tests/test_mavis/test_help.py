import sys
from unittest.mock import patch

from mavis_config.constants import SUBCOMMAND

from mavis.main import main


class TestHelpMenu:
    def test_main(self):
        with patch.object(sys, 'argv', ['mavis', '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_pipeline(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SETUP, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_cluster(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CLUSTER, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_validate(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.VALIDATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_annotate(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.ANNOTATE, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_pairing(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.PAIR, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_summary(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SUMMARY, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_convert(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.CONVERT, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_overlay(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.OVERLAY, '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0

    def test_bad_option(self):
        with patch.object(sys, 'argv', ['mavis', SUBCOMMAND.SETUP, '--blargh']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code != 0
            else:
                assert returncode != 0

    def test_ref_alt_count(self):
        with patch.object(sys, 'argv', ['calculate_ref_alt_counts', '-h']):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode == 0
