import json
import os
import shutil
import sys
import tempfile
from unittest.mock import patch

import pytest
from mavis_config.constants import SUBCOMMAND

from mavis.main import main

from ..util import get_data, glob_exists

ANNOTATIONS = get_data('annotations_subsample.json')
BAM = get_data('mock_reads_for_events.sorted.bam')


@pytest.fixture
def output_dir():
    temp_output = tempfile.mkdtemp()
    yield temp_output
    shutil.rmtree(temp_output)


@pytest.fixture(scope='module')
def config_json():
    _, p = tempfile.mkstemp()
    print(p)
    with open(p, 'w') as fh:
        fh.write(json.dumps({'reference.annotations': [ANNOTATIONS]}))
    yield p


class TestOverlayOptions:
    def test_basic(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_marker(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--marker',
                'm',
                '49364900',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_marker_range(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--marker',
                'm',
                '49364900',
                '49365900',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_marker_not_enough_args(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--marker',
                'm',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code != 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                '--config',
                config_json,
                'GAGE4',
                '--output',
                output_dir,
                '--marker',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code != 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_marker_not_int(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--marker',
                'm',
                'k',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code != 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_read_depth_plot(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--read_depth_plot',
                'axis',
                BAM,
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_read_depth_plot_binned(self, config_json, output_dir):
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--read_depth_plot',
                'axis',
                BAM,
                '0.5',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))

    def test_read_depth_plot_not_binned_but_stranded(self, config_json, output_dir):
        # no ymax
        with patch.object(
            sys,
            'argv',
            [
                'mavis',
                SUBCOMMAND.OVERLAY,
                'GAGE4',
                '--config',
                config_json,
                '--output',
                output_dir,
                '--read_depth_plot',
                'axis',
                BAM,
                '1',
                'none',
                'True',
            ],
        ):
            try:
                returncode = main()
            except SystemExit as err:
                assert err.code == 0
            else:
                assert returncode is None
                assert glob_exists(os.path.join(output_dir, '*GAGE4*.svg'))
