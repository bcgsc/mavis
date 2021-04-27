import json
import os
import shutil
import sys
import tempfile
from unittest.mock import patch

import pytest

from snakemake import main as snakemake_main

from ..util import glob_exists, long_running_test, package_relative_file


@pytest.fixture
def output_dir():
    temp_output = tempfile.mkdtemp()

    os.makedirs(os.path.join(temp_output, 'mavis/schemas'))

    with open(package_relative_file('tests/mini-tutorial.config.json'), 'r') as fh:
        config = json.load(fh)
    config['output_dir'] = os.path.join(temp_output, 'output_dir')
    with open(os.path.join(temp_output, 'mini-tutorial.config.json'), 'w') as fh:
        fh.write(json.dumps(config))
    yield temp_output
    shutil.rmtree(temp_output)


@long_running_test
def test_workflow(output_dir):
    argv = [
        'snakemake',
        '-s',
        package_relative_file('Snakefile'),
        '-j',
        '1',
        '--configfile',
        os.path.join(output_dir, 'mini-tutorial.config.json'),
        '-d',
        package_relative_file(),
    ]
    with patch.object(sys, 'argv', argv):
        try:
            snakemake_main()
            assert glob_exists(os.path.join(output_dir, 'summary', 'MAVIS.COMPLETE'))
            assert glob_exists(os.path.join(output_dir, 'pairing', 'MAVIS.COMPLETE'))
            assert glob_exists(os.path.join(output_dir, 'mock-A47933', 'cluster', 'MAVIS.COMPLETE'))
            assert glob_exists(
                os.path.join(output_dir, 'mock-A47933', 'validate', '*', 'MAVIS.COMPLETE')
            )
            assert glob_exists(
                os.path.join(output_dir, 'mock-A47933', 'annotate', '*', 'MAVIS.COMPLETE')
            )
            assert glob_exists(os.path.join(output_dir, 'mock-A36971', 'cluster', 'MAVIS.COMPLETE'))
            assert glob_exists(
                os.path.join(output_dir, 'mock-A36971', 'validate', '*', 'MAVIS.COMPLETE')
            )
            assert glob_exists(
                os.path.join(output_dir, 'mock-A36971', 'annotate', '*', 'MAVIS.COMPLETE')
            )
        except SystemExit as err:
            if err.code != 0:
                raise err
