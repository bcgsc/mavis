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
def blat_output_dir():
    temp_output = tempfile.mkdtemp()

    with open(package_relative_file('tests/mini-tutorial.config.json'), 'r') as fh:
        config = json.load(fh)
    config['output_dir'] = os.path.join(temp_output, 'output_dir')
    config['validate.aligner'] = 'blat'
    with open(os.path.join(temp_output, 'mini-tutorial.config.json'), 'w') as fh:
        fh.write(json.dumps(config))
    yield temp_output
    shutil.rmtree(temp_output)


@pytest.fixture
def bwa_output_dir():
    temp_output = tempfile.mkdtemp()

    with open(package_relative_file('tests/mini-tutorial.config.json'), 'r') as fh:
        config = json.load(fh)
    config['output_dir'] = os.path.join(temp_output, 'output_dir')
    config['validate.aligner'] = 'bwa mem'
    config['reference.aligner_reference'] = config['reference.reference_genome']
    with open(os.path.join(temp_output, 'mini-tutorial.config.json'), 'w') as fh:
        fh.write(json.dumps(config))
    yield temp_output
    shutil.rmtree(temp_output)


@pytest.fixture
def annotate_only_output_dir():
    temp_output = tempfile.mkdtemp()

    with open(package_relative_file('tests/mini-tutorial.annotate_only.config.json'), 'r') as fh:
        config = json.load(fh)
    config['output_dir'] = os.path.join(temp_output, 'output_dir')
    with open(os.path.join(temp_output, 'mini-tutorial.config.json'), 'w') as fh:
        fh.write(json.dumps(config))
    yield temp_output
    shutil.rmtree(temp_output)


@pytest.fixture
def output_dir(request):
    return request.getfixturevalue(request.param)


@long_running_test
@pytest.mark.parametrize('output_dir', ['blat_output_dir', 'bwa_output_dir'], indirect=True)
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

        except SystemExit as err:
            if err.code != 0:
                raise err

    for expected_file in [
        os.path.join('summary', 'MAVIS.COMPLETE'),
        os.path.join('pairing', 'MAVIS.COMPLETE'),
        os.path.join('mock-A47933', 'cluster', 'MAVIS.COMPLETE'),
        os.path.join('mock-A47933', 'annotate', 'batch-*', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'cluster', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'annotate', 'batch-*', 'MAVIS.COMPLETE'),
        os.path.join('mock-A47933', 'validate', 'batch-*', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'validate', 'batch-*', 'MAVIS.COMPLETE'),
    ]:
        if not glob_exists(os.path.join(output_dir, 'output_dir', expected_file)):
            raise AssertionError(f'{expected_file} does not exist')


@long_running_test
@pytest.mark.parametrize('output_dir', ['annotate_only_output_dir'], indirect=True)
def test_no_validate_worflow(output_dir):
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

        except SystemExit as err:
            if err.code != 0:
                raise err

    for expected_file in [
        os.path.join('summary', 'MAVIS.COMPLETE'),
        os.path.join('pairing', 'MAVIS.COMPLETE'),
        os.path.join('mock-A47933', 'cluster', 'MAVIS.COMPLETE'),
        os.path.join('mock-A47933', 'annotate', 'batch-*', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'cluster', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'annotate', 'batch-*', 'MAVIS.COMPLETE'),
    ]:
        if not glob_exists(os.path.join(output_dir, 'output_dir', expected_file)):
            raise AssertionError(f'{expected_file} does not exist')

    for unexpected_file in [
        os.path.join('mock-A47933', 'validate', 'batch-*', 'MAVIS.COMPLETE'),
        os.path.join('mock-A36971', 'validate', 'batch-*', 'MAVIS.COMPLETE'),
    ]:
        if glob_exists(os.path.join(output_dir, 'output_dir', unexpected_file)):
            raise AssertionError(f'{unexpected_file} exists')
