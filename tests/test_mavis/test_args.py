import json
import sys
import tempfile
from unittest.mock import Mock, patch

import pytest
from mavis.cluster import main as cluster_main
from mavis.main import main as mavis_main
from mavis.validate import main as validate_main

from ..util import get_data


@pytest.fixture
def output_dir():
    temp_output = tempfile.mkdtemp()
    yield temp_output


@pytest.fixture
def configpath(tmp_path):
    p = tmp_path / "config.json"
    return p


def expect_error(testcase, func, catchtype=None):
    try:
        func()
    except (SystemExit, Exception) as err:
        if catchtype is None or isinstance(err, catchtype):
            return err
        raise AssertionError('Did not throw the expected error', catchtype)


class TestCluster:
    def test_trans_multiple_annotations_no_masking(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'transcriptome',
                            'assign': [get_data('mock_sv_events.tsv')],
                        }
                    },
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'cluster',
            '--library',
            'translib',
            '--inputs',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(cluster_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_trans_multiple_annotations_with_masking(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'transcriptome',
                            'assign': [get_data('mock_sv_events.tsv')],
                        }
                    },
                    'cluster.uninformative_filter': True,
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'reference.masking': [get_data('mock_masking.tab')],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'cluster',
            '--library',
            'translib',
            '--inputs',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(cluster_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_error_missing_annotations_translib_uninform(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'transcriptome',
                            'assign': [get_data('mock_sv_events.tsv')],
                        }
                    },
                    'cluster.uninformative_filter': True,
                    'output_dir': output_dir,
                }
            )
        )
        args = ['mavis', 'cluster', '--library', 'translib', '--output', output_dir]
        with patch.object(cluster_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)


class TestValidate:
    def test_error_missing_annotations_translib(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'transcriptome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'cluster.uninformative_filter': True,
                    'reference.reference_genome': [get_data('mock_reference_genome.fa')],
                    'reference.aligner_reference': [get_data('mock_reference_genome.fa')],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)

    def test_ok_multi_ref_genome(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'genome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'cluster.uninformative_filter': True,
                    'reference.reference_genome': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'reference.aligner_reference': [get_data('mock_reference_genome.fa')],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_error_multi_aligner_ref(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'genome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'cluster.uninformative_filter': True,
                    'reference.reference_genome': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'reference.aligner_reference': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)

    def test_error_missing_aligner_ref(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'genome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'cluster.uninformative_filter': True,
                    'reference.reference_genome': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)

    def test_error_missing_reference_genome(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'genome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'cluster.uninformative_filter': True,
                    'reference.aligner_reference': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)

    def test_error_bad_aligner_ref(self, configpath, output_dir):
        configpath.write_text(
            json.dumps(
                {
                    'libraries': {
                        'translib': {
                            'disease_status': 'diseased',
                            'protocol': 'genome',
                            'assign': [get_data('mock_sv_events.tsv')],
                            'bam_file': get_data('mock_trans_reads_for_events.sorted.bam'),
                            'read_length': 125,
                            'median_fragment_size': 200,
                            'stdev_fragment_size': 50,
                        }
                    },
                    'reference.annotations': [
                        get_data('example_genes.json'),
                        get_data('mock_annotations.json'),
                    ],
                    'cluster.uninformative_filter': True,
                    'reference.reference_genome': [
                        get_data('mock_reference_genome.fa'),
                        get_data('example_genes.fa'),
                    ],
                    'reference.aligner_reference': [
                        'fake_path',
                    ],
                    'output_dir': output_dir,
                }
            )
        )
        args = [
            'mavis',
            'validate',
            '--library',
            'translib',
            '--input',
            get_data('mock_sv_events.tsv'),
            '--output',
            output_dir,
            '--config',
            str(configpath),
        ]
        with patch.object(validate_main, 'main', Mock()):
            with patch.object(sys, 'argv', args):
                expect_error(self, mavis_main)
