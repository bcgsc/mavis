import argparse
import os
import unittest
from unittest.mock import patch
import sys
from mavis.main import main as mavis_main
from mavis.cluster import main as cluster_main
from mavis.validate import main as validate_main
from mavis import util

from . import ARGUMENT_ERROR
from ..util import get_data


def expect_error(testcase, func, catchtype):
    try:
        func()
    except catchtype as err:
        return err
    else:
        raise AssertionError('Did not throw the expected error', catchtype)


class TestCluster(unittest.TestCase):

    def test_trans_multiple_annotations_no_masking(self):
        args = [
            'mavis', 'cluster',
            '--annotations', get_data('example_genes.json'), get_data('mock_annotations.json'),
            '--library', 'translib',
            '--protocol', 'transcriptome',
            '--disease_status', 'diseased',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir'
        ]
        with patch.object(cluster_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_trans_multiple_annotations_with_masking(self):
        args = [
            'mavis', 'cluster',
            '--annotations', get_data('example_genes.json'), get_data('mock_annotations.json'),
            '--library', 'translib',
            '--protocol', 'transcriptome',
            '--disease_status', 'diseased',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--masking', get_data('mock_masking.tab')
        ]
        with patch.object(cluster_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_error_missing_annotations_translib_uninform(self):
        args = [
            'mavis', 'cluster',
            '--library', 'translib',
            '--protocol', 'transcriptome',
            '--disease_status', 'diseased',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--uninformative_filter', 'True'
        ]
        with patch.object(cluster_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_ok_missing_annotations_translib_nofilter(self):
        args = [
            'mavis', 'cluster',
            '--library', 'translib',
            '--protocol', 'transcriptome',
            '--disease_status', 'diseased',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir'
        ]
        with patch.object(cluster_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                mavis_main()


class TestValidate(unittest.TestCase):

    def test_error_missing_annotations_translib(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'transcriptome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--aligner_reference', get_data('mock_reference_genome.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_ok_missing_annotations_genome(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--aligner_reference', get_data('mock_reference_genome.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_ok_multi_ref_genome(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'), get_data('example_genes.fa'),
            '--aligner_reference', get_data('mock_reference_genome.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                mavis_main()

    def test_error_multi_aligner_ref(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--aligner_reference', get_data('mock_reference_genome.fa'), get_data('example_genes.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_error_missing_aligner_ref(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_error_missing_reference_genome(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--aligner_reference', get_data('mock_reference_genome.fa'),
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_error_bad_aligner_ref(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--aligner_reference', 'bad',
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)

    def test_error_none_aligner_ref(self):
        args = [
            'mavis', 'validate',
            '--library', 'translib',
            '--protocol', 'genome',
            '--bam_file', get_data('mock_trans_reads_for_events.sorted.bam'),
            '--stdev_fragment_size', '50',
            '--median_fragment_size', '200',
            '--input', get_data('mock_sv_events.tsv'),
            '--output', 'outdir',
            '--reference_genome', get_data('mock_reference_genome.fa'),
            '--aligner_reference', 'none',
            '--read_length', '125'
        ]
        with patch.object(validate_main, 'main', util.DEVNULL):
            with patch.object(sys, 'argv', args):
                err = expect_error(self, mavis_main, SystemExit)
                self.assertEqual(ARGUMENT_ERROR, err.code)
