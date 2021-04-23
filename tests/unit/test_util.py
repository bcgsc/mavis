import os

import pytest
from mavis.constants import COLUMNS, ORIENT, STRAND
from mavis.error import NotSpecifiedError
from mavis.util import (
    ENV_VAR_PREFIX,
    cast,
    get_connected_components,
    get_env_variable,
    read_bpp_from_input_file,
)

from .mock import Mock


class TestGetConnectedComponents:
    def test_no_nodes(self):
        assert get_connected_components({}) == []

    def test_no_connections(self):
        graph = {1: {}, 2: {}, 3: {}}
        components = get_connected_components(graph)
        assert len(components) == 3

    def test_fully_connected(self):
        graph = {1: {2, 3, 1}, 2: {1, 2, 2}, 3: {3, 2}}
        components = get_connected_components(graph)
        assert len(components) == 1
        assert components == [{1, 2, 3}]

    def test_multiple_components(self):
        graph = {1: {2}, 2: {3}, 3: {4}, 6: {7, 8}}
        components = get_connected_components(graph)
        assert len(components) == 2
        assert components[0] == {1, 2, 3, 4}
        assert components[1] == {6, 7, 8}


class TestCast:
    def test_float(self):
        assert type(cast('1', float)) == type(1.0)
        assert type(cast('1', int)) != type(1.0)

    def test_boolean(self):
        assert type(cast('f', bool)) == type(False)
        assert type(cast('false', bool)) == type(False)
        assert not cast('f', bool)
        assert not cast('false', bool)
        assert not cast('0', bool)
        assert not cast('F', bool)


def mock_file_content(row):
    header = [c for c in row]
    line = [row[c] for c in header]
    lines = ['\t'.join(header), '\t'.join([str(v) for v in line])]
    return '\n'.join(lines)


class TestReadBreakpointPairsFromFile:
    def test_break1_strand_ns(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"

        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.POS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        with pytest.raises(NotSpecifiedError):
            bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=True)
            for b in bpps:
                print(b)

    def test_stranded_no_expand_error(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.POS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=True)
        assert len(bpps) == 1
        assert bpps[0].break1.strand == STRAND.POS
        assert bpps[0].break2.strand == STRAND.POS

    def test_break2_strand_ns(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: False,
                }
            )
        )

        with pytest.raises(NotSpecifiedError) as err:
            print(err)
            bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=False)

        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=True)
        assert len(bpps) == 1
        assert bpps[0].break2.strand == STRAND.POS

    def test_stranded_expand_strands_and_orient(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=False)
        assert len(bpps) == 2
        assert bpps[0].break1.strand == STRAND.POS
        assert bpps[0].break2.strand == STRAND.POS
        assert bpps[1].break1.strand == STRAND.NEG
        assert bpps[1].break2.strand == STRAND.NEG

    def test_expand_strands_and_orient(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: False,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=False)
        assert len(bpps) == 1
        assert bpps[0].break1.strand == STRAND.NS
        assert bpps[0].break2.strand == STRAND.NS

    def test_stranded_expand_strands_not_orient(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=True)
        assert len(bpps) == 2
        assert bpps[0].break1.strand == STRAND.POS
        assert bpps[0].break2.strand == STRAND.POS
        assert bpps[1].break1.strand == STRAND.NEG
        assert bpps[1].break2.strand == STRAND.NEG

    def test_expand_orient_not_strand(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.NS,
                    COLUMNS.break1_orientation: ORIENT.LEFT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: False,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=True)
        assert len(bpps) == 1
        assert bpps[0].break1.strand == STRAND.NS
        assert bpps[0].break2.strand == STRAND.NS

    def test_break1_orient_ns(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.NS,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.POS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: False,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=True)
        assert len(bpps) == 1
        assert bpps[0].break1.orient == ORIENT.LEFT

    @pytest.mark.skip(reason='TODO')
    def test_break2_orient_ns(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.NS,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.POS,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: False,
                    COLUMNS.opposing_strands: False,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=True)
        assert len(bpps) == 1
        assert bpps[0].break1.orient == ORIENT.LEFT

    @pytest.mark.skip(reason='TODO')
    def test_both_break_orient_ns(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"

    def test_base_case(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.RIGHT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NEG,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: True,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=False)
        assert len(bpps) == 1
        assert bpps[0].break1.orient == ORIENT.RIGHT
        assert bpps[0].opposing_strands == True

    def test_unstranded_with_strand_calls(self, tmp_path):
        input_file = tmp_path / "inputs.tsv"
        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.RIGHT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NEG,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: False,
                    COLUMNS.opposing_strands: True,
                }
            )
        )
        bpps = read_bpp_from_input_file(input_file, expand_strand=False, expand_orient=False)
        assert len(bpps) == 1
        assert bpps[0].break1.strand == STRAND.NS
        assert bpps[0].break2.strand == STRAND.NS

        input_file = tmp_path / "inputs2.tsv"

        input_file.write_text(
            mock_file_content(
                {
                    COLUMNS.break1_chromosome: '1',
                    COLUMNS.break1_position_start: 1,
                    COLUMNS.break1_position_end: 1,
                    COLUMNS.break1_strand: STRAND.POS,
                    COLUMNS.break1_orientation: ORIENT.RIGHT,
                    COLUMNS.break2_chromosome: '1',
                    COLUMNS.break2_position_start: 10,
                    COLUMNS.break2_position_end: 10,
                    COLUMNS.break2_strand: STRAND.NEG,
                    COLUMNS.break2_orientation: ORIENT.RIGHT,
                    COLUMNS.stranded: True,
                    COLUMNS.opposing_strands: True,
                }
            )
        )

        bpps = read_bpp_from_input_file(input_file, expand_strand=True, expand_orient=False)
        assert len(bpps) == 1
        assert bpps[0].break1.strand == STRAND.POS
        assert bpps[0].break2.strand == STRAND.NEG
