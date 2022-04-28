import pytest
from mavis.annotate.file_io import load_reference_genome
from mavis.bam import cigar as _cigar
from mavis.bam.cache import BamCache
from mavis.bam.read import SamRead
from mavis.breakpoint import Breakpoint
from mavis.constants import NA_MAPPING_QUALITY, ORIENT, PYSAM_READ_FLAGS
from mavis.validate.base import Evidence
from mavis.validate.evidence import GenomeEvidence
from mavis.validate.gather import collect_split_read, collect_flanking_pair, load_evidence
from mavis_config import DEFAULTS

from ...util import get_data, long_running_test
from ..mock import MockLongString, MockObject, MockRead, mock_read_pair

REFERENCE_GENOME = None


def setUpModule():
    global REFERENCE_GENOME
    REFERENCE_GENOME = load_reference_genome(get_data('mock_reference_genome.fa'))
    if (
        'CTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGT'
        != REFERENCE_GENOME['fake'].seq[0:50].upper()
    ):
        raise AssertionError('fake genome file does not have the expected contents')
    global BAM_CACHE
    BAM_CACHE = BamCache(get_data('mini_mock_reads_for_events.sorted.bam'))
    global FULL_BAM_CACHE
    FULL_BAM_CACHE = BamCache(get_data('mock_reads_for_events.sorted.bam'))
    global READS
    READS = {}
    for read in BAM_CACHE.fetch('reference3', 1, 8000):
        if read.qname not in READS:
            READS[read.qname] = [None, None]
        if read.is_supplementary:
            continue
        if read.is_read1:
            READS[read.qname][0] = read
        else:
            READS[read.qname][1] = read
    # add a check to determine if it is the expected bam file


@long_running_test
class TestFullEvidenceGathering:
    # need to make the assertions more specific by checking the actual names of the reads found in each bin
    # rather than just the counts.
    def genome_evidence(self, break1, break2, opposing_strands):
        ge = GenomeEvidence(
            break1,
            break2,
            FULL_BAM_CACHE,
            REFERENCE_GENOME,
            opposing_strands=opposing_strands,
            read_length=125,
            stdev_fragment_size=100,
            median_fragment_size=380,
            config={
                'validate.stdev_count_abnormal': 3,
                'validate.min_flanking_pairs_resolution': 3,
                'validate.max_sc_preceeding_anchor': 3,
                'validate.outer_window_min_event_size': 0,
                'validate.min_mapping_quality': 20,
            },
        )
        print(ge.min_expected_fragment_size, ge.max_expected_fragment_size)
        print(ge.break1.chr, ge.outer_window1)
        print(ge.break1.chr, ge.inner_window1)
        print(ge.break2.chr, ge.outer_window2)
        print(ge.break2.chr, ge.inner_window2)
        return ge

    def print_evidence(self, ev):
        print('evidence for', ev)
        print('flanking pairs')
        for pair, mate in ev.flanking_pairs:
            print(
                pair.query_name,
                pair.reference_name,
                ':',
                pair.reference_start,
                _cigar.convert_cigar_to_string(pair.cigar),
            )
            print(
                mate.query_name,
                mate.reference_name,
                ':',
                mate.reference_start,
                _cigar.convert_cigar_to_string(mate.cigar),
            )

        print('first breakpoint split reads')
        for read in ev.split_reads[0]:
            print(
                read.query_name,
                read.reference_name,
                ':',
                read.reference_start,
                _cigar.convert_cigar_to_string(read.cigar),
            )

        print('second breakpoint split reads')
        for read in ev.split_reads[1]:
            print(
                read.query_name,
                read.reference_name,
                ':',
                read.reference_start,
                _cigar.convert_cigar_to_string(read.cigar),
            )
        print()

    def count_original_reads(self, reads):
        count = 0
        for read in sorted(reads, key=lambda x: x.query_name):
            if not read.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                count += 1
            elif not read.get_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT):
                count += 1
        return count

    def test_load_evidence_translocation(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 520, orient=ORIENT.RIGHT),
            Breakpoint('reference19', 964, orient=ORIENT.LEFT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 14
        assert self.count_original_reads(ev1.split_reads[1]) == 20
        assert len(ev1.flanking_pairs) == 21

        # second example
        ev1 = self.genome_evidence(
            Breakpoint('reference2', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference4', 2000, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 21
        # one of the reads that appears to look good in the bam is too low quality % match
        assert self.count_original_reads(ev1.split_reads[1]) == 40
        assert len(ev1.flanking_pairs) == 57

    def test_load_evidence_inversion(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
            Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
            opposing_strands=True,
        )

        load_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 54
        assert self.count_original_reads(ev1.split_reads[1]) == 20
        assert len(ev1.flanking_pairs) == 104

        # second example
        ev1 = self.genome_evidence(
            Breakpoint('reference7', 15000, orient=ORIENT.RIGHT),
            Breakpoint('reference7', 19000, orient=ORIENT.RIGHT),
            opposing_strands=True,
        )
        load_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[1]) == 15
        assert self.count_original_reads(ev1.split_reads[0]) == 27
        assert len(ev1.flanking_pairs) == 52

    def test_load_evidence_duplication(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference7', 5000, orient=ORIENT.RIGHT),
            Breakpoint('reference7', 11000, orient=ORIENT.LEFT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 35
        assert self.count_original_reads(ev1.split_reads[1]) == 11
        assert len(ev1.flanking_pairs) == 64

    def test_load_evidence_deletion1(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference20', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference20', 6000, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert len(ev1.flanking_pairs) == 49
        assert self.count_original_reads(ev1.split_reads[0]) == 22
        assert self.count_original_reads(ev1.split_reads[1]) == 14

    def test_load_evidence_deletion2(self):
        # second example
        ev1 = self.genome_evidence(
            Breakpoint('referenceX', 2000, orient=ORIENT.LEFT),
            Breakpoint('referenceX', 6000, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 4
        assert self.count_original_reads(ev1.split_reads[1]) == 10
        assert len(ev1.flanking_pairs) == 27

    def test_load_evidence_deletion3(self):
        # third example
        ev1 = self.genome_evidence(
            Breakpoint('referenceX', 10000, orient=ORIENT.LEFT),
            Breakpoint('referenceX', 14000, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 8
        assert self.count_original_reads(ev1.split_reads[1]) == 9
        assert len(ev1.flanking_pairs) == 26

    def test_load_evidence_deletion4(self):
        # forth example
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 3609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 3818, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert self.count_original_reads(ev1.split_reads[0]) == 20
        assert self.count_original_reads(ev1.split_reads[1]) == 18
        assert len(ev1.flanking_pairs) == 40

    def test_load_evidence_small_deletion1(self):
        # first example
        ev1 = self.genome_evidence(
            Breakpoint('reference11', 6000, orient=ORIENT.LEFT),
            Breakpoint('reference11', 6003, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))

        assert self.count_original_reads(ev1.split_reads[0]) == 5
        assert self.count_original_reads(ev1.split_reads[1]) == 3
        assert len(ev1.spanning_reads) == 20
        assert len(ev1.flanking_pairs) == 6

    def test_load_evidence_small_deletion2(self):
        # second example
        ev1 = self.genome_evidence(
            Breakpoint('reference11', 10000, orient=ORIENT.LEFT),
            Breakpoint('reference11', 10030, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read, mate in ev1.flanking_pairs:
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 27
        assert self.count_original_reads(ev1.split_reads[1]) == 52
        assert len(ev1.spanning_reads) == 19
        assert len(ev1.flanking_pairs) == 7

    def test_load_evidence_small_deletion_test1(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 2001, orient=ORIENT.LEFT),
            Breakpoint('reference12', 2120, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read, mate in ev1.flanking_pairs:
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 18
        assert self.count_original_reads(ev1.split_reads[1]) == 16
        assert len(ev1.spanning_reads) == 0
        assert len(ev1.flanking_pairs) == 22

    def test_load_evidence_small_deletion_test2(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 3609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 3818, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        assert self.count_original_reads(ev1.split_reads[0]) == 20
        assert self.count_original_reads(ev1.split_reads[1]) == 18
        assert len(ev1.spanning_reads) == 0
        assert len(set(ev1.flanking_pairs)) == 40

    def test_load_evidence_small_deletion_test3(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 8609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 8927, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 27
        assert self.count_original_reads(ev1.split_reads[1]) == 5
        assert len(ev1.spanning_reads) == 0
        assert len(ev1.flanking_pairs) == 53

    def test_load_evidence_small_deletion_test4(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 12609, orient=ORIENT.LEFT),
            Breakpoint('reference10', 13123, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        assert self.count_original_reads(ev1.split_reads[0]) == 33
        assert self.count_original_reads(ev1.split_reads[1]) == 6
        assert len(ev1.spanning_reads) == 0
        assert len(ev1.flanking_pairs) == 77

    def test_load_evidence_small_deletion_test5(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 17109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 17899, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        assert self.count_original_reads(ev1.split_reads[0]) == 19
        assert self.count_original_reads(ev1.split_reads[1]) == 11
        assert len(ev1.spanning_reads) == 0
        assert len(ev1.flanking_pairs) == 48

    def test_load_evidence_small_deletion_test6(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 22109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 24330, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        assert self.count_original_reads(ev1.split_reads[0]) == 18
        assert self.count_original_reads(ev1.split_reads[1]) == 13
        assert len(ev1.flanking_pairs) == 53

    def test_load_evidence_small_deletion_test7(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 28109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 31827, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        assert self.count_original_reads(ev1.split_reads[0]) == 39
        assert self.count_original_reads(ev1.split_reads[1]) == 13
        assert len(ev1.flanking_pairs) == 49

    def test_load_evidence_small_deletion_test8(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference10', 36109, orient=ORIENT.LEFT),
            Breakpoint('reference10', 42159, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        assert self.count_original_reads(ev1.split_reads[0]) == 59
        assert self.count_original_reads(ev1.split_reads[1]) == 8
        assert len(ev1.flanking_pairs) == 59

    @pytest.mark.skip(reason='skip because too complex')
    def test_load_evidence_complex_deletion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 6001, orient=ORIENT.LEFT),
            Breakpoint('reference12', 6016, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        print(self.count_original_reads(ev1.split_reads[0]))
        print(self.count_original_reads(ev1.split_reads[1]))

        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 76
        assert self.count_original_reads(ev1.split_reads[1]) == 83
        assert len(ev1.spanning_reads) == 1
        assert len(ev1.flanking_pairs) == 2

    @pytest.mark.skip(reason='skip because high coverage')
    def test_load_evidence_small_insertion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference1', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference1', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 17
        assert self.count_original_reads(ev1.split_reads[1]) == 17
        assert len(ev1.spanning_reads) == 48
        assert len(ev1.flanking_pairs) == 4

    @pytest.mark.skip(reason='skip because too high coverage')
    def test_load_evidence_small_insertion_high_coverage(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference9', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference9', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 37
        assert self.count_original_reads(ev1.split_reads[1]) == 52
        assert len(ev1.spanning_reads) == 37
        assert len(ev1.flanking_pairs) == 9

        ev1 = self.genome_evidence(
            Breakpoint('reference16', 2000, orient=ORIENT.LEFT),
            Breakpoint('reference16', 2001, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 27
        assert self.count_original_reads(ev1.split_reads[1]) == 52
        assert len(ev1.spanning_reads) == 19
        assert len(ev1.flanking_pairs) == 9

    def test_load_evidence_small_duplication(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference12', 10000, orient=ORIENT.RIGHT),
            Breakpoint('reference12', 10021, orient=ORIENT.LEFT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 29
        assert self.count_original_reads(ev1.split_reads[1]) == 51
        assert len(ev1.spanning_reads) == 0
        assert len(ev1.flanking_pairs) == 0

        # Example 2
        ev1 = self.genome_evidence(
            Breakpoint('reference17', 1974, orient=ORIENT.RIGHT),
            Breakpoint('reference17', 2020, orient=ORIENT.LEFT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)

        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs), len(ev1.spanning_reads))
        print(len(ev1.spanning_reads))
        for read in sorted(ev1.spanning_reads, key=lambda x: x.query_name):
            print(read.query_name)

        assert self.count_original_reads(ev1.split_reads[0]) == 25
        assert self.count_original_reads(ev1.split_reads[1]) == 56
        assert len(ev1.spanning_reads) == 3
        assert len(ev1.flanking_pairs) == 0

    def test_load_evidence_low_qual_deletion(self):
        ev1 = self.genome_evidence(
            Breakpoint('reference19', 4847, 4847, orient=ORIENT.LEFT),
            Breakpoint('reference19', 5219, 5219, orient=ORIENT.RIGHT),
            opposing_strands=False,
        )
        load_evidence(ev1)
        self.print_evidence(ev1)
        print(len(ev1.spanning_reads))
        print(len(ev1.split_reads[0]), len(ev1.flanking_pairs))
        assert len(ev1.split_reads[0]) == 0
        assert len(ev1.split_reads[1]) == 0
        assert len(ev1.flanking_pairs) == 0


@pytest.fixture
def ev_gathering_setup():
    return GenomeEvidence(
        Breakpoint('reference3', 1114, orient=ORIENT.RIGHT),
        Breakpoint('reference3', 2187, orient=ORIENT.RIGHT),
        BAM_CACHE,
        REFERENCE_GENOME,
        opposing_strands=True,
        read_length=125,
        stdev_fragment_size=100,
        median_fragment_size=380,
        config={
            'validate.stdev_count_abnormal': 3,
            'validate.min_flanking_pairs_resolution': 3,
            'validate.assembly_min_edge_trim_weight': 3,
        },
    )


class TestEvidenceGathering:
    def test_collect_split_read(self, ev_gathering_setup):
        ev1_sr = MockRead(
            query_name='HISEQX1_11:3:1105:15351:25130:split',
            reference_id=1,
            cigar=[(4, 68), (7, 82)],
            reference_start=1114,
            reference_end=1154,
            query_alignment_start=110,
            query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
            query_alignment_end=150,
            flag=113,
            next_reference_id=1,
            next_reference_start=2341,
        )
        collect_split_read(ev_gathering_setup, ev1_sr, True)
        assert list(ev_gathering_setup.split_reads[0])[0] == ev1_sr

    def test_collect_split_read_failure(self, ev_gathering_setup):
        # wrong cigar string
        ev1_sr = MockRead(
            query_name='HISEQX1_11:4:1203:3062:55280:split',
            reference_id=1,
            cigar=[(7, 110), (7, 40)],
            reference_start=1114,
            reference_end=1154,
            query_alignment_start=110,
            query_sequence='CTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATG',
            query_alignment_end=150,
            flag=371,
            next_reference_id=1,
            next_reference_start=2550,
        )
        assert not collect_split_read(ev_gathering_setup, ev1_sr, True)

    def test_collect_flanking_pair(self, ev_gathering_setup):
        collect_flanking_pair(
            ev_gathering_setup,
            MockRead(
                reference_id=1,
                reference_start=2214,
                reference_end=2364,
                is_reverse=True,
                next_reference_id=1,
                next_reference_start=1120,
                mate_is_reverse=True,
            ),
            MockRead(
                reference_id=1,
                reference_start=1120,
                reference_end=2364,
                is_reverse=True,
                next_reference_id=1,
                next_reference_start=1120,
                mate_is_reverse=True,
                is_read1=False,
            ),
        )
        assert len(ev_gathering_setup.flanking_pairs) == 1

    def test_collect_flanking_pair_not_overlapping_evidence_window(self, ev_gathering_setup):
        # first read in pair does not overlap the first evidence window
        # therefore this should return False and not add to the flanking_pairs
        pair = mock_read_pair(
            MockRead(reference_id=1, reference_start=1903, reference_end=2053, is_reverse=True),
            MockRead(reference_id=1, reference_start=2052, reference_end=2053, is_reverse=True),
        )
        assert not collect_flanking_pair(ev_gathering_setup, *pair)
        assert len(ev_gathering_setup.flanking_pairs) == 0

    def test_load_evidence(self, ev_gathering_setup):
        print(ev_gathering_setup)
        load_evidence(ev_gathering_setup)
        print(ev_gathering_setup.spanning_reads)
        assert (
            len(
                [
                    r
                    for r in ev_gathering_setup.split_reads[0]
                    if not r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT)
                ]
            )
            == 2
        )
        assert len(ev_gathering_setup.flanking_pairs) == 7
        assert (
            len(
                [
                    r
                    for r in ev_gathering_setup.split_reads[1]
                    if not r.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT)
                ]
            )
            == 2
        )

    def test_assemble_split_reads(self, ev_gathering_setup):
        sr1 = MockRead(
            query_name='HISEQX1_11:3:1105:15351:25130:split',
            query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
            flag=113,
        )
        sr2 = MockRead(
            query_sequence='GTCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTT',
            flag=121,
        )
        sr3 = MockRead(
            query_sequence='TCGTGAGTGGCAGGTGCCATCGTGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTG',
            flag=113,
        )
        sr7 = MockRead(
            query_sequence='TGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATA',
            flag=113,
        )
        sr9 = MockRead(
            query_sequence='TGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGC',
            flag=113,
        )
        sr12 = MockRead(
            query_sequence='GATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAG',
            flag=113,
        )
        sr15 = MockRead(
            query_sequence='GTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACA',
            flag=113,
        )
        sr19 = MockRead(
            query_sequence='TGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCA',
            flag=113,
        )
        sr24 = MockRead(
            query_sequence='CTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTT',
            flag=113,
        )
        ev_gathering_setup.split_reads = (
            {sr1},
            {sr1, sr3, sr7, sr9, sr12, sr15, sr19, sr24},
        )  # subset needed to make a contig
        #        ev_gathering_setup.split_reads=([],[sr1,sr3,sr5,sr6,sr7,sr8,sr9,sr10,sr11,sr12,sr13,sr14,sr15,sr16,sr17,sr18,sr19,sr20,sr21,sr22,sr23,sr24]) #full set of reads produces different contig from subset.
        # full contig with more read support should be
        # CTGAGCATGAAAGCCCTGTAAACACAGAATTTGGATTCTTTCCTGTTTGGTTCCTGGTCGTGAGTGGCAGGTGCCATCATGTTTCATTCTGCCTGAGAGCAGTCTACCTAAATATATAGCTCTGCTCACAGTTTCCCTGCAATGCATAATTAAAATAGCACTATGCAGTTGCTTACACTTCAGATAATGGCTTCCTACATATTGTTGGTTATGAAATTTCAGGGTTTTCATTTCTGTATGTTAAT
        ev_gathering_setup.half_mapped = (set(), {sr2})
        ev_gathering_setup.assemble_contig()
        print(ev_gathering_setup.contigs)
        exp = 'CAACAATATGTAGGAAGCCATTATCTGAAGTGTAAGCAACTGCATAGTGCTATTTTAATTATGCATTGCAGGGAAACTGTGAGCAGAGCTATATATTTAGGTAGACTGCTCTCAGGCAGAATGAAACATGATGGCACCTGCCACTCACGACCAGGAACCAAACAGGAAAGAATC'
        assert ev_gathering_setup.contigs[0].seq == exp


class TestStandardizeRead:
    def test_bwa_mem(self):
        mock_evidence = MockObject(
            reference_genome={
                '1': MockObject(
                    seq=MockLongString(
                        'TGGGTATCAGACACACTGGGTAGCTGAGTGCTCAGAGGAAGATGCGAGGTATTCAGGGAAAGTGTCAGTGGGGTCTCCCAGTGCCTGTTTGGTCCACAGTTAGGAGA'
                        'GGCCCTGCTTGCACTTCTAATACAGTCCCGGAAAGACGGGGCCAGAACTTAGGAGGGGAGCGCTTTGCAGCAACTTTTCAAGAAAAGGGGAAAATTTAAGCACCATA'
                        'CTGTTATGTGGTCCTTGTACCCAGAGGCCCTGTTCAGCTCCAGTGATCAGCTCTCTTAGGGCACACCCTCCAAGGTGCCTAAATGCCATCCCAGGATTGGTTCCAGT'
                        'GTCTATTATCTGTTTGACTCCAAATGGCCAAACACCTGACTTCCTCTCTGGTAGCCTGGCTTTTATCTTCTAGGACATCCAGGGCCCCTCTCTTTGCCTTCCCCTCT'
                        'TTCTTCCTTCTACTGCTTAGATCAAGTCTTCAGCAGACATCATGTGACCTTGAGGATGGATGTCACATGCTGGAGGAAACAGAAGGCCGAAACCCTGATGACTTCAC'
                        'AGAGCTGCCAAAACAGTTCCTGACTGTTTATTCCGGGTCTTTAACAAAGTGATGAAAAGAAATCCTTGCAGTATGAAAACAACTTTTCTATTCCATGGAGCCAAACC'
                        'TCATTATAACAGATAACGTGACCCTCAGCGATATCCCAAGTATTTTCCTGTTCTCATCTATACTATGGCAAAGGGGCAAATACCTCTCAGTAAAGAAAGAAATAACA'
                        'ACTTCTATCTTGGGCGAGGCATTTCTTCTGTTAGAACTTTGTACACGGAATAAAATAGATCTGTTTGTGCTTATCTTTCTCCTTAGAATTATTGAATTTGAAGTCTT'
                        'TCCCAGGGTGGGGGTGGAGTGAAGCTGGGGTTTCATAAGCACATAGATAGTAGTG',
                        offset=224646450,
                    )
                )
            },
            bam_cache=MockObject(get_read_reference_name=lambda x: x.reference_name),
            config={
                'validate.contig_aln_merge_inner_anchor': 10,
                'validate.contig_aln_merge_outer_anchor': 20,
                **DEFAULTS,
            },
        )
        # SamRead(1:224646710-224646924, 183=12D19=, TCAGCTCTCT...) TCAGCTCTCTTAGGGCACACCCTCCAAGGTGCCTAAATGCCATCCCAGGATTGGTTCCAGTGTCTATTATCTGTTTGACTCCAAATGGCCAAACACCTGACTTCCTCTCTGGTAGCCTGGCTTTTATCTTCTAGGACATCCAGGGCCCCTCTCTTTGCCTTCCCCTCTTTCTTCCTTCTACTGCTTCAGCAGACATCATGTG
        # std SamRead(1:224646710-224646924, 183=12D19=, TCAGCTCTCT...) TCAGCTCTCTTAGGGCACACCCTCCAAGGTGCCTAAATGCCATCCCAGGATTGGTTCCAGTGTCTATTATCTGTTTGACTCCAAATGGCCAAACACCTGACTTCCTCTCTGGTAGCCTGGCTTTTATCTTCTAGGACATCCAGGGCCCCTCTCTTTGCCTTCCCCTCTTTCTTCCTTCTACTGCTTCAGCAGACATCATGTG
        # > BPP(Breakpoint(1:224646893L-), Breakpoint(1:224646906R-), opposing=False, seq='')
        read = SamRead(reference_name='1')
        read.query_sequence = 'TCAGCTCTCTTAGGGCACACCCTCCAAGGTGCCTAAATGCCATCCCAGGATTGGTTCCAGTGTCTATTATCTGTTTGACTCCAAATGGCCAAACACCTGACTTCCTCTCTGGTAGCCTGGCTTTTATCTTCTAGGACATCCAGGGCCCCTCTCTTTGCCTTCCCCTCTTTCTTCCTTCTACTGCTTCAGCAGACATCATGTG'
        read.reference_start = 224646710
        read.reference_id = 0
        print(_cigar.convert_string_to_cigar('183=12D19='))
        read.cigar = _cigar.join(_cigar.convert_string_to_cigar('183=12D19='))
        read.query_name = 'name'
        read.mapping_quality = NA_MAPPING_QUALITY
        std_read = Evidence.standardize_read(mock_evidence, read)
        assert std_read.cigar == _cigar.convert_string_to_cigar('186=12D16=')
        assert std_read.reference_start == read.reference_start


class MockEvidence:
    def __init__(self, ref=None):
        self.HUMAN_REFERENCE_GENOME = ref
