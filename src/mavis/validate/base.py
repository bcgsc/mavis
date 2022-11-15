import itertools
from abc import abstractmethod
from typing import Dict, List, Optional, Set, Tuple

import pysam
from mavis_config import DEFAULTS

from ..bam import cigar as _cigar
from ..bam import read as _read
from ..bam.cache import BamCache
from ..breakpoint import Breakpoint, BreakpointPair, classify_breakpoint_pair
from ..constants import COLUMNS, ORIENT, PYSAM_READ_FLAGS, STRAND, SVTYPE, reverse_complement
from ..error import NotSpecifiedError
from ..interval import Interval
from ..util import logger
from .assemble import Contig, assemble


class Evidence(BreakpointPair):
    assembly_max_kmer_size: int
    bam_cache: BamCache
    classification: Optional[str]
    compatible_flanking_pairs: Set
    compatible_window1: Optional[Interval]
    compatible_window2: Optional[Interval]
    config: Dict
    contigs: List[Contig]
    counts: List[int]
    flanking_pairs: Set
    half_mapped: Tuple[Set, Set]
    median_fragment_size: int
    read_length: int
    reference_genome: Dict
    spanning_reads: Set
    split_reads: Tuple[Set, Set]
    stdev_fragment_size: int
    strand_determining_read: int
    # abstract properties
    inner_window1: Interval
    inner_window2: Interval
    outer_window1: Interval
    outer_window2: Interval

    @property
    def min_expected_fragment_size(self):
        # cannot be negative
        return int(
            round(
                max(
                    [
                        self.median_fragment_size
                        - self.stdev_fragment_size * self.config['validate.stdev_count_abnormal'],
                        0,
                    ]
                ),
                0,
            )
        )

    @property
    def max_expected_fragment_size(self):
        return int(
            round(
                self.median_fragment_size
                + self.stdev_fragment_size * self.config['validate.stdev_count_abnormal'],
                0,
            )
        )

    @property
    @abstractmethod
    def min_mapping_quality(self):
        pass

    def __init__(
        self,
        break1,
        break2,
        bam_cache,
        reference_genome,
        read_length,
        stdev_fragment_size,
        median_fragment_size,
        stranded=False,
        opposing_strands=None,
        untemplated_seq=None,
        classification=None,
        config=DEFAULTS,
        assembly_max_kmer_size=None,
        strand_determining_read=2,
        **kwargs,
    ):
        """
        Args:
            breakpoint_pair (BreakpointPair): the breakpoint pair to collect evidence for
            bam_cache (BamCache): the bam cache (and assc file) to collect evidence from
            reference_genome (Dict[str,Bio.SeqRecord]):
              dict of reference sequence by template/chr name
            data (dict): a dictionary of data to associate with the evidence object
            classification (SVTYPE): the event type
            protocol (PROTOCOL): genome or transcriptome
        """
        # initialize the breakpoint pair
        self.bam_cache = bam_cache
        self.stranded = stranded and bam_cache.stranded
        self.config = dict(**DEFAULTS)
        self.config.update(config)
        BreakpointPair.__init__(
            self,
            break1,
            break2,
            stranded=stranded,
            opposing_strands=opposing_strands,
            untemplated_seq=untemplated_seq,
            **kwargs,
        )
        # check that the breakpoints are within the reference length
        if reference_genome:
            if self.break1.start < 1 or self.break1.end > len(
                reference_genome[self.break1.chr].seq
            ):
                raise ValueError(
                    'Breakpoint {}-{} is outside the range of the reference sequence {} (1-{})'.format(
                        self.break1.start,
                        self.break1.end,
                        self.break1.chr,
                        len(reference_genome[self.break1.chr].seq),
                    )
                )
            if self.break2.start < 1 or self.break2.end > len(
                reference_genome[self.break2.chr].seq
            ):
                raise ValueError(
                    'Breakpoint {}-{} is outside the range of the reference sequence {} (1-{})'.format(
                        self.break2.start,
                        self.break2.end,
                        self.break2.chr,
                        len(reference_genome[self.break2.chr].seq),
                    )
                )
        self.assembly_max_kmer_size = (
            assembly_max_kmer_size if assembly_max_kmer_size is not None else int(read_length * 0.7)
        )
        self.bam_cache = bam_cache
        self.classification = classification
        self.compatible_window1 = None
        self.compatible_window2 = None
        self.median_fragment_size = median_fragment_size
        self.read_length = read_length
        self.reference_genome = reference_genome
        self.stdev_fragment_size = stdev_fragment_size
        self.strand_determining_read = strand_determining_read

        if self.classification is not None and self.classification not in classify_breakpoint_pair(
            self
        ):
            raise AttributeError(
                'breakpoint pair improper classification',
                classify_breakpoint_pair(self),
                self.classification,
            )

        if self.break1.orient == ORIENT.NS or self.break2.orient == ORIENT.NS:
            raise NotSpecifiedError(
                'input breakpoint pair must specify strand and orientation. Cannot be \'not specified'
                '\' for evidence gathering'
            )

        self.split_reads = (set(), set())
        self.flanking_pairs = set()
        self.compatible_flanking_pairs = set()
        self.spanning_reads = set()
        # for each breakpoint stores the number of reads that were read from the associated
        # bamfile for the window surrounding the breakpoint
        self.counts = [0, 0]  # has to be a list to assign
        self.contigs = []

        self.half_mapped = (set(), set())

        try:
            self.compute_fragment_size(None, None)
        except NotImplementedError:
            raise NotImplementedError('abstract class cannot be initialized')
        except BaseException:
            pass

    @staticmethod
    def distance(start: int, end: int):
        return Interval(abs(end - start))

    @staticmethod
    def traverse(start: int, distance: int, direction: str) -> Interval:
        if direction == ORIENT.LEFT:
            return Interval(start - distance)
        return Interval(start + distance)

    def collect_from_outer_window(self):
        """
        determines if evidence should be collected from the outer window (looking for flanking evidence)
        or should be limited to the inner window (split/spanning/contig only)

        Returns:
            bool: True or False
        """
        if self.interchromosomal:
            return True
        elif len(self.break1 | self.break2) >= self.outer_window_min_event_size:
            return True
        return False

    def standardize_read(self, read):
        # recomputing to standardize b/c split reads can be used to call breakpoints exactly
        read.set_tag(PYSAM_READ_FLAGS.RECOMPUTED_CIGAR, 1, value_type='i')
        # recalculate the read cigar string to ensure M is replaced with = or X
        cigar = _cigar.recompute_cigar_mismatch(
            read, self.reference_genome[self.bam_cache.get_read_reference_name(read)].seq
        )
        prefix = 0
        try:
            cigar, prefix = _cigar.extend_softclipping(
                cigar, self.config['validate.min_anchor_exact']
            )
        except AttributeError:
            pass
        read.cigar = _cigar.join(cigar)
        read.cigar = _cigar.merge_internal_events(
            read.cigar,
            inner_anchor=self.config['validate.contig_aln_merge_inner_anchor'],
            outer_anchor=self.config['validate.contig_aln_merge_outer_anchor'],
        )
        read.reference_start = read.reference_start + prefix

        # makes sure all indels are called as far 'right' as possible
        read.cigar = _cigar.hgvs_standardize_cigar(
            read, self.reference_genome[self.bam_cache.get_read_reference_name(read)].seq
        )
        return read

    def putative_event_types(self):
        """
        Returns:
            List[mavis.constants.SVTYPE]: list of the possible classifications
        """
        if self.classification:
            return {self.classification}
        return classify_breakpoint_pair(self)

    @property
    def compatible_type(self):
        if SVTYPE.INS in self.putative_event_types():
            return SVTYPE.DUP
        elif SVTYPE.DUP in self.putative_event_types():
            return SVTYPE.INS
        return None

    def compute_fragment_size(self, read: pysam.AlignedSegment, mate: pysam.AlignedSegment):
        """
        Returns:
            Interval: interval representing the range of possible fragment sizes for this read pair
        """
        raise NotImplementedError('abstract method must be overridden')

    def supporting_reads(self):
        """
        convenience method to return all flanking, split and spanning reads associated with an evidence object
        """
        result = set()
        for read_set in self.flanking_pairs:
            result.update(read_set)
        for read_set in self.split_reads:
            result.update(read_set)
        result.update(self.spanning_reads)
        return result

    def decide_sequenced_strand(self, reads: Set[pysam.AlignedSegment]):
        """
        given a set of reads, determines the sequenced strand (if possible) and then returns the majority
        strand found

        Args:
            reads (Set[pysam.AlignedSegment)]: set of reads

        Returns:
            STRAND: the sequenced strand

        Raises:
            ValueError: input was an empty set or the ratio was not sufficient to decide on a strand
        """
        if not reads:
            raise ValueError('cannot determine the strand of a set of reads if the set is empty')

        strand_calls = {STRAND.POS: 0, STRAND.NEG: 0}
        for read in reads:
            try:
                strand = _read.sequenced_strand(read, self.strand_determining_read)
                strand_calls[strand] = strand_calls.get(strand, 0) + 1
            except ValueError:
                pass
        if sum(strand_calls.values()) == 0:
            raise ValueError('Could not determine strand. Insufficient mapped reads')
        if strand_calls[STRAND.POS] == 0:
            return STRAND.NEG
        elif strand_calls[STRAND.NEG] == 0:
            return STRAND.POS
        else:
            ratio = strand_calls[STRAND.POS] / (strand_calls[STRAND.NEG] + strand_calls[STRAND.POS])
            neg_ratio = 1 - ratio
            if ratio >= self.config['validate.assembly_strand_concordance']:
                return STRAND.POS
            elif neg_ratio >= self.config['validate.assembly_strand_concordance']:
                return STRAND.NEG
            raise ValueError(
                'Could not determine the strand. Equivocal POS/(NEG + POS) ratio',
                ratio,
                strand_calls,
            )

    def assemble_contig(self):
        """
        uses the split reads and the partners of the half mapped reads to create a contig
        representing the sequence across the breakpoints

        if it is not strand specific then sequences are sorted alphanumerically and only the
        first of a pair is kept (paired by sequence)
        """
        # gather reads for the putative assembly
        assembly_sequences = {}
        # add split reads
        for read in list(itertools.chain.from_iterable(self.split_reads)) + list(
            self.spanning_reads
        ):
            # ignore targeted realignments
            if read.has_tag(PYSAM_READ_FLAGS.TARGETED_ALIGNMENT) and read.get_tag(
                PYSAM_READ_FLAGS.TARGETED_ALIGNMENT
            ):
                continue
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

        # add half-mapped reads
        for read in itertools.chain.from_iterable(self.half_mapped):
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

        # add flanking reads
        for read, mate in self.flanking_pairs:
            assembly_sequences.setdefault(read.query_sequence, set()).add(read)
            rqs_comp = reverse_complement(read.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(read)

            assembly_sequences.setdefault(mate.query_sequence, set()).add(mate)
            rqs_comp = reverse_complement(mate.query_sequence)
            assembly_sequences.setdefault(rqs_comp, set()).add(mate)

        logger.info(f'assembly size of {len(assembly_sequences) // 2} sequences')

        kmer_size = self.read_length * self.config['validate.assembly_kmer_size']
        remap_min_overlap = max(
            self.read_length - self.config['validate.assembly_min_exact_match_to_remap'], kmer_size
        )

        contigs = assemble(
            assembly_sequences,
            kmer_size,
            min_edge_trim_weight=self.config['validate.assembly_min_edge_trim_weight'],
            assembly_max_paths=self.config['validate.assembly_max_paths'],
            min_contig_length=self.read_length,
            remap_min_overlap=remap_min_overlap,
            remap_min_exact_match=self.config['validate.assembly_min_exact_match_to_remap'],
            assembly_min_uniq=self.config['validate.assembly_min_uniq'],
            min_complexity=self.config['validate.min_call_complexity'],
        )

        # add the input reads
        # drop any contigs without reads from both breakpoints
        filtered_contigs = []
        for ctg in contigs:
            for read_seq in ctg.remapped_sequences:
                ctg.input_reads.update(assembly_sequences[read_seq.query_sequence])
            break1_reads = {
                r.query_sequence
                for r in self.split_reads[0] | self.half_mapped[0] | self.spanning_reads
            }
            break2_reads = {
                r.query_sequence
                for r in self.split_reads[1] | self.half_mapped[1] | self.spanning_reads
            }
            for read, mate in self.flanking_pairs | self.compatible_flanking_pairs:
                break1_reads.add(read.query_sequence)
                break2_reads.add(mate.query_sequence)

            ctg_reads = {r.query_sequence for r in ctg.input_reads}
            ctg_reads.update({reverse_complement(r) for r in ctg_reads})
            if (ctg_reads & break1_reads and ctg_reads & break2_reads) or (
                not self.interchromosomal and len(self.break1 | self.break2) < self.read_length
            ):
                filtered_contigs.append(ctg)
        logger.info(
            f'filtered contigs from {len(contigs)} to {len(filtered_contigs)} based on remapped reads from both breakpoints'
        )
        contigs = filtered_contigs

        # now determine the strand from the remapped reads if possible
        if self.stranded and self.bam_cache.stranded:  # strand specific
            for contig in contigs:
                build_strand = {STRAND.POS: 0, STRAND.NEG: 0}  # if neg will have to flip
                for read_seq in contig.remapped_sequences:
                    for read in assembly_sequences[read_seq.query_sequence]:
                        if read.is_unmapped:
                            continue
                        flip = False
                        if read.query_sequence != read_seq.query_sequence:
                            flip = not flip
                        try:
                            seq_strand = _read.sequenced_strand(read, self.strand_determining_read)
                            if seq_strand == STRAND.NEG:
                                flip = not flip
                            build_strand[STRAND.NEG if flip else STRAND.POS] += 1
                        except ValueError:
                            pass
                if sum(build_strand.values()) == 0:
                    continue
                elif build_strand[STRAND.POS] == 0:
                    flipped_build = True
                elif build_strand[STRAND.NEG] == 0:
                    flipped_build = False
                else:
                    ratio = build_strand[STRAND.POS] / (
                        build_strand[STRAND.NEG] + build_strand[STRAND.POS]
                    )
                    neg_ratio = 1 - ratio
                    if ratio >= self.config['validate.assembly_strand_concordance']:
                        flipped_build = False
                    elif neg_ratio >= self.config['validate.assembly_strand_concordance']:
                        flipped_build = True
                    else:
                        continue
                if flipped_build:
                    contig.seq = reverse_complement(contig.seq)
                contig.strand_specific = True

        filtered_contigs = {}
        # sort so that the function is deterministic
        for contig in sorted(contigs, key=lambda x: (x.remap_score() * -1, x.seq)):
            # filter on evidence level
            if (
                contig.remap_score() < self.config['validate.assembly_min_remapped_seq']
                or contig.remap_coverage() < self.config['validate.assembly_min_remap_coverage']
            ):
                continue
            if self.stranded and self.bam_cache.stranded:
                filtered_contigs.setdefault(contig.seq, contig)
            else:
                rseq = reverse_complement(contig.seq)
                if contig.seq not in filtered_contigs and rseq not in filtered_contigs:
                    filtered_contigs[contig.seq] = contig
        self.contigs = sorted(
            list(filtered_contigs.values()), key=lambda x: (x.remap_score() * -1, x.seq)
        )

    def copy(self):
        raise NotImplementedError('not appropriate for copy of evidence')

    def flatten(self):
        row = BreakpointPair.flatten(self)
        row.update(
            {
                COLUMNS.raw_flanking_pairs: len(self.flanking_pairs),
                COLUMNS.raw_spanning_reads: len(self.spanning_reads),
                COLUMNS.raw_break1_split_reads: len(self.split_reads[0]),
                COLUMNS.raw_break2_split_reads: len(self.split_reads[1]),
                COLUMNS.raw_break1_half_mapped_reads: len(self.half_mapped[0]),
                COLUMNS.raw_break2_half_mapped_reads: len(self.half_mapped[1]),
                COLUMNS.protocol: self.protocol,
                COLUMNS.event_type: ';'.join(sorted(self.putative_event_types())),
                COLUMNS.contigs_assembled: len(self.contigs),
                COLUMNS.break1_ewindow: '{}-{}'.format(*self.outer_window1),
                COLUMNS.break2_ewindow: '{}-{}'.format(*self.outer_window2),
                COLUMNS.break1_ewindow_count: self.counts[0],
                COLUMNS.break2_ewindow_count: self.counts[1],
                COLUMNS.contigs_assembled: len(self.contigs),
            }
        )
        return row

    def get_bed_repesentation(self):
        bed = []
        name = self.data.get(COLUMNS.cluster_id, None)
        bed.append((self.break1.chr, self.outer_window1[0] - 1, self.outer_window1[1], name))
        bed.append((self.break1.chr, self.inner_window1[0] - 1, self.inner_window1[1], name))
        bed.append((self.break2.chr, self.outer_window2[0] - 1, self.outer_window2[1], name))
        bed.append((self.break2.chr, self.inner_window2[0] - 1, self.inner_window2[1], name))
        return bed

    def generate_window(self, breakpoint: Breakpoint) -> Interval:
        """
        given some input breakpoint uses the current evidence setting to determine an
        appropriate window/range of where one should search for supporting reads

        Args:
            breakpoint (Breakpoint): the breakpoint we are generating the evidence window for
            read_length (int): the read length
            call_error (int):
                adds a buffer to the calculations if confidence in the breakpoint calls is low can increase this
        Returns:
            Interval: the range where reads should be read from the bam looking for evidence for this event
        """
        call_error = self.config['validate.call_error']
        start = breakpoint.start - self.max_expected_fragment_size - call_error + 1
        end = breakpoint.end + self.max_expected_fragment_size + call_error - 1

        if breakpoint.orient == ORIENT.LEFT:
            end = breakpoint.end + call_error + self.read_length - 1
        elif breakpoint.orient == ORIENT.RIGHT:
            start = breakpoint.start - call_error - self.read_length + 1
        return Interval(max([1, start]), max([end, 1]))
