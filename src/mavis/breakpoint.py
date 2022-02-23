from copy import copy as _copy
from typing import Callable, Dict, List, Optional, Set, Tuple

from .constants import COLUMNS, DNA_ALPHABET, ORIENT, STRAND, SVTYPE, reverse_complement
from .error import InvalidRearrangement, NotSpecifiedError
from .interval import Interval
from .types import ReferenceGenome


class Breakpoint(Interval):
    """
    class for storing information about a SV breakpoint
    coordinates are given as 1-indexed
    """

    orient: str
    chr: str
    strand: str
    seq: str

    @property
    def key(self):
        return (self.chr, self.start, self.end, self.orient, self.strand)

    def __init__(
        self,
        chr: str,
        start: int,
        end: Optional[int] = None,
        orient=ORIENT.NS,
        strand=STRAND.NS,
        seq: Optional[str] = None,
    ):
        """
        Args:
            chr: the chromosome
            start: the genomic position of the breakpoint
            end: if the breakpoint is uncertain (a range) then specify the end of the range here
            orient (ORIENT): the orientation (which side is retained at the break)
            strand (STRAND): the strand
            seq: the seq

        Examples:
            >>> Breakpoint('1', 1, 2)
            >>> Breakpoint('1', 1)
            >>> Breakpoint('1', 1, 2, 'R', )
            >>> Breakpoint('1', 1, orient='R')
        """
        from .annotate.base import ReferenceName

        Interval.__init__(self, start, end)
        self.orient = ORIENT.enforce(orient)
        self.chr = ReferenceName(chr)
        self.strand = STRAND.enforce(strand)
        self.seq = seq

    def __repr__(self):
        strand = '' if self.strand == STRAND.NS else self.strand
        orient = '' if self.orient == ORIENT.NS else self.orient

        return 'Breakpoint({0}:{1}{2}{3}{4})'.format(
            self.chr,
            self.start,
            '-' + str(self.end) if self.end != self.start else '',
            orient,
            strand,
        )

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def to_dict(self):
        return {
            'chr': self.chr,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'seq': self.seq,
            'orientation': self.orient,
            'type': self.__class__.__name__,
        }


class BreakpointPair:
    break1: Breakpoint
    break2: Breakpoint
    stranded: bool
    opposing_strands: bool
    untemplated_seq: Optional[str]
    data: Dict

    def __getitem__(self, index):
        try:
            index = int(index)
        except ValueError:
            raise IndexError('index input accessor must be an integer', index)
        if index == 0:
            return self.break1
        elif index == 1:
            return self.break2
        raise IndexError('index input accessor is out of bounds: 1 or 2 only', index)

    def __hash__(self):
        return hash(
            (self.break1, self.break2, self.opposing_strands, self.stranded, self.untemplated_seq)
        )

    def __eq__(self, other):
        for attr in ['break1', 'break2', 'opposing_strands', 'stranded', 'untemplated_seq']:
            if not hasattr(other, attr):
                return False
            elif getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __lt__(self, other):
        if self.break1.key < other.break1.key:
            return True
        elif self.break1.key > other.break1.key:
            return False
        elif self.break2.key < other.break2.key:
            return True
        elif self.break2.key > other.break2.key:
            return False
        elif self.untemplated_seq == other.untemplated_seq:
            return False
        elif self.untemplated_seq is None:
            return False
        elif other.untemplated_seq is None:
            return True
        return self.untemplated_seq < other.untemplated_seq

    @property
    def library(self) -> Optional[str]:
        return self.data.get(COLUMNS.library)

    @property
    def cdna_synon(self) -> Optional[bool]:
        return self.data.get(COLUMNS.cdna_synon)

    @property
    def contig_remapped_reads(self) -> Optional[int]:
        return self.data.get(COLUMNS.contig_remapped_reads)

    @property
    def disease_status(self) -> Optional[str]:
        return self.data.get(COLUMNS.disease_status)

    @property
    def event_type(self) -> Optional[str]:
        return self.data.get(COLUMNS.event_type)

    @property
    def inferred_pairing(self) -> Optional[str]:
        return self.data.get(COLUMNS.inferred_pairing)

    @property
    def pairing(self) -> Optional[str]:
        return self.data.get(COLUMNS.pairing)

    @property
    def protocol(self) -> Optional[str]:
        return self.data.get(COLUMNS.protocol)

    @property
    def fusion_cdna_coding_start(self) -> Optional[int]:
        return self.data.get(COLUMNS.fusion_cdna_coding_start)

    @property
    def fusion_cdna_coding_end(self) -> Optional[int]:
        return self.data.get(COLUMNS.fusion_cdna_coding_end)

    @property
    def fusion_sequence_fasta_id(self) -> Optional[str]:
        return self.data.get(COLUMNS.fusion_sequence_fasta_id)

    @property
    def fusion_splicing_pattern(self) -> Optional[str]:
        return self.data.get(COLUMNS.fusion_splicing_pattern)

    @property
    def linking_split_reads(self) -> Optional[int]:
        return self.data.get(COLUMNS.linking_split_reads)

    @property
    def repeat_count(self) -> Optional[int]:
        return self.data.get(COLUMNS.repeat_count)

    @property
    def tracking_id(self) -> Optional[str]:
        return self.data.get(COLUMNS.tracking_id)

    @property
    def cluster_id(self) -> Optional[str]:
        return self.data.get(COLUMNS.cluster_id)

    @property
    def annotation_id(self) -> Optional[str]:
        return self.data.get(COLUMNS.annotation_id)

    @property
    def interchromosomal(self) -> bool:
        """bool: True if the breakpoints are on different chromosomes, False otherwise"""
        if self.break1.chr == self.break2.chr:
            return False
        return True

    @property
    def LL(self) -> bool:
        return self.break1.orient == ORIENT.LEFT and self.break2.orient == ORIENT.LEFT

    @property
    def LR(self) -> bool:
        return self.break1.orient == ORIENT.LEFT and self.break2.orient == ORIENT.RIGHT

    @property
    def RL(self) -> bool:
        return self.break1.orient == ORIENT.RIGHT and self.break2.orient == ORIENT.LEFT

    @property
    def RR(self) -> bool:
        return self.break1.orient == ORIENT.RIGHT and self.break2.orient == ORIENT.RIGHT

    def copy(self):
        temp = _copy(self)
        temp.data = {}
        temp.data.update(self.data)
        temp.break1 = _copy(self.break1)
        temp.break2 = _copy(self.break2)
        return temp

    def __init__(
        self,
        b1: Breakpoint,
        b2: Breakpoint,
        stranded: bool = False,
        opposing_strands: Optional[bool] = None,
        untemplated_seq: Optional[str] = None,
        **kwargs
    ):
        """
        Args:
            b1: the first breakpoint
            b2: the second breakpoint
            stranded: if not stranded then +/- is equivalent to -/+
            opposing_strands: are the strands at the breakpoint opposite? i.e. +/- instead of +/+
            untemplated_seq: seq between the breakpoints that is not part of either breakpoint
            data: optional dictionary of attributes associated with this pair

        Note:
            untemplated_seq should always be given wrt to the positive/forward reference strand

        Example:
            >>> BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
            >>> BreakpointPair(Breakpoint('1', 1, strand='+'), Breakpoint('1', 9999, strand='-'))
        """

        if b1.key[:3] > b2.key[:3]:
            self.break1 = b2
            self.break2 = b1
        else:
            self.break1 = b1
            self.break2 = b2
        self.stranded = stranded
        self.opposing_strands = opposing_strands  # type: ignore

        if self.break1.orient != ORIENT.NS and self.break2.orient != ORIENT.NS:
            if self.opposing_strands is not None:
                if (self.break1.orient == self.break2.orient and not self.opposing_strands) or (
                    self.break1.orient != self.break2.orient and self.opposing_strands
                ):
                    raise InvalidRearrangement(
                        'invalid breakpoint pair cannot form a valid combination',
                        b1,
                        b2,
                        self.opposing_strands,
                    )
            else:
                self.opposing_strands = self.break1.orient == self.break2.orient
        # between break1 and break2 not in either
        self.untemplated_seq = untemplated_seq
        self.data = kwargs

        if self.break1.strand != STRAND.NS and self.break2.strand != STRAND.NS:
            opposing = self.break1.strand != self.break2.strand
            if self.opposing_strands is None:
                self.opposing_strands = opposing
            elif self.opposing_strands != opposing:
                raise AssertionError(
                    'conflict in input arguments, opposing_strands must agree with input breakpoints'
                    ' when the strand has been specified'
                )

        if self.opposing_strands is None:
            raise NotSpecifiedError('must specify if opposing_strands')
        if self.stranded and STRAND.NS in [self.break1.strand, self.break2.strand]:
            raise NotSpecifiedError(
                'if stranded is specified, breakpoint strands cannot be \'not specified\''
            )

        # try classifying to make sure it's a valid combination
        BreakpointPair.classify(self)

    def column(self, colname: str):
        return self.data.get(COLUMNS[colname])

    def __str__(self):
        return 'BPP({}, {}{}{})'.format(
            str(self.break1),
            str(self.break2),
            ', opposing={}'.format(self.opposing_strands) if not self.stranded else '',
            ', seq=' + repr(self.untemplated_seq) if self.untemplated_seq is not None else '',
        )

    def flatten(self):
        """
        returns the key-value self for the breakpoint self information as
        can be written directly as a tab row
        """
        row = {}
        row.update(self.data)
        temp = {
            COLUMNS.break1_chromosome: self.break1.chr,
            COLUMNS.break1_position_start: self.break1.start,
            COLUMNS.break1_position_end: self.break1.end,
            COLUMNS.break1_orientation: self.break1.orient,
            COLUMNS.break1_strand: self.break1.strand,
            COLUMNS.break2_chromosome: self.break2.chr,
            COLUMNS.break2_position_start: self.break2.start,
            COLUMNS.break2_position_end: self.break2.end,
            COLUMNS.break2_orientation: self.break2.orient,
            COLUMNS.break2_strand: self.break2.strand,
            COLUMNS.opposing_strands: self.opposing_strands,
            COLUMNS.stranded: self.stranded,
            COLUMNS.untemplated_seq: self.untemplated_seq,
            COLUMNS.break1_seq: self.break1.seq,
            COLUMNS.break2_seq: self.break2.seq,
        }
        for col in temp:
            temp[col] = str(temp[col])
        row.update(temp)
        return row

    @classmethod
    def classify(cls, pair, distance: Optional[Callable] = None) -> Set[str]:
        """
        uses the chr, orientations and strands to determine the
        possible structural_variant types that this pair could support

        Args:
            pair (BreakpointPair): the pair to classify
            distance: if defined, will be passed to net size to use in narrowing the list of putative types (del vs ins)
        Returns:
            a list of possible SVTYPE

        Example:
            >>> bpp = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
            >>> BreakpointPair.classify(bpp)
            ['inversion']
            >>> bpp = BreakpointPair(Breakpoint('1', 1, orient='L'), Breakpoint('1', 9999, orient='R'), opposing_strands=False)
            >>> BreakpointPair.classify(bpp)
            {'deletion', 'insertion'}

        Note:
            see [related theory documentation](/background/theory/#classifying-events)
        """
        if not pair.interchromosomal:  # intrachromosomal
            if pair.opposing_strands:
                if pair.LR or pair.RL:
                    raise InvalidRearrangement(pair)
                return {SVTYPE.INV}
            else:
                if pair.LL or pair.RR:
                    raise InvalidRearrangement(pair)
                elif pair.break1.orient == ORIENT.LEFT or pair.break2.orient == ORIENT.RIGHT:
                    if (
                        len(pair.break1) == 1
                        and len(pair.break2) == 1
                        and abs(pair.break1.start - pair.break2.start) < 2
                    ):
                        if pair.untemplated_seq == '':
                            return set()
                        return {SVTYPE.INS}
                    elif pair.untemplated_seq == '':
                        return {SVTYPE.DEL}
                    elif distance:
                        try:
                            if pair.net_size(distance).start > 0:
                                return {SVTYPE.INS}
                            elif pair.net_size(distance).end < 0:
                                return {SVTYPE.DEL}
                        except ValueError:
                            pass
                    return {SVTYPE.DEL, SVTYPE.INS}
                elif pair.break1.orient == ORIENT.RIGHT or pair.break2.orient == ORIENT.LEFT:
                    return {SVTYPE.DUP}
                return {SVTYPE.DEL, SVTYPE.INS, SVTYPE.DUP}
        else:  # interchromosomal
            if pair.opposing_strands:
                if pair.LR or pair.RL:
                    raise InvalidRearrangement(pair)
                return {SVTYPE.ITRANS}
            else:
                if pair.LL or pair.RR:
                    raise InvalidRearrangement(pair)
                return {SVTYPE.TRANS}

    def net_size(self, distance=lambda x, y: Interval(abs(x - y))) -> Interval:
        """
        Returns the size of the event for a given pair. Mainly applicable to indels
        """
        if self.untemplated_seq is None or len(self.break1) > 1 or len(self.break2) > 1:
            raise ValueError('cannot determine net size of a non-specific breakpoint pair')
        if self.interchromosomal:
            return Interval(0)
        size = Interval(len(self.untemplated_seq))

        min_dist, max_dist = distance(self.break1.start, self.break2.start)
        # if it is a duplication the net change is the untemplated_seq as well as the duplicated region
        if self.RL:
            size.start += min_dist + 1
            size.end += max_dist + 1
        elif not self.opposing_strands:
            size.start -= max_dist - 1
            size.end -= min_dist - 1
        return size

    @property
    def is_putative_indel(self) -> bool:
        if self.interchromosomal or self.opposing_strands or self.break1.orient == ORIENT.RIGHT:
            return False
        return True

    def breakpoint_sequence_homology(self, reference_genome: ReferenceGenome):
        """
        for a given set of breakpoints matches the sequence opposite the partner breakpoint
        this sequence comparison is done with reference to a reference genome and does not
        use novel or untemplated sequence in the comparison. For this reason, insertions
        will never return any homologous sequence


            small duplication event CTT => CTTCTT

            GATACATTTCTTCTTGAAAA reference
            ---------<========== first breakpoint
            ===========>-------- second breakpoint
            ---------CT-CT------ first break homology
            -------TT-TT-------- second break homology

        Args:
            reference_genome: dict of reference sequence by template/chr name

        Returns:
            Tuple[str,str]: homologous sequence at the first breakpoint and second breakpoints

        Raises:
            AttributeError: for non specific breakpoints
        """

        b1_refseq = reference_genome[self.break1.chr].seq
        b2_refseq = reference_genome[self.break2.chr].seq
        if len(self.break1) > 1 or len(self.break2) > 1:
            raise AttributeError('cannot call shared sequence for non-specific breakpoints')

        # find for the first breakpoint
        pos1 = self.break1.start - 1
        pos2 = self.break2.start - 1
        shift1 = -1 if self.break1.orient == ORIENT.LEFT else 1
        shift2 = -1 if self.break2.orient == ORIENT.RIGHT else 1
        pos2 += shift2

        first_seq = []
        while all(
            [
                pos1 >= 0,
                pos1 < len(b1_refseq),
                pos2 >= 0,
                pos2 < len(b2_refseq),
                self.interchromosomal or pos2 > self.break1.start - 1,
                self.interchromosomal or pos1 < self.break2.start - 1,
            ]
        ):
            break1_bp = (
                b1_refseq[pos1]
                if self.break1.strand == STRAND.POS
                else reverse_complement(b1_refseq[pos1])
            )
            break2_bp = (
                b2_refseq[pos2]
                if self.break2.strand == STRAND.POS
                else reverse_complement(b2_refseq[pos2])
            )
            if DNA_ALPHABET.match(break1_bp, break2_bp):
                first_seq.append(b1_refseq[pos1])
            else:
                break
            pos1 += shift1
            pos2 += shift2

        if shift1 < 0:
            first_seq.reverse()
        # now go over the second breakpoint
        pos1 = self.break1.start - 1  # reset the start points
        pos2 = self.break2.start - 1
        shift1 *= -1  # flip the directions
        shift2 *= -1
        pos1 += shift1

        second_seq = []
        while all(
            [
                pos1 >= 0,
                pos1 < len(b1_refseq),
                pos2 >= 0,
                pos2 < len(b2_refseq),
                self.interchromosomal or pos2 > self.break1.start - 1,
                self.interchromosomal or pos1 < self.break2.start - 1,
            ]
        ):
            break1_bp = (
                b1_refseq[pos1]
                if self.break1.strand == STRAND.POS
                else reverse_complement(b1_refseq[pos1])
            )
            break2_bp = (
                b2_refseq[pos2]
                if self.break2.strand == STRAND.POS
                else reverse_complement(b2_refseq[pos2])
            )
            if DNA_ALPHABET.match(break1_bp, break2_bp):
                second_seq.append(b2_refseq[pos2])
            else:
                break
            pos1 += shift1
            pos2 += shift2

        if shift1 < 0:
            second_seq.reverse()

        return ''.join(first_seq).upper(), ''.join(second_seq).upper()

    @staticmethod
    def _breakpoint_utemp_shift(breakpoint, untemplated_seq, reference_genome):
        """
        Calculates the amount of the untemplated sequence that could be aligned further.
        Gives an idea of how much different alignments might differ.
        """
        break_shift = 0
        lutemp = len(untemplated_seq)
        if breakpoint.orient == ORIENT.LEFT:
            break_seq = reference_genome[breakpoint.chr].seq[
                breakpoint.start - lutemp : breakpoint.start
            ]
            for i in range(1, lutemp + 1):
                if break_seq[lutemp - i] == untemplated_seq[lutemp - i]:
                    break_shift += 1
                else:
                    break
        else:
            break_seq = reference_genome[breakpoint.chr].seq[
                breakpoint.start - 1 : breakpoint.start + lutemp - 1
            ]
            for i in range(0, lutemp):
                if break_seq[i] == untemplated_seq[i]:
                    break_shift += 1
                else:
                    break
        return break_shift

    def untemplated_shift(self, reference_genome):
        """gives a range for each breakpoint on the possible alignment range in the shifting the untemplated
        sequence"""
        if any(
            [
                self.untemplated_seq is None,
                len(self.break1) + len(self.break2) > 2,
                self.break1.orient == ORIENT.NS,
                self.break2.orient == ORIENT.NS,
            ]
        ):
            raise AttributeError(
                'Cannot calculate the shift for non specific breakpoint calls', self
            )
        break1_shift = BreakpointPair._breakpoint_utemp_shift(
            self.break1, self.untemplated_seq, reference_genome
        )
        break2_shift = BreakpointPair._breakpoint_utemp_shift(
            self.break2, self.untemplated_seq, reference_genome
        )
        return (break2_shift, break1_shift)

    def get_bed_repesentation(self) -> List[Tuple[str, int, int, Optional[str]]]:
        bed = []
        if self.interchromosomal:
            bed.append(
                (
                    self.break1.chr,
                    self.break1.start,
                    self.break1.end,
                    self.data.get(COLUMNS.cluster_id, None),
                )
            )
            bed.append(
                (
                    self.break2.chr,
                    self.break2.start,
                    self.break2.end,
                    self.data.get(COLUMNS.cluster_id, None),
                )
            )
        else:
            bed.append(
                (
                    self.break1.chr,
                    self.break1.start,
                    self.break2.end,
                    self.data.get(COLUMNS.cluster_id, None),
                )
            )
        return bed
