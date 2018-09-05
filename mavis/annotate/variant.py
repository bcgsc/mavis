import itertools
import json
from shortuuid import uuid

from .fusion import determine_prime, FusionTranscript
from .genomic import IntergenicRegion
from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import COLUMNS, GENE_PRODUCT_TYPE, PROTOCOL, STOP_AA, STRAND, SVTYPE
from ..error import NotSpecifiedError
from ..interval import Interval
from ..util import DEVNULL


class Annotation(BreakpointPair):
    """
    a fusion of two transcripts created by the associated breakpoint_pair
    will also hold the other annotations for overlapping and encompassed and nearest genes
    """

    def __init__(
        self, bpp,
        transcript1=None,
        transcript2=None,
        proximity=5000,
        data=None,
        **kwargs
    ):
        """
        Holds a breakpoint call and a set of transcripts, other information is gathered relative to these

        Args:
            bpp (BreakpointPair): the breakpoint pair call. Will be adjusted and then stored based on the transcripts
            transcript1 (Transcript): transcript at the first breakpoint
            transcript2 (Transcript): Transcript at the second breakpoint
            data (dict): optional dictionary to hold related attributes
            event_type (SVTYPE): the type of event
        """
        # narrow the breakpoint windows by the transcripts being used for annotation
        temp = bpp.break1 if transcript1 is None else bpp.break1 & transcript1
        b1 = Breakpoint(bpp.break1.chr, temp[0], temp[1], strand=bpp.break1.strand, orient=bpp.break1.orient)

        temp = bpp.break2 if transcript2 is None else bpp.break2 & transcript2
        b2 = Breakpoint(bpp.break2.chr, temp[0], temp[1], strand=bpp.break2.strand, orient=bpp.break2.orient)
        BreakpointPair.__init__(
            self,
            b1,
            b2,
            opposing_strands=bpp.opposing_strands,
            stranded=bpp.stranded,
            untemplated_seq=bpp.untemplated_seq
        )
        self.data.update(bpp.data)
        if data is not None:
            conflicts = set(kwargs.keys()) & set(data.keys())
            self.data.update(data)
            if conflicts:
                raise TypeError('got multiple values for data elements:', conflicts)
        self.data.update(kwargs)

        self.transcript1 = transcript1
        self.transcript2 = transcript2

        self.encompassed_genes = set()
        self.genes_proximal_to_break1 = set()
        self.genes_proximal_to_break2 = set()
        self.genes_overlapping_break1 = set()
        self.genes_overlapping_break2 = set()

        SVTYPE.enforce(self.event_type)
        PROTOCOL.enforce(self.protocol)
        self.proximity = proximity
        self.fusion = None

    def add_gene(self, input_gene):
        """
        adds a input_gene to the current set of annotations. Checks which set it should be added to

        Args:
            input_gene (input_gene): the input_gene being added
        """
        if input_gene.chr not in [self.break1.chr, self.break2.chr]:
            raise AttributeError('cannot add input_gene not on the same chromosome as either breakpoint')

        if not self.interchromosomal:
            try:
                encompassment = Interval(self.break1.end + 1, self.break2.start - 1)
                if input_gene in encompassment:
                    self.encompassed_genes.add(input_gene)
            except AttributeError:
                pass
        if Interval.overlaps(input_gene, self.break1) and input_gene.chr == self.break1.chr \
                and input_gene != self.transcript1.reference_object:
            self.genes_overlapping_break1.add(input_gene)
        if Interval.overlaps(input_gene, self.break2) and input_gene.chr == self.break2.chr \
                and input_gene != self.transcript2.reference_object:
            self.genes_overlapping_break2.add(input_gene)

        if input_gene in self.genes_overlapping_break1 or input_gene in self.genes_overlapping_break2 or \
                input_gene in self.encompassed_genes or input_gene == self.transcript1.reference_object or \
                input_gene == self.transcript2.reference_object:
            return

        dist1 = Interval.dist(input_gene, self.break1)
        dist2 = Interval.dist(input_gene, self.break2)

        if self.interchromosomal:
            if input_gene.chr == self.break1.chr:
                self.genes_proximal_to_break1.add((input_gene, dist1))
            elif input_gene.chr == self.break2.chr:
                self.genes_proximal_to_break2.add((input_gene, dist2))
        else:
            if dist1 < 0:
                self.genes_proximal_to_break1.add((input_gene, dist1))
            if dist2 > 0:
                self.genes_proximal_to_break2.add((input_gene, dist2))

        if self.genes_proximal_to_break1:
            temp = set()
            tgt = min([abs(d) for g, d in self.genes_proximal_to_break1])

            for gene, dist in self.genes_proximal_to_break1:
                if self.proximity is None:
                    if abs(dist) == tgt:
                        temp.add((gene, dist))
                elif abs(dist) <= self.proximity:
                    temp.add((gene, dist))

            self.genes_proximal_to_break1 = temp

        if self.genes_proximal_to_break2:
            temp = set()
            tgt = min([abs(d) for g, d in self.genes_proximal_to_break2])

            for gene, dist in self.genes_proximal_to_break2:
                if self.proximity is None:
                    if abs(dist) == tgt:
                        temp.add((gene, dist))
                elif abs(dist) <= self.proximity:
                    temp.add((gene, dist))

            self.genes_proximal_to_break2 = temp

    def flatten(self):
        """
        generates a dictionary of the annotation information as strings

        Returns:
            :class:`dict` of :class:`str` by :class:`str`: dictionary of attribute names and values
        """
        row = BreakpointPair.flatten(self)
        row.update({
            COLUMNS.genes_proximal_to_break1: self.genes_proximal_to_break1,
            COLUMNS.genes_proximal_to_break2: self.genes_proximal_to_break2,
            COLUMNS.gene1_direction: None,
            COLUMNS.gene2_direction: None,
            COLUMNS.gene_product_type: None,
            COLUMNS.gene1: None,
            COLUMNS.gene2: None,
            COLUMNS.transcript1: '{}:{}_{}{}'.format(
                self.transcript1.reference_object,
                self.transcript1.start,
                self.transcript1.end,
                self.transcript1.get_strand()),
            COLUMNS.transcript2: '{}:{}_{}{}'.format(
                self.transcript2.reference_object,
                self.transcript2.start,
                self.transcript2.end,
                self.transcript2.get_strand()),
            COLUMNS.genes_encompassed: ';'.join(sorted([x.name for x in self.encompassed_genes])),
            COLUMNS.genes_overlapping_break1: ';'.join(sorted([x.name for x in self.genes_overlapping_break1])),
            COLUMNS.genes_overlapping_break2: ';'.join(sorted([x.name for x in self.genes_overlapping_break2])),
            COLUMNS.genes_proximal_to_break1: ';'.join(
                sorted(['{}({})'.format(x[0].name, x[1]) for x in self.genes_proximal_to_break1])),
            COLUMNS.genes_proximal_to_break2: ';'.join(
                sorted(['{}({})'.format(x[0].name, x[1]) for x in self.genes_proximal_to_break2])),
            COLUMNS.event_type: self.event_type,
            COLUMNS.gene1_aliases: None,
            COLUMNS.gene2_aliases: None
        })
        if hasattr(self.transcript1, 'gene'):
            row[COLUMNS.gene1] = self.transcript1.gene.name
            row[COLUMNS.transcript1] = self.transcript1.name
            if self.transcript1.gene.aliases:
                row[COLUMNS.gene1_aliases] = ';'.join(sorted(self.transcript1.gene.aliases))
            try:
                row[COLUMNS.gene1_direction] = str(determine_prime(self.transcript1, self.break1))
            except NotSpecifiedError:
                pass
        if hasattr(self.transcript2, 'gene'):
            row[COLUMNS.gene2] = self.transcript2.gene.name
            row[COLUMNS.transcript2] = self.transcript2.name
            if self.transcript2.gene.aliases:
                row[COLUMNS.gene2_aliases] = ';'.join(sorted(self.transcript2.gene.aliases))
            try:
                row[COLUMNS.gene2_direction] = str(determine_prime(self.transcript2, self.break2))
                if row[COLUMNS.gene1_direction] is not None:
                    if row[COLUMNS.gene1_direction] == row[COLUMNS.gene2_direction]:
                        row[COLUMNS.gene_product_type] = GENE_PRODUCT_TYPE.ANTI_SENSE
                    else:
                        row[COLUMNS.gene_product_type] = GENE_PRODUCT_TYPE.SENSE
            except NotSpecifiedError:
                pass
        return row

    def single_transcript(self):
        return bool(self.transcript1 == self.transcript2 and self.transcript1)


def flatten_fusion_translation(translation):
    """
    for a given fusion product (translation) gather the information to be output to the tabbed files

    Args:
        translation (Translation): the translation which is on the fusion transcript
    Returns:
        dict: the dictionary of column names to values
    """
    row = dict()
    row[COLUMNS.fusion_cdna_coding_start] = translation.start
    row[COLUMNS.fusion_cdna_coding_end] = translation.end

    # select the exon that has changed
    domains = []
    for dom in translation.domains:
        m, t = dom.score_region_mapping()
        temp = {
            'name': dom.name,
            'sequences': dom.get_seqs(),
            'regions': [
                {'start': dr.start, 'end': dr.end} for dr in sorted(dom.regions, key=lambda x: x.start)
            ],
            'mapping_quality': round(m * 100 / t, 0),
            'matches': m
        }
        domains.append(temp)
    row[COLUMNS.fusion_mapped_domains] = json.dumps(domains)
    return row


class IndelCall:
    def __init__(self, refseq, mutseq):
        """
        Given two sequences, Assuming there exists a single difference between the two
        call an indel which accounts for the change

        Args:
            refseq (str): The reference (amino acid) sequence
            mutseq (str): The mutated (amino acid) sequence

        Attributes:
            nterm_aligned (int): the number of characters aligned consecutively from the start of both strings
            cterm_aligned (int): the number of characters aligned consecutively from the end of both strings
            is_dup (bool): flag to indicate a duplication
            ref_seq (str): the reference sequence
            mut_seq (str): the mutated sequence
            ins_seq (str): the inserted sequence
            del_seq (str): the deleted sequence
            terminates (bool): both sequences end in stop AAs
        """
        self.nterm_aligned = 0
        self.cterm_aligned = 0
        self.ref_seq = refseq
        self.mut_seq = mutseq
        self.is_dup = False
        if self.ref_seq[-1] == STOP_AA and self.mut_seq[-1] == STOP_AA:
            self.terminates = True
            self.ref_seq = self.ref_seq[:-1]
            self.mut_seq = self.mut_seq[:-1]
        else:
            self.terminates = False

        min_len = min(len(self.mut_seq), len(self.ref_seq))

        if self.mut_seq[:min_len] == self.ref_seq[:min_len]:
            # full n-terminal match
            self.nterm_aligned = min_len
            if len(self.mut_seq) > len(self.ref_seq):
                self.ins_seq = self.mut_seq[min_len:]
                self.del_seq = ''
            else:
                self.del_seq = self.ref_seq[min_len:]
                self.ins_seq = ''

        elif self.mut_seq[-1 * min_len:] == self.ref_seq[-1 * min_len:]:
            # full c-terminal match
            self.cterm_aligned = min_len
            if len(self.mut_seq) > len(self.ref_seq):
                self.ins_seq = self.mut_seq[:0 - min_len]
                self.del_seq = ''
            else:
                self.del_seq = self.ref_seq[:0 - min_len]
                self.ins_seq = ''
        else:
            for pos in range(0, min_len):
                if self.ref_seq[pos] != self.mut_seq[pos]:
                    break
                self.nterm_aligned = pos + 1

            for pos in range(0, min_len):
                if self.ref_seq[-1 - pos] != self.mut_seq[-1 - pos]:
                    break
                self.cterm_aligned = pos + 1

            if not self.cterm_aligned:
                self.del_seq = self.ref_seq[self.nterm_aligned:]
                self.ins_seq = self.mut_seq[self.nterm_aligned:]

            elif not self.nterm_aligned:
                self.del_seq = self.ref_seq[:0 - self.cterm_aligned]
                self.ins_seq = self.mut_seq[:0 - self.cterm_aligned]

            elif len(self.ref_seq) - self.cterm_aligned + 1 <= self.nterm_aligned:

                # repeat region
                diff = len(self.mut_seq) - len(self.ref_seq)
                if diff > 0:
                    ins_length = diff
                    del_length = 0
                else:
                    del_length = abs(diff)
                    ins_length = 0
                self.ins_seq = mutseq[self.nterm_aligned:self.nterm_aligned + ins_length]
                self.del_seq = refseq[self.nterm_aligned:self.nterm_aligned + del_length]

                if self.ins_seq:
                    repeat_start = max(
                        len(self.ref_seq) - self.cterm_aligned,
                        self.nterm_aligned - len(self.ins_seq)
                    )
                    dupped_refseq = self.ref_seq[repeat_start:self.nterm_aligned]

                    if dupped_refseq == self.ins_seq:
                        self.is_dup = True

            else:
                # regular indel
                self.del_seq = self.ref_seq[self.nterm_aligned:0 - self.cterm_aligned]
                self.ins_seq = self.mut_seq[self.nterm_aligned:0 - self.cterm_aligned]

    def hgvs_protein_notation(self):
        """
        returns the HGVS protein notation for an indel call
        """
        if not self.ins_seq and not self.del_seq:  # synonymous variant
            return None

        if not self.cterm_aligned and self.ins_seq:
            if not self.del_seq:
                # c-terminal extension
                notation = 'p.{}{}ext{}{}'.format(
                    self.ref_seq[-1] if not self.terminates else STOP_AA,
                    len(self.ref_seq) + self.terminates,
                    STOP_AA if self.mut_seq[-1] == STOP_AA or self.terminates else '',
                    len(self.mut_seq) - len(self.ref_seq)
                )
            else:
                # frameshift indel
                notation = 'p.{}{}{}fs'.format(
                    self.ref_seq[self.nterm_aligned],
                    self.nterm_aligned + 1,
                    self.ins_seq[0]
                )
                if STOP_AA in self.ins_seq:
                    notation += '*{}'.format(self.ins_seq.index(STOP_AA) + 1)
                elif self.terminates:
                    notation += '*{}'.format(len(self.ins_seq) + 1)
        elif not self.nterm_aligned and self.ins_seq and not self.del_seq:
            # n-terminal extension
            notation = 'p.{}1ext-{}'.format(
                self.ref_seq[0],
                len(self.mut_seq) - len(self.ref_seq)
            )
        elif self.is_dup:
            if self.del_seq:
                raise NotImplementedError('duplication/deletion no supported', self)

            dup_start = self.nterm_aligned - len(self.ins_seq) + 1
            if dup_start == self.nterm_aligned:
                notation = 'p.{}{}dup{}'.format(self.ref_seq[self.nterm_aligned - 1], self.nterm_aligned, self.ins_seq)
            else:
                notation = 'p.{}{}_{}{}dup{}'.format(
                    self.ref_seq[dup_start - 1], dup_start, self.ref_seq[self.nterm_aligned - 1], self.nterm_aligned, self.ins_seq)
        else:
            if self.del_seq:  # indel
                notation = 'p.{}{}'.format(
                    self.ref_seq[self.nterm_aligned],
                    self.nterm_aligned + 1
                )
                if len(self.del_seq) > 1:
                    notation += '_{}{}'.format(
                        self.ref_seq[self.nterm_aligned + len(self.del_seq) - 1],
                        self.nterm_aligned + len(self.del_seq)
                    )
                notation += 'del{}'.format(self.del_seq)

                if self.ins_seq:
                    notation += 'ins{}'.format(self.ins_seq)

            else:  # insertion
                notation = 'p.{}{}_{}{}ins{}'.format(
                    self.ref_seq[self.nterm_aligned - 1],
                    self.nterm_aligned,
                    self.ref_seq[self.nterm_aligned],
                    self.nterm_aligned + 1,
                    self.ins_seq
                )

        return notation

    def __str__(self):
        return 'IndelCall({})'.format(', '.join(['{}={}'.format(k, repr(v)) for k, v in sorted(self.__dict__.items())]))


def call_protein_indel(ref_translation, fusion_translation, reference_genome=None):
    """
    compare the fusion protein/aa sequence to the reference protein/aa sequence and
    return an hgvs notation indel call

    Args:
        ref_translation (Translation): the reference protein/translation
        fusion_translation (Translation): the fusion protein/translation
        reference_genome: the reference genome object used to fetch the reference translation AA sequence
    Returns:
        str: the :term:`HGVS` protein indel notation
    """
    ref_aa_seq = ref_translation.get_aa_seq(reference_genome)
    call = IndelCall(ref_aa_seq, fusion_translation.get_aa_seq())
    notation = call.hgvs_protein_notation()
    if not notation:
        return None
    name = ref_translation.name
    curr = ref_translation
    while name is None:
        curr = curr.reference_object
        name = curr.name
    return '{}:{}'.format(name, notation)


def flatten_fusion_transcript(spliced_fusion_transcript):
    row = {}
    five_prime_exons = []
    three_prime_exons = []
    fusion_transcript = spliced_fusion_transcript.unspliced_transcript
    for ex in spliced_fusion_transcript.exons:
        try:
            src_exon = fusion_transcript.exon_mapping[ex.position]
            number = src_exon.transcript.exon_number(src_exon)
            if ex.end <= fusion_transcript.break1:
                five_prime_exons.append(number)
            elif ex.start >= fusion_transcript.break2:
                three_prime_exons.append(number)
            else:
                raise AssertionError(
                    'exon should not be mapped if not within a break region',
                    ex, fusion_transcript.break1, fusion_transcript.break2
                )
        except KeyError:  # novel exon
            for us_exon, src_exon in sorted(fusion_transcript.exon_mapping.items()):
                if Interval.overlaps(ex, us_exon):
                    number = src_exon.transcript.exon_number(src_exon)
                    if us_exon.end <= fusion_transcript.break1:
                        five_prime_exons.append(number)
                    elif us_exon.start >= fusion_transcript.break2:
                        three_prime_exons.append(number)
                    else:
                        raise AssertionError(
                            'exon should not be mapped if not within a break region',
                            us_exon, fusion_transcript.break1. fusion_transcript.break2
                        )
    row[COLUMNS.exon_last_5prime] = five_prime_exons[-1]
    row[COLUMNS.exon_first_3prime] = three_prime_exons[0]
    row[COLUMNS.fusion_splicing_pattern] = spliced_fusion_transcript.splicing_pattern.splice_type
    return row


def overlapping_transcripts(ref_ann, breakpoint):
    """
    Args:
        ref_ann (:class:`dict` of :class:`list` of :any:`Gene` by :class:`str`): the reference list of genes split
            by chromosome
        breakpoint (Breakpoint): the breakpoint in question
    Returns:
        :class:`list` of :any:`PreTranscript`: a list of possible transcripts
    """
    putative_annotations = set()
    for gene in ref_ann.get(breakpoint.chr, []):
        for transcript in gene.transcripts:
            if breakpoint.strand != STRAND.NS and transcript.get_strand() != STRAND.NS \
                    and transcript.get_strand() != breakpoint.strand:
                continue
            if Interval.overlaps(breakpoint, transcript):
                putative_annotations.add(transcript)
    return putative_annotations


def _gather_breakpoint_annotations(ref_ann, breakpoint):
    """
    Args:
        ref_ann (:class:`dict` of :class:`list` of :class:`Gene` by :class:`str`): the reference annotations split
            into lists of genes by chromosome
        breakpoint (Breakpoint): the breakpoint annotations are to be gathered for

    Returns:
        tuple: tuple contains

            - :class:`list` of (:class:`PreTranscript` or :class:`IntergenicRegion`): transcripts or intergenic regions
              overlapping the breakpoint on the positive strand
            - :class:`list` of (:class:`PreTranscript` or :class:`IntergenicRegion`): transcripts or intergenic regions
              overlapping the breakpoint on the negative strand

    .. todo::

        Support for setting the transcript in the annotation when the breakpoint is just ahead of the transcript
        and the transcript would be 3'. Then assuming the splicing model takes the 2nd exon onward
    """

    pos_overlapping_transcripts = []
    neg_overlapping_transcripts = []
    for gene in ref_ann.get(breakpoint.chr, []):
        for t in gene.transcripts:
            if Interval.overlaps(t, breakpoint):
                if STRAND.compare(t.get_strand(), STRAND.POS):
                    pos_overlapping_transcripts.append(t)
                if STRAND.compare(t.get_strand(), STRAND.NEG):
                    neg_overlapping_transcripts.append(t)

    pos_intervals = Interval.min_nonoverlapping(*pos_overlapping_transcripts)
    neg_intervals = Interval.min_nonoverlapping(*neg_overlapping_transcripts)

    temp = []
    # before the first?
    if len(pos_intervals) > 0:
        first = pos_intervals[0]
        last = pos_intervals[-1]
        if breakpoint.start < first.start:
            temp.append(IntergenicRegion(breakpoint.chr, breakpoint[0], first[0] - 1, STRAND.POS))
        if breakpoint.end > last.end:
            temp.append(IntergenicRegion(breakpoint.chr, last[1] + 1, breakpoint[1], STRAND.POS))

        for i, curr in enumerate(pos_intervals):
            if i > 0:
                prev = pos_intervals[i - 1]
                try:
                    temp.append(IntergenicRegion(breakpoint.chr, prev[1] + 1, curr[0] - 1, STRAND.POS))
                except AttributeError:
                    pass
    else:
        temp.append(IntergenicRegion(breakpoint.chr, breakpoint.start, breakpoint.end, STRAND.POS))
    pos_overlapping_transcripts.extend(temp)

    temp = []
    # before the first?
    if len(neg_intervals) > 0:
        first = neg_intervals[0]
        last = neg_intervals[-1]
        if breakpoint.start < first.start:
            temp.append(IntergenicRegion(breakpoint.chr, breakpoint[0], first[0] - 1, STRAND.NEG))
        if breakpoint.end > last.end:
            temp.append(IntergenicRegion(breakpoint.chr, last[1] + 1, breakpoint[1], STRAND.NEG))

        for i, curr in enumerate(neg_intervals):
            if i > 0:
                prev = neg_intervals[i - 1]
                try:
                    temp.append(IntergenicRegion(breakpoint.chr, prev[1] + 1, curr[0] - 1, STRAND.NEG))
                except AttributeError:
                    pass
    else:
        temp.append(IntergenicRegion(breakpoint.chr, breakpoint.start, breakpoint.end, STRAND.NEG))
    neg_overlapping_transcripts.extend(temp)

    if len(pos_overlapping_transcripts) == 0:
        raise AssertionError('neither strand group should ever be empty', pos_overlapping_transcripts)
    if len(neg_overlapping_transcripts) == 0:
        raise AssertionError('neither strand group should ever be empty', neg_overlapping_transcripts)
    return (
        sorted(pos_overlapping_transcripts, key=lambda x: x.position),
        sorted(neg_overlapping_transcripts, key=lambda x: x.position))


def _gather_annotations(ref, bp, proximity=None):
    """
    each annotation is defined by the annotations selected at the breakpoints
    the other annotations are given relative to this
    the annotation at the breakpoint can be a transcript or an intergenic region

    Args:
        ref (:class:`dict` of :class:`list` of :any:`Gene` by :class:`str`): the list of reference genes hashed
            by chromosomes
        breakpoint_pairs (:class:`list` of :any:`BreakpointPair`): breakpoint pairs we wish to annotate as events

    Returns:
        :class:`list` of :class:`Annotation`: The annotations
    """
    annotations = dict()
    break1_pos, break1_neg = _gather_breakpoint_annotations(ref, bp.break1)
    break2_pos, break2_neg = _gather_breakpoint_annotations(ref, bp.break2)

    combinations = []

    if bp.stranded:
        if bp.break1.strand == STRAND.POS:
            if bp.break2.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_pos, break2_pos))
            else:
                combinations.extend(itertools.product(break1_pos, break2_neg))
        else:
            if bp.break2.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_neg, break2_pos))
            else:
                combinations.extend(itertools.product(break1_neg, break2_neg))
    else:
        # single transcript starts ....
        for t in (set(break1_pos) | set(break1_neg)) & (set(break2_pos) | set(break2_neg)):
            try:
                t.gene
            except AttributeError:
                pass
            else:
                combinations.append((t, t))
        if bp.opposing_strands:
            combinations.extend(itertools.product(break1_pos, break2_neg))
            combinations.extend(itertools.product(break1_neg, break2_pos))
        else:
            combinations.extend(itertools.product(break1_pos, break2_pos))
            combinations.extend(itertools.product(break1_neg, break2_neg))

    same = set()
    for a1, a2 in combinations:
        """
        if a1 != a2 and hasattr(a1, 'exons') != hasattr(a2, 'exons') and not bp.interchromosomal:
            # one is a transcript, the other an intergenic region
            # take the transcript if it covers both breakpoints
            # this is due to the special case 'single transcript inversion'
            if hasattr(a1, 'exons'):
                if Interval.overlaps(bp.break1, a1) and Interval.overlaps(bp.break2, a1):
                    a2 = a1
            else:
                if Interval.overlaps(bp.break1, a2) and Interval.overlaps(bp.break2, a2):
                    a1 = a2
        """
        if (a1, a2) in annotations:  # ignore duplicates
            continue
        try:
            if a1.gene == a2.gene and a1 != a2:
                continue
        except AttributeError:
            pass

        if a1 == a2 and hasattr(a1, 'exons'):
            same.add(a1)

        b1_itvl = bp.break1 & a1
        b2_itvl = bp.break2 & a2
        bpp = BreakpointPair.copy(bp)
        bpp.break1.start = b1_itvl[0]
        bpp.break1.end = b1_itvl[1]
        bpp.break2.start = b2_itvl[0]
        bpp.break2.end = b2_itvl[1]

        a = Annotation(bpp, a1, a2, proximity=proximity)

        for gene in ref.get(bp.break1.chr, []):
            a.add_gene(gene)
        if bp.interchromosomal:
            for gene in ref.get(bp.break2.chr, []):
                a.add_gene(gene)
        annotations[(a1, a2)] = a
    filtered = []  # remove any inter-gene/inter-region annotations where a same transcript was found
    for pair, ann in annotations.items():
        a1, a2 = pair
        if (a1 in same or a2 in same) and a1 != a2:
            pass
        else:
            filtered.append(ann)
    return filtered


def choose_more_annotated(ann_list):
    """
    for a given set of annotations if there are annotations which contain transcripts and
    annotations that are simply intergenic regions, discard the intergenic region annotations

    similarly if there are annotations where both breakpoints fall in a transcript and
    annotations where one or more breakpoints lands in an intergenic region, discard those
    that land in the intergenic region

    Args:
        ann_list (list of :class:`Annotation`): list of input annotations

    Warning:
        input annotations are assumed to be the same event (the same validation_id)
        the logic used would not apply to different events

    Returns:
        list of :class:`Annotation`: the filtered list
    """
    two_transcript = []
    one_transcript = []
    intergenic = []

    for ann in ann_list:
        if isinstance(ann.transcript1, IntergenicRegion) and isinstance(ann.transcript2, IntergenicRegion):
            intergenic.append(ann)
        elif isinstance(ann.transcript1, IntergenicRegion) or isinstance(ann.transcript2, IntergenicRegion):
            one_transcript.append(ann)
        else:
            two_transcript.append(ann)

    if len(two_transcript) > 0:
        return two_transcript
    elif len(one_transcript) > 0:
        return one_transcript
    else:
        return intergenic


def choose_transcripts_by_priority(ann_list):
    """
    for each set of annotations with the same combinations of genes, choose the
    annotation with the most "best_transcripts" or most "alphanumeric" choices
    of transcript. Throw an error if they are identical

    Args:
        ann_list (list of :class:`Annotation`): input annotations

    Warning:
        input annotations are assumed to be the same event (the same validation_id)
        the logic used would not apply to different events

    Returns:
        list of :class:`Annotation`: the filtered list
    """
    annotations_by_gene_combination = {}
    genes = set()

    for ann in ann_list:
        gene1 = None
        gene2 = None
        try:
            gene1 = ann.transcript1.gene
            genes.add(gene1)
        except AttributeError:
            pass
        try:
            gene2 = ann.transcript2.gene
            genes.add(gene2)
        except AttributeError:
            pass
        annotations_by_gene_combination.setdefault((gene1, gene2), []).append(ann)

    filtered_annotations = []
    for g, sublist in annotations_by_gene_combination.items():
        gene1, gene2 = g
        if gene1 is None and gene2 is None:
            filtered_annotations.extend(sublist)
        elif gene2 is None:
            ann = min(sublist, key=lambda a: gene1.transcript_priority(a.transcript1))
            filtered_annotations.append(ann)
        elif gene1 is None:
            ann = min(sublist, key=lambda a: gene2.transcript_priority(a.transcript2))
            filtered_annotations.append(ann)
        else:
            ann = min(sublist, key=lambda a: (
                gene1.transcript_priority(a.transcript1) + gene2.transcript_priority(a.transcript2),
                gene1.transcript_priority(a.transcript1),
                gene2.transcript_priority(a.transcript2)
            ))
            filtered_annotations.append(ann)
    return filtered_annotations


def annotate_events(
    bpps,
    annotations,
    reference_genome,
    max_proximity=5000,
    min_orf_size=200,
    min_domain_mapping_match=0.95,
    max_orf_cap=3,
    log=DEVNULL,
    filters=None
):
    """
    Args:
        bpps (list of :class:`~mavis.breakpoint.BreakpointPair`): list of events
        annotations: reference annotations
        reference_genome (dict of string by string): dictionary of reference sequences by name
        max_proximity (int): see :term:`max_proximity`
        min_orf_size (int): see :term:`min_orf_size`
        min_domain_mapping_match (float): see :term:`min_domain_mapping_match`
        max_orf_cap (int): see :term:`max_orf_cap`
        log (callable): callable function to take in strings and time_stamp args
        filters (list of callable): list of functions taking in a list and returning a list for filtering

    Returns:
        list of :class:`Annotation`: list of the putative annotations
    """
    if filters is None:
        filters = [choose_more_annotated, choose_transcripts_by_priority]
    results = []
    total = len(bpps)
    for i, bpp in enumerate(bpps):
        log('({} of {}) gathering annotations for'.format(i + 1, total), bpp)
        bpp.data[COLUMNS.validation_id] = bpp.data.get(COLUMNS.validation_id, str(uuid()))
        ann_list = _gather_annotations(
            annotations,
            bpp,
            proximity=max_proximity
        )
        for f in filters:
            ann_list = f(ann_list)  # apply the filter
        results.extend(ann_list)
        for j, ann in enumerate(ann_list):
            ann.data[COLUMNS.annotation_id] = '{}-a{}'.format(ann.validation_id, j + 1)
            if ann.untemplated_seq is None:
                if len(ann.break1) == 1 and len(ann.break2) == 1 and ann.event_type != SVTYPE.INS:
                    ann.untemplated_seq = ''
                    ann.data[COLUMNS.assumed_untemplated] = True
            else:
                ann.data[COLUMNS.assumed_untemplated] = False
            # try building the fusion product
            try:
                ft = FusionTranscript.build(
                    ann, reference_genome,
                    min_orf_size=min_orf_size,
                    max_orf_cap=max_orf_cap,
                    min_domain_mapping_match=min_domain_mapping_match
                )
                ann.fusion = ft
            except NotSpecifiedError:
                pass  # shouldn't build fusions for non-specific calls anyway
            except AttributeError:
                pass  # will be thrown when transcript1/2 are intergenic ranges and not actual transcripts
            except NotImplementedError:
                pass  # anti-sense fusions will throw this error
            except KeyError as e:
                log('warning. could not build fusion product', repr(e))
        log('generated', len(ann_list), 'annotations', time_stamp=False)
    return results
