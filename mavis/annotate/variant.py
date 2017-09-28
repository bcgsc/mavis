import itertools
import json
import uuid

from .genomic import Exon, IntergenicRegion, Transcript, UsTranscript
from .protein import calculate_orf, Domain, Translation
from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import COLUMNS, GENE_PRODUCT_TYPE, ORIENT, PRIME, PROTOCOL, reverse_complement, STRAND, SVTYPE
from ..error import NotSpecifiedError
from ..interval import Interval, IntervalMapping
from ..util import devnull


def determine_prime(transcript, breakpoint):
    """
    determine the side of the transcript 5' or 3' which is 'kept' given the breakpoint

    Args:
        transcript (Transcript): the transcript
        breakpoint (Breakpoint): the breakpoint

    Returns:
        PRIME: 5' or 3'

    Raises:
        AttributeError: if the orientation of the breakpoint or the strand of the transcript is not specified
    """
    if transcript.get_strand() == STRAND.POS:
        if breakpoint.orient == ORIENT.LEFT:
            return PRIME.FIVE
        elif breakpoint.orient == ORIENT.RIGHT:
            return PRIME.THREE
        else:
            raise NotSpecifiedError('cannot determine_prime if the orient of the breakpoint has not been specified')
    elif transcript.get_strand() == STRAND.NEG:
        if breakpoint.orient == ORIENT.LEFT:
            return PRIME.THREE
        elif breakpoint.orient == ORIENT.RIGHT:
            return PRIME.FIVE
        else:
            raise NotSpecifiedError('cannot determine_prime if the orient of the breakpoint has not been specified')
    else:
        raise NotSpecifiedError('cannot determine prime if the strand of the transcript has not been specified')


class FusionTranscript(UsTranscript):

    def __init__(self):
        self.exon_mapping = {}
        self.exons = []
        self.seq = None
        self.spliced_transcripts = []
        self.position = None
        self.strand = STRAND.POS  # always built on the positive strand
        self.reference_object = None
        self.name = None
        # interval to interval mapping of transcript coord to their original genome coord
        self.mapping_to_genome = IntervalMapping()
        self.mapping_to_chrs = dict()  # keeps track of what chromosome per interval
        self.break1 = None  # first breakpoint position in the fusion transcript
        self.break2 = None  # second breakpoint position in the fusion transcript

    def exon_number(self, exon):
        """
        Args:
            exon (Exon): the exon to be numbered

        Returns:
            int: the number of the exon in the original transcript (prior to fusion)
        """
        old_exon = self.exon_mapping[exon.position]
        return old_exon.transcript.exon_number(old_exon)

    def map_region_to_genome(self, chr, interval_on_fusion, genome_interval, flipped=False):
        self.mapping_to_genome.add(interval_on_fusion, genome_interval, flipped)
        self.mapping_to_chrs[Interval(interval_on_fusion[0], interval_on_fusion[1])] = chr

    @classmethod
    def _build_single_gene_inversion(cls, ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match):
        """
        builds a fusion transcript for a single gene inversion. Note that this is an incomplete
        fusion transcript and still requires translations and domain information to be added
        """
        assert(ann.event_type == SVTYPE.INV)
        assert(ann.transcript1 == ann.transcript2)
        ft = FusionTranscript()
        # pull the exons from either size of the breakpoints window
        if ann.break1.orient == ORIENT.RIGHT:
            window = Interval(ann.break1.end, ann.break2.end - 1)
        else:
            window = Interval(ann.break1.end + 1, ann.break2.end)
        window_seq = reference_genome[ann.break1.chr].seq[window.start - 1:window.end]
        # now create 'pseudo-deletion' breakpoints
        b1 = Breakpoint(ann.break1.chr, window.start - 1, orient=ORIENT.LEFT)
        b2 = Breakpoint(ann.break2.chr, window.end + 1, orient=ORIENT.RIGHT)

        seq1, ex1 = cls._pull_exons(ann.transcript1, b1, reference_genome[b1.chr].seq)
        seq2, ex2 = cls._pull_exons(ann.transcript2, b2, reference_genome[b2.chr].seq)
        useq = ann.untemplated_seq

        if ann.transcript1.get_strand() == STRAND.POS:
            window_seq = reverse_complement(window_seq)  # b/c inversion should be opposite
            ft.map_region_to_genome(b1.chr, (1, len(seq1)), (ann.transcript1.start, b1.start))
            ft.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + len(window_seq) + 1, len(seq1) + len(useq) + len(window_seq) + len(seq2)),
                (b2.start, ann.transcript1.end)
            )
        else:
            useq = reverse_complement(useq)  # exon sequence will already be revcomp from pull exons method
            seq1, seq2 = seq2, seq1
            ex1, ex2 = ex2, ex1
            ft.map_region_to_genome(b1.chr, (1, len(seq1)), (b2.start, ann.transcript1.end), True)
            ft.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + len(window_seq) + 1, len(seq1) + len(useq) + len(window_seq) + len(seq2)),
                (ann.transcript1.start, b1.start), True
            )

        ft.seq = seq1
        ft.break1 = len(seq1)
        ft.break2 = len(seq1) + len(useq) + len(window_seq) + 1

        if ann.break1.orient == ORIENT.LEFT:
            ft.seq += useq + window_seq
            ft.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(window_seq)),
                window,
                True if ann.transcript1.get_strand() == STRAND.POS else False
            )
        else:
            ft.seq += window_seq + useq
            ft.map_region_to_genome(
                b1.chr,
                (len(seq1) + 1, len(seq1) + len(window_seq)),
                window,
                True if ann.transcript1.get_strand() == STRAND.POS else False
            )

        ft.last_five_prime_exon = ex1[-1][0]
        ft.first_three_prime_exon = ex2[0][0]

        for ex, old_ex in ex1:
            ft.exons.append(ex)
            ex.reference_object = ft
            ft.exon_mapping[ex.position] = old_ex

        offset = len(ft.seq)

        if ann.protocol == PROTOCOL.TRANS:
            ft.exons[-1].end_splice_site.intact = False
            ex2[0][0].start_splice_site.intact = False
            novel_exon_start = ft.exons[-1].end + 1
            novel_exon_end = offset + ex2[0][0].start - 1
            if novel_exon_end >= novel_exon_start:  # create a novel exon
                e = Exon(
                    novel_exon_start, novel_exon_end, ft,
                    intact_start_splice=False,
                    intact_end_splice=False
                )
                ft.exons.append(e)
        ft.seq += seq2
        for ex, old_ex in ex2:
            e = Exon(
                ex.start + offset, ex.end + offset, ft,
                intact_start_splice=ex.start_splice_site.intact,
                intact_end_splice=ex.end_splice_site.intact
            )
            ft.exons.append(e)
            ft.exon_mapping[e.position] = old_ex
        return ft

    @classmethod
    def _build_single_gene_duplication(cls, ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match):
        """
        builds a fusion transcript for a single gene duplication. Note that this is an incomplete
        fusion transcript and still requires translations and domain information to be added
        """
        assert(ann.event_type == SVTYPE.DUP)
        assert(ann.transcript1 == ann.transcript2)
        ft = FusionTranscript()

        seq1, ex1 = cls._pull_exons(ann.transcript1, ann.break1, reference_genome[ann.break1.chr].seq)
        seq2, ex2 = cls._pull_exons(ann.transcript2, ann.break2, reference_genome[ann.break2.chr].seq)
        useq = ann.untemplated_seq
        front = Interval(ann.transcript1.start, ann.break2.start)
        back = Interval(ann.break1.start, ann.transcript1.end)
        flipped = False
        if ann.transcript1.get_strand() == STRAND.NEG:
            seq1, seq2 = (seq2, seq1)
            ex1, ex2 = (ex2, ex1)
            front, back = (back, front)
            useq = reverse_complement(useq)
            flipped = True
        ft.map_region_to_genome(ann.break1.chr, Interval(1, len(seq2)), front, flipped)
        ft.map_region_to_genome(
            ann.break1.chr,
            Interval(len(seq2) + len(useq) + 1, len(seq1) + len(seq2) + len(useq)),
            back, flipped
        )
        ft.seq = seq2 + useq
        ft.break1 = len(seq2)
        ft.break2 = len(seq2) + len(useq) + 1
        for ex, old_ex in ex2:
            ft.exons.append(ex)
            ex.reference_object = ft
            ft.exon_mapping[ex.position] = old_ex
        ft.last_five_prime_exon = ex2[-1][0]
        ft.first_three_prime_exon = ex1[0][0]
        offset = len(ft.seq)
        if ann.protocol == PROTOCOL.TRANS:
            ft.exons[-1].end_splice_site.intact = False
            ex1[0][0].start_splice_site.intact = False
            novel_exon_start = ft.exons[-1].end + 1
            novel_exon_end = offset + ex1[0][0].start - 1
            if novel_exon_end >= novel_exon_start:  # create a novel exon
                e = Exon(
                    novel_exon_start, novel_exon_end, ft,
                    intact_start_splice=False,
                    intact_end_splice=False
                )
                ft.exons.append(e)

        for ex, old_ex in ex1:
            e = Exon(
                ex.start + offset, ex.end + offset, ft,
                intact_start_splice=ex.start_splice_site.intact,
                intact_end_splice=ex.end_splice_site.intact
            )
            ft.exons.append(e)
            ft.exon_mapping[e.position] = old_ex
        ft.seq += seq1
        return ft

    @classmethod
    def build(cls, ann, reference_genome, min_orf_size=None, max_orf_cap=None, min_domain_mapping_match=None):
        """
        Args:
            ann (Annotation): the annotation object we want to build a FusionTranscript for
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            FusionTranscript: the newly built fusion transcript

        """
        if not ann.transcript1 or not ann.transcript2:
            raise NotSpecifiedError('cannot produce fusion transcript for non-annotated fusions')
        elif not ann.event_type and ann.transcript1 == ann.transcript2:
            raise NotSpecifiedError('event_type must be specified to produce a fusion transcript')
        elif ann.untemplated_seq is None:
            raise NotSpecifiedError(
                'cannot build a fusion transcript where the untemplated sequence has not been specified')
        elif len(ann.break1) > 1 or len(ann.break2) > 1:
            raise NotSpecifiedError('cannot build fusion transcripts for non-specific breakpoints')

        ft = FusionTranscript()

        if ann.transcript1 == ann.transcript2 and ann.event_type not in [SVTYPE.DEL, SVTYPE.INS]:
            # single transcript events are special cases if the breakpoints face each other
            # as is the case for duplications and inversions
            if ann.event_type == SVTYPE.DUP:
                ft = cls._build_single_gene_duplication(
                    ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match)
            elif ann.event_type == SVTYPE.INV:
                ft = cls._build_single_gene_inversion(
                    ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match)
            else:
                raise AttributeError('unrecognized event type')
        else:
            t1 = determine_prime(ann.transcript1, ann.break1)
            t2 = determine_prime(ann.transcript2, ann.break2)

            if t1 == t2:
                raise NotImplementedError('do not produce fusion transcript for anti-sense fusions')
            seq1, ex1 = cls._pull_exons(ann.transcript1, ann.break1, reference_genome[ann.break1.chr].seq)
            seq2, ex2 = cls._pull_exons(ann.transcript2, ann.break2, reference_genome[ann.break2.chr].seq)
            useq = ann.untemplated_seq
            if t1 == PRIME.FIVE:
                ft.break1 = len(seq1)
                ft.break2 = len(seq1) + len(useq) + 1
                if ann.transcript1.strand == STRAND.NEG:
                    useq = reverse_complement(useq)
                    ft.map_region_to_genome(
                        ann.break1.chr, (1, len(seq1)), (ann.break1.start, ann.transcript1.end), True
                    )
                else:
                    ft.map_region_to_genome(
                        ann.break1.chr, (1, len(seq1)), (ann.transcript1.start, ann.break1.start)
                    )
                if ann.transcript2.strand == STRAND.NEG:  # t2 is the 3' transcript
                    ft.map_region_to_genome(
                        ann.break2.chr,
                        (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.transcript2.start, ann.break2.start), True
                    )
                else:
                    ft.map_region_to_genome(
                        ann.break2.chr,
                        (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.break2.start, ann.transcript2.end)
                    )
            else:
                ft.break1 = len(seq2)
                ft.break2 = len(seq2) + len(useq) + 1
                if ann.transcript2.strand == STRAND.NEG:
                    useq = reverse_complement(useq)
                    ft.map_region_to_genome(
                        ann.break2.chr, (1, len(seq2)), (ann.break2.start, ann.transcript2.end), True
                    )
                else:
                    ft.map_region_to_genome(
                        ann.break2.chr, (1, len(seq2)), (ann.transcript2.start, ann.break2.start)
                    )
                if ann.transcript1.strand == STRAND.NEG:  # t1 is the 3' transcript
                    ft.map_region_to_genome(
                        ann.break1.chr,
                        (len(seq2) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.transcript1.start, ann.break1.start)
                    )
                else:
                    ft.map_region_to_genome(
                        ann.break1.chr,
                        (len(seq2) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.break1.start, ann.transcript1.end)
                    )
                seq1, seq2 = seq2, seq1
                ex1, ex2 = ex2, ex1

            ft.seq = seq1 + useq

            for ex, old_ex in ex1:
                ft.exons.append(ex)
                ex.reference_object = ft
                ft.exon_mapping[ex.position] = old_ex
            offset = len(ft.seq)

            if ann.protocol == PROTOCOL.TRANS:
                ft.exons[-1].end_splice_site.intact = False
                ex2[0][0].start_splice_site.intact = False
                novel_exon_start = ft.exons[-1].end + 1
                novel_exon_end = offset + ex2[0][0].start - 1
                if novel_exon_end >= novel_exon_start:  # create a novel exon
                    e = Exon(
                        novel_exon_start, novel_exon_end, ft,
                        intact_start_splice=False,
                        intact_end_splice=False
                    )
                    ft.exons.append(e)
            for ex, old_ex in ex2:
                e = Exon(
                    ex.start + offset, ex.end + offset, ft,
                    intact_start_splice=ex.start_splice_site.intact,
                    intact_end_splice=ex.end_splice_site.intact
                )
                ft.exons.append(e)
                ft.exon_mapping[e.position] = old_ex
            ft.seq += seq2

        ft.position = Interval(1, len(ft.seq))
        # add all splice variants
        for spl_patt in ft.generate_splicing_patterns():
            t = Transcript(ft, spl_patt)
            ft.spliced_transcripts.append(t)

            # calculate the putative open reading frames
            orfs = calculate_orf(t.get_seq(), min_orf_size=min_orf_size)
            # limit the length to either only the longest ORF or anything longer than the input translations
            min_orf_length = max([len(o) for o in orfs] + [min_orf_size if min_orf_size else 0])
            for ref_tx in [ann.transcript1, ann.transcript2]:
                for tlx in ref_tx.translations:
                    min_orf_length = min(min_orf_length, len(tlx))

            # filter the orfs based on size
            orfs = [o for o in orfs if len(o) >= min_orf_length]

            # if there are still too many filter to reasonable number
            if max_orf_cap and len(orfs) > max_orf_cap:  # limit the number of orfs returned
                orfs = sorted(orfs, key=lambda x: len(x), reverse=True)
                orfs = orfs[0:max_orf_cap]
            # create the translations
            for orf in orfs:
                tl = Translation(orf.start - t.start + 1, orf.end - t.start + 1, t)
                t.translations.append(tl)

        # remap the domains from the original translations to the current translations
        for t in ft.spliced_transcripts:
            for new_tl in t.translations:
                aa_seq = new_tl.get_aa_seq()
                assert(aa_seq[0] == 'M')
                translations = ann.transcript1.translations[:]
                if ann.transcript1 != ann.transcript2:
                    translations += ann.transcript2.translations
                for tl in translations:
                    for dom in tl.domains:
                        try:
                            match, total, regions = dom.align_seq(aa_seq, reference_genome)
                            if min_domain_mapping_match is None or match / total >= min_domain_mapping_match:
                                new_dom = Domain(dom.name, regions, new_tl)
                                new_tl.domains.append(new_dom)
                        except UserWarning:
                            pass
        return ft

    def get_seq(self, reference_genome=None, ignore_cache=False):
        return UsTranscript.get_seq(self)

    def get_spliced_cdna_seq(self, splicing_pattern, reference_genome=None, ignore_cache=False):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): the list of splicing positions
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference seq
                by template/chr name

        Returns:
            str: the spliced cDNA seq
        """
        return UsTranscript.get_spliced_cdna_seq(self, splicing_pattern)

    @classmethod
    def _pull_exons(cls, transcript, breakpoint, reference_sequence):
        """
        given a transcript and breakpoint returns the exons and sequence expected
        the exons are returned in the order wrt to strand (i.e. reversed from a genomic sort
        if they are on the negative strand). The positions of the exons are wrt to the sequence
        being returned (starts at 1).
        """
        if len(breakpoint) > 1:
            raise AttributeError('cannot pull exons on non-specific breakpoints')
        new_exons = []
        s = ''
        exons = sorted(transcript.exons, key=lambda x: x.start)
        if breakpoint.orient == ORIENT.LEFT:  # five prime
            for i, exon in enumerate(exons):
                intact_start_splice = True
                intact_end_splice = True
                if breakpoint.start < exon.start:  # =====----|----|----
                    if i > 0:  # add intron
                        temp = reference_sequence[exons[i - 1].end:breakpoint.start]
                        s += temp
                    break
                else:
                    if i > 0:  # add intron
                        temp = reference_sequence[exons[i - 1].end:exon.start - 1]
                        s += temp
                    if breakpoint.start <= exon.end_splice_site.end:
                        intact_end_splice = False
                    if breakpoint.start <= exon.start_splice_site.end:
                        intact_start_splice = False
                    t = min(breakpoint.start, exon.end)
                    e = Exon(
                        len(s) + 1, len(s) + t - exon.start + 1,
                        intact_start_splice=intact_start_splice,
                        intact_end_splice=intact_end_splice,
                        strand=STRAND.POS
                    )
                    temp = reference_sequence[exon.start - 1:t]
                    s += temp
                    new_exons.append((e, exon))
        elif breakpoint.orient == ORIENT.RIGHT:  # three prime
            for i, exon in enumerate(exons):
                intact_start_splice = True
                intact_end_splice = True

                if breakpoint.start < exon.start:  # --==|====|====
                    if i > 0:  # add last intron
                        t = max(exons[i - 1].end + 1, breakpoint.start)
                        temp = reference_sequence[t - 1:exon.start - 1]
                        s += temp
                    if Interval.overlaps(breakpoint, exon.start_splice_site):
                        intact_start_splice = False
                    # add the exon
                    e = Exon(
                        len(s) + 1, len(s) + len(exon),
                        intact_start_splice=intact_start_splice,
                        intact_end_splice=intact_end_splice,
                        strand=STRAND.POS
                    )
                    temp = reference_sequence[exon.start - 1:exon.end]
                    assert(len(temp) == len(e))
                    s += temp
                    new_exons.append((e, exon))
                elif breakpoint.start <= exon.end:  # --|-====|====
                    intact_start_splice = False
                    temp = reference_sequence[breakpoint.start - 1:exon.end]
                    if Interval.overlaps(breakpoint, exon.end_splice_site):
                        intact_end_splice = False
                    # add the exon
                    e = Exon(
                        len(s) + 1, len(s) + len(temp),
                        intact_start_splice=intact_start_splice,
                        intact_end_splice=intact_end_splice,
                        strand=STRAND.POS
                    )
                    s += temp
                    new_exons.append((e, exon))
        else:
            raise NotSpecifiedError('breakpoint orientation must be specified to pull exons')
        if transcript.get_strand() == STRAND.NEG:
            # reverse complement the sequence and reverse the exons
            temp = new_exons
            new_exons = []

            for ex, old_exon in temp[::-1]:
                e = Exon(
                    len(s) - ex.end + 1,
                    len(s) - ex.start + 1,
                    intact_start_splice=ex.end_splice_site.intact,
                    intact_end_splice=ex.start_splice_site.intact,
                    strand=STRAND.POS
                )
                new_exons.append((e, old_exon))
            s = reverse_complement(s)
        elif transcript.get_strand() != STRAND.POS:
            raise NotSpecifiedError('transcript strand must be specified to pull exons')

        return str(s), new_exons


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
            if len(conflicts) > 0:
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

    def add_gene(self, gene):
        """
        adds a gene to the current set of annotations. Checks which set it should be added to

        Args:
            gene (Gene): the gene being added
        """
        if gene.chr not in [self.break1.chr, self.break2.chr]:
            raise AttributeError('cannot add gene not on the same chromosome as either breakpoint')

        if not self.interchromosomal:
            try:
                encompassment = Interval(self.break1.end + 1, self.break2.start - 1)
                if gene in encompassment:
                    self.encompassed_genes.add(gene)
            except AttributeError:
                pass
        if Interval.overlaps(gene, self.break1) and gene.chr == self.break1.chr \
                and gene != self.transcript1.reference_object:
            self.genes_overlapping_break1.add(gene)
        if Interval.overlaps(gene, self.break2) and gene.chr == self.break2.chr \
                and gene != self.transcript2.reference_object:
            self.genes_overlapping_break2.add(gene)

        if gene in self.genes_overlapping_break1 or gene in self.genes_overlapping_break2 or \
                gene in self.encompassed_genes or gene == self.transcript1.reference_object or \
                gene == self.transcript2.reference_object:
            return

        d1 = Interval.dist(gene, self.break1)
        d2 = Interval.dist(gene, self.break2)

        if self.interchromosomal:
            if gene.chr == self.break1.chr:
                self.genes_proximal_to_break1.add((gene, d1))
            elif gene.chr == self.break2.chr:
                self.genes_proximal_to_break2.add((gene, d2))
        else:
            if d1 < 0:
                self.genes_proximal_to_break1.add((gene, d1))
            if d2 > 0:
                self.genes_proximal_to_break2.add((gene, d2))

        if len(self.genes_proximal_to_break1) > 0:
            temp = set()
            tgt = min([abs(d) for g, d in self.genes_proximal_to_break1])

            for gene, dist in self.genes_proximal_to_break1:
                if self.proximity is None:
                    if abs(dist) == tgt:
                        temp.add((gene, dist))
                elif abs(dist) <= self.proximity:
                    temp.add((gene, dist))

            self.genes_proximal_to_break1 = temp

        if len(self.genes_proximal_to_break2) > 0:
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
        :class:`list` of :any:`UsTranscript`: a list of possible transcripts
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

            - :class:`list` of (:class:`UsTranscript` or :class:`IntergenicRegion`): transcripts or intergenic regions
              overlapping the breakpoint on the positive strand
            - :class:`list` of (:class:`UsTranscript` or :class:`IntergenicRegion`): transcripts or intergenic regions
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
    log=devnull,
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
        bpp.data[COLUMNS.validation_id] = bpp.data.get(COLUMNS.validation_id, str(uuid.uuid4()))
        ann = _gather_annotations(
            annotations,
            bpp,
            proximity=max_proximity
        )
        for f in filters:
            ann = f(ann)  # apply the filter
        results.extend(ann)
        for j, a in enumerate(ann):
            a.data[COLUMNS.annotation_id] = '{}-a{}'.format(a.validation_id, j + 1)
            # try building the fusion product
            try:
                ft = FusionTranscript.build(
                    a, reference_genome,
                    min_orf_size=min_orf_size,
                    max_orf_cap=max_orf_cap,
                    min_domain_mapping_match=min_domain_mapping_match
                )
                a.fusion = ft
            except (NotSpecifiedError, AttributeError, NotImplementedError):
                pass
            except KeyError as e:
                log('warning. could not build fusion product', repr(e))
        log('generated', len(ann), 'annotations', time_stamp=False)
    return results
