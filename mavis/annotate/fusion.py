from .genomic import Exon, Transcript, PreTranscript
from .protein import calculate_orf, Domain, Translation
from ..breakpoint import Breakpoint
from ..constants import ORIENT, PRIME, PROTOCOL, reverse_complement, STRAND, SVTYPE
from ..error import NotSpecifiedError
from ..interval import Interval, IntervalMapping


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


class FusionTranscript(PreTranscript):
    """
    FusionTranscript is a PreTranscript built from two parent PreTranscripts. It has most of the
    same functionality as a regular PreTranscript except that it will not have a parent gene and
    retains a mapping of the new exons to the exons in the PreTranscript they originated from

    Additionally the FusionTranscript is always constructed on the positive strand.

    The preferred way to construct a FusionTranscript is through the build method.
    """
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
        assert ann.event_type == SVTYPE.INV
        assert ann.transcript1 == ann.transcript2
        fusion_pre_transcript = FusionTranscript()
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
            fusion_pre_transcript.map_region_to_genome(b1.chr, (1, len(seq1)), (ann.transcript1.start, b1.start))
            fusion_pre_transcript.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + len(window_seq) + 1, len(seq1) + len(useq) + len(window_seq) + len(seq2)),
                (b2.start, ann.transcript1.end)
            )
        else:
            useq = reverse_complement(useq)  # exon sequence will already be revcomp from pull exons method
            seq1, seq2 = seq2, seq1
            ex1, ex2 = ex2, ex1
            fusion_pre_transcript.map_region_to_genome(b1.chr, (1, len(seq1)), (b2.start, ann.transcript1.end), True)
            fusion_pre_transcript.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + len(window_seq) + 1, len(seq1) + len(useq) + len(window_seq) + len(seq2)),
                (ann.transcript1.start, b1.start), True
            )

        fusion_pre_transcript.seq = seq1
        fusion_pre_transcript.break1 = len(seq1)
        fusion_pre_transcript.break2 = len(seq1) + len(useq) + len(window_seq) + 1

        if ann.break1.orient == ORIENT.LEFT:
            fusion_pre_transcript.seq += useq + window_seq
            fusion_pre_transcript.map_region_to_genome(
                b1.chr,
                (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(window_seq)),
                window,
                True if ann.transcript1.get_strand() == STRAND.POS else False
            )
        else:
            fusion_pre_transcript.seq += window_seq + useq
            fusion_pre_transcript.map_region_to_genome(
                b1.chr,
                (len(seq1) + 1, len(seq1) + len(window_seq)),
                window,
                True if ann.transcript1.get_strand() == STRAND.POS else False
            )

        fusion_pre_transcript.last_five_prime_exon = ex1[-1][0]
        fusion_pre_transcript.first_three_prime_exon = ex2[0][0]

        for ex, old_ex in ex1:
            fusion_pre_transcript.exons.append(ex)
            ex.reference_object = fusion_pre_transcript
            fusion_pre_transcript.exon_mapping[ex.position] = old_ex

        offset = len(fusion_pre_transcript.seq)

        if ann.protocol == PROTOCOL.TRANS:
            fusion_pre_transcript.exons[-1].end_splice_site.intact = False
            ex2[0][0].start_splice_site.intact = False
            novel_exon_start = fusion_pre_transcript.exons[-1].end + 1
            novel_exon_end = offset + ex2[0][0].start - 1
            if novel_exon_end >= novel_exon_start:  # create a novel exon
                exon = Exon(
                    novel_exon_start, novel_exon_end, fusion_pre_transcript,
                    intact_start_splice=False,
                    intact_end_splice=False,
                    seq=fusion_pre_transcript.seq[novel_exon_start - 1:novel_exon_end]
                )
                fusion_pre_transcript.exons.append(exon)
        fusion_pre_transcript.seq += seq2
        for ex, old_ex in ex2:
            exon = Exon(
                ex.start + offset, ex.end + offset, fusion_pre_transcript,
                intact_start_splice=ex.start_splice_site.intact,
                intact_end_splice=ex.end_splice_site.intact,
                seq=ex.seq
            )
            fusion_pre_transcript.exons.append(exon)
            fusion_pre_transcript.exon_mapping[exon.position] = old_ex
        return fusion_pre_transcript

    @classmethod
    def _build_single_gene_duplication(cls, ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match):
        """
        builds a fusion transcript for a single gene duplication. Note that this is an incomplete
        fusion transcript and still requires translations and domain information to be added
        """
        assert ann.event_type == SVTYPE.DUP
        assert ann.transcript1 == ann.transcript2
        fusion_pre_transcript = FusionTranscript()

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
        fusion_pre_transcript.map_region_to_genome(ann.break1.chr, Interval(1, len(seq2)), front, flipped)
        fusion_pre_transcript.map_region_to_genome(
            ann.break1.chr,
            Interval(len(seq2) + len(useq) + 1, len(seq1) + len(seq2) + len(useq)),
            back, flipped
        )
        fusion_pre_transcript.seq = seq2 + useq
        fusion_pre_transcript.break1 = len(seq2)
        fusion_pre_transcript.break2 = len(seq2) + len(useq) + 1
        for ex, old_ex in ex2:
            fusion_pre_transcript.exons.append(ex)
            ex.reference_object = fusion_pre_transcript
            fusion_pre_transcript.exon_mapping[ex.position] = old_ex
        fusion_pre_transcript.last_five_prime_exon = ex2[-1][0]
        fusion_pre_transcript.first_three_prime_exon = ex1[0][0]
        offset = len(fusion_pre_transcript.seq)
        if ann.protocol == PROTOCOL.TRANS:
            fusion_pre_transcript.exons[-1].end_splice_site.intact = False
            ex1[0][0].start_splice_site.intact = False
            novel_exon_start = fusion_pre_transcript.exons[-1].end + 1
            novel_exon_end = offset + ex1[0][0].start - 1
            if novel_exon_end >= novel_exon_start:  # create a novel exon
                e = Exon(
                    novel_exon_start, novel_exon_end, fusion_pre_transcript,
                    intact_start_splice=False,
                    intact_end_splice=False,
                    seq=fusion_pre_transcript.seq[novel_exon_start - 1:novel_exon_end]
                )
                fusion_pre_transcript.exons.append(e)

        for ex, old_ex in ex1:
            e = Exon(
                ex.start + offset, ex.end + offset, fusion_pre_transcript,
                intact_start_splice=ex.start_splice_site.intact,
                intact_end_splice=ex.end_splice_site.intact,
                seq=ex.seq
            )
            fusion_pre_transcript.exons.append(e)
            fusion_pre_transcript.exon_mapping[e.position] = old_ex
        fusion_pre_transcript.seq += seq1
        return fusion_pre_transcript

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

        fusion_pre_transcript = FusionTranscript()

        if ann.transcript1 == ann.transcript2 and ann.event_type not in [SVTYPE.DEL, SVTYPE.INS]:
            # single transcript events are special cases if the breakpoints face each other
            # as is the case for duplications and inversions
            if ann.event_type == SVTYPE.DUP:
                fusion_pre_transcript = cls._build_single_gene_duplication(
                    ann, reference_genome, min_orf_size, max_orf_cap, min_domain_mapping_match)
            elif ann.event_type == SVTYPE.INV:
                fusion_pre_transcript = cls._build_single_gene_inversion(
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
                fusion_pre_transcript.break1 = len(seq1)
                fusion_pre_transcript.break2 = len(seq1) + len(useq) + 1
                if ann.transcript1.strand == STRAND.NEG:
                    useq = reverse_complement(useq)
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break1.chr, (1, len(seq1)), (ann.break1.start, ann.transcript1.end), True
                    )
                else:
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break1.chr, (1, len(seq1)), (ann.transcript1.start, ann.break1.start)
                    )
                if ann.transcript2.strand == STRAND.NEG:  # t2 is the 3' transcript
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break2.chr,
                        (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.transcript2.start, ann.break2.start), True
                    )
                else:
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break2.chr,
                        (len(seq1) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.break2.start, ann.transcript2.end)
                    )
            else:
                fusion_pre_transcript.break1 = len(seq2)
                fusion_pre_transcript.break2 = len(seq2) + len(useq) + 1
                if ann.transcript2.strand == STRAND.NEG:
                    useq = reverse_complement(useq)
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break2.chr, (1, len(seq2)), (ann.break2.start, ann.transcript2.end), True
                    )
                else:
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break2.chr, (1, len(seq2)), (ann.transcript2.start, ann.break2.start)
                    )
                if ann.transcript1.strand == STRAND.NEG:  # t1 is the 3' transcript
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break1.chr,
                        (len(seq2) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.transcript1.start, ann.break1.start)
                    )
                else:
                    fusion_pre_transcript.map_region_to_genome(
                        ann.break1.chr,
                        (len(seq2) + len(useq) + 1, len(seq1) + len(useq) + len(seq2)),
                        (ann.break1.start, ann.transcript1.end)
                    )
                seq1, seq2 = seq2, seq1
                ex1, ex2 = ex2, ex1

            fusion_pre_transcript.seq = seq1 + useq

            for ex, old_ex in ex1:
                fusion_pre_transcript.exons.append(ex)
                ex.reference_object = fusion_pre_transcript
                fusion_pre_transcript.exon_mapping[ex.position] = old_ex
            offset = len(fusion_pre_transcript.seq)

            if ann.protocol == PROTOCOL.TRANS:
                fusion_pre_transcript.exons[-1].end_splice_site.intact = False
                ex2[0][0].start_splice_site.intact = False
                novel_exon_start = fusion_pre_transcript.exons[-1].end + 1
                novel_exon_end = offset + ex2[0][0].start - 1
                if novel_exon_end >= novel_exon_start:  # create a novel exon
                    e = Exon(
                        novel_exon_start, novel_exon_end, fusion_pre_transcript,
                        intact_start_splice=False,
                        intact_end_splice=False,
                        seq=fusion_pre_transcript.seq[novel_exon_start - 1:novel_exon_end]
                    )
                    fusion_pre_transcript.exons.append(e)
            for ex, old_ex in ex2:
                e = Exon(
                    ex.start + offset, ex.end + offset, fusion_pre_transcript,
                    intact_start_splice=ex.start_splice_site.intact,
                    intact_end_splice=ex.end_splice_site.intact,
                    seq=ex.seq
                )
                fusion_pre_transcript.exons.append(e)
                fusion_pre_transcript.exon_mapping[e.position] = old_ex
            fusion_pre_transcript.seq += seq2

        fusion_pre_transcript.position = Interval(1, len(fusion_pre_transcript.seq))
        # add all splice variants
        for spl_patt in fusion_pre_transcript.generate_splicing_patterns():
            fusion_spl_tx = Transcript(fusion_pre_transcript, spl_patt)
            fusion_pre_transcript.spliced_transcripts.append(fusion_spl_tx)

            # calculate the putative open reading frames
            orfs = calculate_orf(fusion_spl_tx.get_seq(), min_orf_size=min_orf_size)
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
                translation = Translation(
                    orf.start - fusion_spl_tx.start + 1,
                    orf.end - fusion_spl_tx.start + 1, fusion_spl_tx)
                fusion_spl_tx.translations.append(translation)

        # remap the domains from the original translations to the current translations
        for fusion_spl_tx in fusion_pre_transcript.spliced_transcripts:
            for new_tl in fusion_spl_tx.translations:
                aa_seq = new_tl.get_aa_seq()
                assert aa_seq[0] == 'M'
                translations = ann.transcript1.translations[:]
                if ann.transcript1 != ann.transcript2:
                    translations += ann.transcript2.translations
                for translation in translations:
                    for dom in translation.domains:
                        try:
                            match, total, regions = dom.align_seq(aa_seq, reference_genome)
                            if min_domain_mapping_match is None or match / total >= min_domain_mapping_match:
                                new_dom = Domain(dom.name, regions, new_tl)
                                new_tl.domains.append(new_dom)
                        except UserWarning:
                            pass
        return fusion_pre_transcript

    def get_seq(self, reference_genome=None, ignore_cache=False):
        return PreTranscript.get_seq(self)

    def get_cdna_seq(self, splicing_pattern, reference_genome=None, ignore_cache=False):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): the list of splicing positions
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference seq
                by template/chr name

        Returns:
            str: the spliced cDNA seq
        """
        return PreTranscript.get_cdna_seq(self, splicing_pattern)

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
                    e.seq = str(temp)
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
                    e.seq = str(temp)
                    assert len(temp) == len(e)
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
                    e.seq = str(temp)
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
