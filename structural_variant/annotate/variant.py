from .genomic import Transcript, Exon, IntergenicRegion
from ..constants import STRAND, SVTYPE, reverse_complement, ORIENT, PRIME, COLUMNS
from ..breakpoint import Breakpoint, BreakpointPair
from ..interval import Interval
from .protein import Translation, Domain, calculate_ORF
import itertools


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
    if transcript.strand == STRAND.POS:
        if breakpoint.orient == ORIENT.LEFT:
            return PRIME.FIVE
        elif breakpoint.orient == ORIENT.RIGHT:
            return PRIME.THREE
        else:
            raise AttributeError('cannot determine_prime if the orient of the breakpoint has not been specified')
    elif transcript.strand == STRAND.NEG:
        if breakpoint.orient == ORIENT.LEFT:
            return PRIME.THREE
        elif breakpoint.orient == ORIENT.RIGHT:
            return PRIME.FIVE
        else:
            raise AttributeError('cannot determine_prime if the orient of the breakpoint has not been specified')
    else:
        raise AttributeError('cannot determine prime if the strand of the transcript has not been specified')


class FusionTranscript(Transcript):

    def __init__(self):
        self.exon_mapping = {}
        self.exons = []
        self.sequence = ''
        self.translations = []
        self.domains = []
        self.position = None
        self._strand = STRAND.POS  # always built on the positive strand

    def exon_number(self, exon):
        """
        Args:
            exon (Exon): the exon to be numbered

        Returns:
            int: the number of the exon in the original transcript (prior to fusion)
        """
        old_exon = self.exon_mapping[exon]
        return old_exon.transcript.exon_number(old_exon)

    @classmethod
    def build(cls, ann, REFERENCE_GENOME, min_orf_size=None, max_orf_cap=None, min_domain_mapping_match=None):
        """
        Args:
            ann (Annotation): the annotation object we want to build a FusionTranscript for
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

        Returns:
            FusionTranscript: the newly built fusion transcript

        .. todo::
            support single transcript inversions

        """
        if not ann.transcript1 or not ann.transcript2:
            raise AttributeError('cannot produce fusion transcript for non-annotated fusions')
        elif not ann.event_type and ann.transcript1 == ann.transcript2:
            raise AttributeError('event_type must be specified to produce a fusion transcript')
        elif ann.untemplated_sequence is None:
            raise AttributeError(
                'cannot build a fusion transcript where the untemplated sequence has not been specified')
        
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('BUILDING FUSION:', ann.transcript1.name, ann.transcript2.name, ann)
        print('transcript1 splicing', ann.transcript1.name, ann.transcript1.splicing_patterns())
        print('transcript2 splicing', ann.transcript2.name, ann.transcript2.splicing_patterns())
        ft = FusionTranscript()

        if ann.transcript1 == ann.transcript2 and ann.event_type not in [SVTYPE.DEL, SVTYPE.INS]:
            # single transcript events are special cases if the breakpoints face each other
            # as is the case for duplications and inversions
            if ann.event_type == SVTYPE.DUP:
                seq1, ex1 = cls._pull_exons(ann.transcript1, ann.break1, REFERENCE_GENOME[ann.break1.chr].seq)
                seq2, ex2 = cls._pull_exons(ann.transcript2, ann.break2, REFERENCE_GENOME[ann.break2.chr].seq)
                useq = ann.untemplated_sequence

                if ann.transcript1.strand == STRAND.NEG:
                    seq1, seq2 = (seq2, seq1)
                    ex1, ex2 = (ex2, ex1)
                    useq = reverse_complement(useq)
                ft.sequence = seq2 + useq
                for ex, old_ex in ex2:
                    ft.exons.append(ex)
                    ft.exon_mapping[ex] = old_ex
                offset = len(ft.sequence) + 1
                for ex, old_ex in ex1:
                    e = Exon(
                        ex.start + offset, ex.end + offset,
                        intact_start_splice=ex.intact_start_splice,
                        intact_end_splice=ex.intact_end_splice
                    )
                    ft.exons.append(e)
                    ft.exon_mapping[e] = old_ex
                ft.sequence += seq1
            elif ann.event_type == SVTYPE.INV:
                # pull the exons from either size of the breakpoints window
                window = Interval(ann.break1.end + 1, ann.break2.end)
                if ann.break1.orient == ORIENT.RIGHT:
                    window = Interval(ann.break1.end + 1, ann.break2.start - 1)
                window_seq = REFERENCE_GENOME[ann.break1.chr].seq[window.start - 1:window.end]
                # now create 'pseudo-deletion' breakpoints
                b1 = Breakpoint(ann.break1.chr, window.start - 1, orient=ORIENT.LEFT)
                b2 = Breakpoint(ann.break2.chr, window.end + 1, orient=ORIENT.RIGHT)

                seq1, ex1 = cls._pull_exons(ann.transcript1, b1, REFERENCE_GENOME[b1.chr].seq)
                seq2, ex2 = cls._pull_exons(ann.transcript2, b1, REFERENCE_GENOME[b2.chr].seq)
                useq = ann.untemplated_sequence

                if ann.transcript1.strand == STRAND.POS:
                    ft.sequence = seq1 + useq + window_seq
                    for ex, old_ex in ex1:
                        ft.exons.append(ex)
                        ft.exon_mapping[ex] = old_ex
                    offset = len(ft.sequence)
                    for ex, old_ex in ex2:
                        e = Exon(
                            ex.start + offset, ex.end + offset,
                            intact_start_splice=ex.intact_start_splice,
                            intact_end_splice=ex.intact_end_splice
                        )
                        ft.exons.append(e)
                        ft.exon_mapping[e] = old_ex
                    ft.sequence += seq2
                else:
                    pass
                raise NotImplementedError('single transcript inversion have not yet been implemented')
            else:
                raise AttributeError('unrecognized event type')
        else:
            t1 = determine_prime(ann.transcript1, ann.break1)
            t2 = determine_prime(ann.transcript2, ann.break2)

            if t1 == t2:
                raise NotImplementedError('do not produce fusion transcript for anti-sense fusions')
            seq1, ex1 = cls._pull_exons(ann.transcript1, ann.break1, REFERENCE_GENOME[ann.break1.chr].seq)
            seq2, ex2 = cls._pull_exons(ann.transcript2, ann.break2, REFERENCE_GENOME[ann.break2.chr].seq)

            if t1 == PRIME.FIVE:
                ft.sequence = seq1 + ann.untemplated_sequence if ann.transcript1.strand == STRAND.POS else \
                    reverse_complement(ann.untemplated_sequence)
                for ex, old_ex in ex1:
                    ft.exons.append(ex)
                    ft.exon_mapping[ex] = old_ex
                offset = len(ft.sequence)
                for ex, old_ex in ex2:
                    e = Exon(
                        ex.start + offset, ex.end + offset,
                        intact_start_splice=ex.intact_start_splice,
                        intact_end_splice=ex.intact_end_splice
                    )
                    ft.exons.append(e)
                    ft.exon_mapping[e] = old_ex
                ft.sequence += seq2
            else:
                ft.sequence = seq2 + ann.untemplated_sequence if ann.transcript2.strand == STRAND.POS else \
                    reverse_complement(ann.untemplated_sequence)
                for ex, old_ex in ex2:
                    ft.exons.append(ex)
                    ft.exon_mapping[ex] = old_ex
                offset = len(ft.sequence)
                for ex, old_ex in ex1:
                    e = Exon(
                        ex.start + offset, ex.end + offset,
                        intact_start_splice=ex.intact_start_splice,
                        intact_end_splice=ex.intact_end_splice
                    )
                    ft.exons.append(e)
                    ft.exon_mapping[e] = old_ex
                ft.sequence += seq1
        ft.position = Interval(1, len(ft.sequence))
        # get all splicing patterns
        for spl_patt in ft.splicing_patterns():
            # calculate all ORF for each splicing pattern
            spliced_cdna_seq = ft.get_spliced_cdna_sequence(spl_patt)
            orfs = calculate_ORF(spliced_cdna_seq, min_orf_size=min_orf_size)
            if max_orf_cap and len(orfs) > max_orf_cap:  # limit the number of orfs returned
                orfs = sorted(orfs, key=lambda x: len(x), reverse=True)
                l = len(orfs[max_orf_cap - 1])
                temp = []
                for i, orf in enumerate(orfs):
                    if len(orf) < l:
                        break
                    else:
                        temp.append(orf)
                orfs = temp
            # create the translations
            for orf in orfs:
                s = spliced_cdna_seq[orf.start - 1:orf.end]
                tl = Translation(orf.start, orf.end, ft, spl_patt, sequence=s)

        # remap the domains from the original translations to the current translations
        for new_tl in ft.translations:
            aa_seq = new_tl.get_AA_sequence()
            print('ft translation', new_tl.get_AA_sequence())
            print('ft splicing', new_tl.splicing_pattern)
            for tl in ann.transcript1.translations + ann.transcript2.translations:
                for dom in tl.domains:
                    try:
                        match, total, regions = dom.align_seq(aa_seq, REFERENCE_GENOME)
                        if min_domain_mapping_match is None or match / total >= min_domain_mapping_match:
                            Domain(dom.name, regions, new_tl)
                        else:
                            print('domain failed to map', ann.transcript1.name, 'fused to', ann.transcript2.name)
                            print('domain', dom.name, 'match=', match, 'total=', total, 'min_domain_mapping_match=', min_domain_mapping_match)
                            print('domain original transcript', dom.translation.transcript.name)
                            print('domain coord', [(dr.start, dr.end) for dr in dom.regions], dom.name)
                            print('tranlsation aa seq')
                            print(aa_seq)
                            print('domain region seqs', dom.get_sequences(REFERENCE_GENOME))
                    except UserWarning:
                        pass
        return ft
    
    def get_sequence(self):
        if self.sequence:
            return self.sequence
        elif self.gene and self.gene.sequence:
            # gene has a sequence set
            start = self.start - self.gene.start
            end = self.end - self.gene.end + len(self.gene.sequence)
            if self.strand == STRAND.NEG:
                return reverse_complement(self.gene.sequence[start:end])
            else:
                return self.gene.sequence[start:end]
        else:
            raise AttributeError('insufficient sequence information')
    
    def get_spliced_cdna_sequence(self, splicing_pattern):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): the list of splicing positions
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

        Returns:
            str: the spliced cDNA sequence
        """
        temp = sorted([self.start] + splicing_pattern + [self.end])
        m = min(temp)
        conti = []
        for i in range(0, len(temp) - 1, 2):
            conti.append(Interval(temp[i] - m, temp[i + 1] - m))
        seq = self.get_sequence()
        if self.strand == STRAND.NEG:
            # adjust the continuous intervals for the min and flip if revcomp
            seq = reverse_complement(seq)
        spliced_seq = ''.join([str(seq[i.start:i.end + 1]) for i in conti])
        return spliced_seq if self.strand == STRAND.POS else reverse_complement(spliced_seq)

    @classmethod
    def _pull_exons(cls, transcript, breakpoint, reference_sequence):
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
                        intact_end_splice=intact_end_splice
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
                        intact_end_splice=intact_end_splice
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
                        intact_end_splice=intact_end_splice
                    )
                    s += temp
                    new_exons.append((e, exon))
        else:
            raise AttributeError('breakpoint orientation must be specified to pull exons')
        if transcript.strand == STRAND.NEG:
            # reverse complement the sequence and reverse the exons
            temp = new_exons
            new_exons = []

            for ex, old_exon in temp[::-1]:
                e = Exon(
                    len(s) - ex.end + 1,
                    len(s) - ex.start + 1,
                    intact_start_splice=ex.intact_end_splice,
                    intact_end_splice=ex.intact_start_splice
                )
                new_exons.append((e, old_exon))
            s = reverse_complement(s)
        elif transcript.strand != STRAND.POS:
            raise AttributeError('transcript strand must be specified to pull exons')
        return s, new_exons


class Annotation(BreakpointPair):
    """
    a fusion of two transcripts created by the associated breakpoint_pair
    will also hold the other annotations for overlapping and encompassed and nearest genes
    """
    def __init__(
            self, bpp,
            transcript1=None,
            transcript2=None,
            data={},
            event_type=None,
            proximity=None):
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
            data=bpp.data,
            untemplated_sequence=bpp.untemplated_sequence
        )

        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.data.update(data)

        self.encompassed_genes = set()
        self.genes_proximal_to_break1 = set()
        self.genes_proximal_to_break2 = set()
        self.genes_overlapping_break1 = set()
        self.genes_overlapping_break2 = set()

        self.event_type = event_type if event_type is None else SVTYPE.enforce(event_type)
        self.proximity = proximity

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
            :class:`dict` of :class:`str` and :class:`str`: dictionary of attribute names and values
        """
        row = BreakpointPair.flatten(self)
        row.update({
            COLUMNS.genes_proximal_to_break1: self.genes_proximal_to_break1,
            COLUMNS.genes_proximal_to_break2: self.genes_proximal_to_break2,
            COLUMNS.gene1_direction: 'None',
            COLUMNS.gene2_direction: 'None'
        })
        if hasattr(self.transcript1, 'gene'):
            row[COLUMNS.gene1] = self.transcript1.gene.name
            row[COLUMNS.transcript1] = self.transcript1.name
            try:
                row[COLUMNS.gene1_direction] = str(determine_prime(self.transcript1, self.break1))
            except AttributeError:
                pass
        else:
            row[COLUMNS.gene1] = 'None'
            row[COLUMNS.transcript1] = '{}:{}_{}{}'.format(
                self.transcript1.reference_object,
                self.transcript1.start,
                self.transcript1.end,
                self.transcript1.strand)
        if hasattr(self.transcript2, 'gene'):
            row[COLUMNS.gene2] = self.transcript2.gene.name
            row[COLUMNS.transcript2] = self.transcript2.name
            try:
                row[COLUMNS.gene2_direction] = str(determine_prime(self.transcript2, self.break2))
            except AttributeError:
                pass
        else:
            row[COLUMNS.gene2] = 'None'
            row[COLUMNS.transcript2] = '{}:{}_{}{}'.format(
                self.transcript2.reference_object,
                self.transcript2.start,
                self.transcript2.end,
                self.transcript2.strand)
        row[COLUMNS.genes_encompassed] = ';'.join(sorted([x.name for x in self.encompassed_genes]))
        row[COLUMNS.genes_overlapping_break1] = ';'.join(sorted([x.name for x in self.genes_overlapping_break1]))
        row[COLUMNS.genes_overlapping_break2] = ';'.join(sorted([x.name for x in self.genes_overlapping_break2]))
        row[COLUMNS.genes_proximal_to_break1] = ';'.join(
            sorted(['{}({})'.format(x[0].name, x[1]) for x in self.genes_proximal_to_break1]))
        row[COLUMNS.genes_proximal_to_break2] = ';'.join(
            sorted(['{}({})'.format(x[0].name, x[1]) for x in self.genes_proximal_to_break2]))
        return row


def overlapping_transcripts(ref_ann, breakpoint):
    """
    Args:
        ref_ann (:class:`dict` of :class:`str` and :class:`list` of :any:`Gene`): the reference list of genes split
            by chromosome
        breakpoint (Breakpoint): the breakpoint in question
    Returns:
        :class:`list` of :any:`Transcript`: a list of possible transcripts
    """
    putative_annotations = []
    for gene in ref_ann[breakpoint.chr]:
        for transcript in gene.transcripts:
            if breakpoint.strand != STRAND.NS and transcript.strand != STRAND.NS \
                    and transcript.strand != breakpoint.strand:
                continue
            if Interval.overlaps(breakpoint, transcript):
                putative_annotations.append(transcript)
    return putative_annotations


def gather_breakpoint_annotations(ref_ann, breakpoint):
    """
    Args:
        ref_ann (:class:`dict` of :class:`str` and :class:`list` of :class:`Gene`): the reference annotations split
            into lists of genes by chromosome
        breakpoint (Breakpoint): the breakpoint annotations are to be gathered for
    """

    pos_overlapping_transcripts = []
    neg_overlapping_transcripts = []
    for gene in ref_ann[breakpoint.chr]:
        for t in gene.transcripts:
            if Interval.overlaps(t, breakpoint):
                if STRAND.compare(t.strand, STRAND.POS):
                    pos_overlapping_transcripts.append(t)
                if STRAND.compare(t.strand, STRAND.NEG):
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
        if breakpoint < first:
            temp.append(IntergenicRegion(breakpoint.chr, breakpoint[0], first[0] - 1, STRAND.NEG))
        if breakpoint[1] > last[1]:
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

    return (
        sorted(pos_overlapping_transcripts, key=lambda x: x.position),
        sorted(neg_overlapping_transcripts, key=lambda x: x.position))


def gather_annotations(ref, bp, event_type=None, proximity=None):  # TODO
    """
    each annotation is defined by the annotations selected at the breakpoints
    the other annotations are given relative to this
    the annotation at the breakpoint can be a transcript or an intergenic region

    Args:
        ref (:class:`dict` of :class:`str` and :class:`list` of :any:`Gene`): the list of reference genes hashed
            by chromosomes
        breakpoint_pairs (:class:`list` of :any:`BreakpointPair`): breakpoint pairs we wish to annotate as events

    Returns:
        :class:`list` of :class:`Annotation`: The annotations
    """
    annotations = dict()

    break1_pos, break1_neg = gather_breakpoint_annotations(ref, bp.break1)
    break2_pos, break2_neg = gather_breakpoint_annotations(ref, bp.break2)

    combinations = []

    if bp.stranded:
        if bp.break1.strand == STRAND.POS:
            if bp.break1.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_pos, break2_pos))
            else:
                combinations.extend(itertools.product(break1_pos, break2_neg))
        else:
            if bp.break1.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_neg, break2_pos))
            else:
                combinations.extend(itertools.product(break1_neg, break2_neg))
    else:
        if bp.opposing_strands:
            combinations.extend(itertools.product(break1_pos, break2_neg))
            combinations.extend(itertools.product(break1_neg, break2_pos))
        else:
            combinations.extend(itertools.product(break1_pos, break2_pos))
            combinations.extend(itertools.product(break1_neg, break2_neg))

    for a1, a2 in combinations:
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
        if (a1, a2) in annotations:  # ignore duplicates
            continue
        b1_itvl = bp.break1 & a1
        b2_itvl = bp.break2 & a2
        bpp = BreakpointPair.copy(bp)
        bpp.break1.start = b1_itvl[0]
        bpp.break1.end = b1_itvl[1]
        bpp.break2.start = b2_itvl[0]
        bpp.break2.end = b2_itvl[1]

        a = Annotation(bpp, a1, a2, event_type=event_type, proximity=proximity)

        for gene in ref[bp.break1.chr]:
            a.add_gene(gene)
        if bp.interchromosomal:
            for gene in ref[bp.break2.chr]:
                a.add_gene(gene)
        annotations[(a1, a2)] = a
    return list(annotations.values())
