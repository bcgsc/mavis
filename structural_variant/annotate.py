import TSV
from structural_variant.interval import Interval
from structural_variant.constants import *
import itertools
from structural_variant.error import *
from structural_variant.breakpoint import BreakpointPair, Breakpoint
from Bio import SeqIO
import re


class Annotation(BreakpointPair):
    """
    a fusion of two transcripts created by the associated breakpoint_pair
    will also hold the other annotations for overlapping and encompassed and nearest genes
    """
    def __init__(self, bpp, transcript1=None, transcript2=None, data={}, event_type=None):
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
        self.data = {}
        self.data.update(data)

        self.encompassed_genes = set()
        self.nearest_gene_break1 = set()
        self.nearest_gene_break2 = set()
        self.genes_at_break1 = set()
        self.genes_at_break2 = set()

        self.event_type = event_type if event_type is None else SVTYPE.enforce(event_type)

    def add_gene(self, gene):
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
            self.genes_at_break1.add(gene)
        if Interval.overlaps(gene, self.break2) and gene.chr == self.break2.chr \
                and gene != self.transcript2.reference_object:
            self.genes_at_break2.add(gene)

        if gene in self.genes_at_break1 or gene in self.genes_at_break2 or gene in self.encompassed_genes \
                or gene == self.transcript1.reference_object or gene == self.transcript2.reference_object:
            return

        d1 = Interval.dist(gene, self.break1)
        d2 = Interval.dist(gene, self.break2)

        if self.interchromosomal:
            if gene.chr == self.break1.chr:
                self.nearest_gene_break1.add((gene, d1))
            elif gene.chr == self.break2.chr:
                self.nearest_gene_break2.add((gene, d2))
        else:
            if d1 < 0:
                self.nearest_gene_break1.add((gene, d1))
            if d2 > 0:
                self.nearest_gene_break2.add((gene, d2))

        temp = set()

        tmin = [d for g, d in self.nearest_gene_break1 if d < 0]
        tmax = [d for g, d in self.nearest_gene_break1 if d > 0]
        tmin = 0 if len(tmin) == 0 else max(tmin)
        tmax = 0 if len(tmax) == 0 else min(tmax)

        for gene, dist in self.nearest_gene_break1:
            if tmin != 0 and dist == tmin:
                temp.add((gene, dist))
            elif tmax != 0 and dist == tmax:
                temp.add((gene, dist))

        self.nearest_gene_break1 = temp

        temp = set()

        tmin = [d for g, d in self.nearest_gene_break2 if d < 0]
        tmax = [d for g, d in self.nearest_gene_break2 if d > 0]
        tmin = 0 if len(tmin) == 0 else max(tmin)
        tmax = 0 if len(tmax) == 0 else min(tmax)

        for gene, dist in self.nearest_gene_break2:
            if tmin != 0 and dist == tmin:
                temp.add((gene, dist))
            elif tmax != 0 and dist == tmax:
                temp.add((gene, dist))

        self.nearest_gene_break2 = temp


class BioInterval:
    def __init__(self, reference_object, start, end, name=None):
        self.reference_object = reference_object
        self.name = name
        self.position = Interval(start, end)

    @property
    def start(self):
        return self.position.start

    @property
    def end(self):
        return self.position.end

    def __getitem__(self, index):
        return Interval.__getitem__(self, index)

    def __len__(self):
        return len(self.position)

    @property
    def key(self):
        return self.reference_object, self.position, self.name

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        else:
            return self.key == other.key

    def __hash__(self):
        return hash(self.key)


class IntergenicRegion(BioInterval):
    def __init__(self, chr, start, end, strand):
        BioInterval.__init__(self, chr, start, end)
        self.strand = STRAND.enforce(strand)

    @property
    def key(self):
        return self.reference_object, self.start, self.end, self.strand


class Gene(BioInterval):
    """
    """
    def __init__(self, chr, start, end, name=None, strand=STRAND.NS, aliases=[]):
        """
        Args:
            chr (str): the chromosome
            name (str): the gene name/id i.e. ENSG0001
            strand (STRAND): the genomic strand '+' or '-'
            aliases (List[str]): a list of aliases. For example the hugo name could go here
        Example:
            >>> Gene('X', 1, 1000, 'ENG0001', '+', ['KRAS'])
        """
        BioInterval.__init__(self, name=name, reference_object=chr, start=start, end=end)
        self.transcripts = set()
        self.strand = STRAND.enforce(strand)
        self.aliases = aliases

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    @property
    def key(self):
        return (self.name, self.strand, self.chr, self.start, self.end)


class FusionTranscript:
    def __init__(self):
        self.exon_mapping = {}
        self.exons = []
        self.sequence = ''

    @classmethod
    def determine_prime(cls, transcript, breakpoint):
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

    @classmethod
    def build(cls, ann, REFERENCE_GENOME):
        if not ann.transcript1 or not ann.transcript2:
            raise AttributeError('cannot produce fusion transcript for non-annotated fusions')
        elif not ann.event_type and ann.transcript1 == ann.transcript2:
            raise AttributeError('event_type must be specified to produce a fusion transcript')
        elif ann.untemplated_sequence is None:
            raise AttributeError(
                'cannot build a fusion transcript where the untemplated sequence has not been specified')

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
                window_seq = REFERENCE_GENOME[b1.chr].seq[window.start - 1:window.end]
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
            t1 = cls.determine_prime(ann.transcript1, ann.break1)
            t2 = cls.determine_prime(ann.transcript2, ann.break2)

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
        return ft

    def _splice_patterns(self):
        """
        returns a list of splice sites to be connected as a splicing pattern

        1 or two splice sites are abrogated
        """
        splice_site_sets = [[]]

        i = 1
        while i < len(self.exons):
            exon = self.exons[i]
            prev = self.exons[i - 1]
            # if the last nearest splice site is intact then add this splice site and it
            if not exon.intact_start_splice and not exon.intact_end_splice:
                if i < len(self.exons) - 1:
                    nexxt = self.exons[i + 1]
                    for s in splice_site_sets:
                        s.extend([prev.end, nexxt.start])
                    i += 1
            elif prev.intact_end_splice:
                if exon.intact_start_splice:  # both intact
                    for s in splice_site_sets:
                        s.extend([prev.end, exon.start])
                else:  # previous is intact but the current is not
                    # two options, can skip this exon or retain the previous intron
                    # duplicate the current sets
                    if i < len(self.exons) - 1:
                        temp = []
                        for s in splice_site_sets:
                            # skip this exon
                            nexxt = self.exons[i + 1]
                            skip = s + [prev.end, nexxt.start]
                            temp.append(skip)
                            retain = s + [exon.end, nexxt.start]
                            # retain the previous intron
                            temp.append(retain)
                        i += 1
                        splice_site_sets = temp
            else:
                if exon.intact_start_splice:
                    # two options: retain previous intron, or skip previous exon
                    temp = []
                    for s in splice_site_sets:
                        skip = s[:]
                        skip[-1] = exon.start
                        temp.append(skip)
                        temp.append(s)
                    splice_site_sets = temp
                else:  # both abrogated retain intron, essentially do nothing
                    pass
            i += 1
        return splice_site_sets

    def decompose(self):
        """
        return a list of putative transcripts from the original transcript
        """
        pass

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
                    print('exon', exon.start, exon.end, exon)
                    print('breakpoint', breakpoint)
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
                    print('new exon', e.start, e.end, intact_start_splice, intact_end_splice)
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

    def translate(self, splicing_pattern=[]):
        if len(splicing_pattern) % 2 != 0:
            raise AttributeError('splice sites must be an even number')
        for s in splicing_pattern:
            if s < 1 or s > len(self.sequence):
                raise AttributeError('splice points outside the transcript')
        ranges = [1] + sorted(splicing_pattern) + [len(self.sequence)]
        cdna = []
        exons = []

        for i in range(0, len(ranges), 2):
            s = ranges[i]
            t = ranges[i + 1]
            exons.append((s, t))
            cdna.append(self.sequence[s - 1:t])
        cdna = ''.join(cdna)
        print(exons)
        translations = []  # (cds_start, cds_end)
        for i in range(0, CODON_SIZE):
            aa_sequence = translate(cdna, i)
            # now calc the open reading frames
            temp = []
            last_met_pos = -1
            for aa_num, char in enumerate(aa_sequence):
                if char == START_AA:
                    if last_met_pos < 0:
                        last_met_pos = aa_num
                elif char == STOP_AA:
                    if last_met_pos >= 0:
                        temp.append((last_met_pos, aa_num))
                        last_met_pos = -1

            for start, end in temp:
                translations.append((start * CODON_SIZE + i + 1, (end + 1) * CODON_SIZE + i))

        transcripts = []
        for cds_start, cds_end in sorted(translations):
            print(cds_start, cds_end)
            t = Transcript(cds_start, cds_end, gene=self, exons=exons, strand=STRAND.POS)
            transcripts.append(t)
        return transcripts


class Transcript(BioInterval):
    """
    """
    def __init__(
        self,
        cds_start=None,
        cds_end=None,
        genomic_start=None,
        genomic_end=None,
        gene=None,
        name=None,
        strand=None,
        exons=[],
        domains=[]
    ):
        """ creates a new transcript object

        Args:
            gene (Gene, optional): the gene this transcript belongs to
            name (str, optional): the name of the transcript or external db id. For example ENTS0001
            cds_start (int): the position (wrt the first exon) where translation would begin
            cds_cdna = self.transcript.convert_genomic_to_cdna(self.start)
        cdna = abs(self.transcript.cds_start - 1 - cdna)
        return cdna % 3end (int): the position (wrt the first exon) where translation would terminate
            strand (STRAND, optional): the strand the transcript occurs on

        """
        if genomic_start is None and len(exons) > 0:
            genomic_start = min([e[0] for e in exons])
        if genomic_end is None and len(exons) > 0:
            genomic_end = max([e[1] for e in exons])

        BioInterval.__init__(self, reference_object=gene, name=name, start=genomic_start, end=genomic_end)

        self.exons = []
        self.domains = set()
        self._strand = strand

        if self._strand and self.gene and hasattr(self.gene, 'strand') and self.gene.strand != self._strand:
            raise AttributeError('strand does not match reference object')

        try:
            self.cds_start = int(cds_start) if cds_start is not None else None
            self.cds_end = int(cds_end) if cds_end is not None else None
            if cds_start is not None and cds_end is not None:
                if self.cds_start > self.cds_end:
                    raise AttributeError('cds_end must be >= cds_start')
        except ValueError:
            raise AttributeError('cds_end and/or cds_start must be integers')

        for e in exons:
            if isinstance(e, Exon):
                if e not in self.exons:
                    self.exons.append(e)
                    e.reference_object = self
            else:
                Exon(e[0], e[1], transcript=self)

        if self.gene is not None and hasattr(self.gene, 'transcripts'):
            self.gene.transcripts.add(self)

        exons = sorted(self.exons, key=lambda x: x.start)
        for i, exon in enumerate(exons):
            if i == 0:
                continue
            previous = exons[i - 1]
            if Interval.overlaps((previous.start, previous.end), (exon.start, exon.end)):
                raise AttributeError('exons cannot overlap within a transcript')

        for d in domains:
            d.reference_object = self
            self.domains.add(d)

    def reading_frame(self):
        """
        returns 0 if the reading frame begins on the first base of the sequence
        1 if the second, and so on up to 2
        """
        if self.cds_start is None:
            raise AttributeError('cannot calculate reading frame if the translation start has not been given')
        cdna = 1
        cdna = abs(self.cds_start - 1 - cdna)
        p = (cdna % CODON_SIZE) + CODON_SIZE - 1
        return p % CODON_SIZE

    @property
    def gene(self):
        """(:class:`~structural_variant.annotate.Gene`): the gene this transcript belongs to"""
        return self.reference_object

    @property
    def strand(self):
        if self._strand is not None:
            return self._strand
        return self.gene.strand

    def genomic_length(self):
        return len(self.position)

    @property
    def genomic_start(self):
        return self.start

    def genomic_utr_regions(self):
        utr = []
        if self.strand not in [STRAND.POS, STRAND.NEG]:
            raise AttributeError('strand must be positive or negative to calculate regions')

        exons = sorted(self.exons, key=lambda x: x.start)

        if self.cds_start is not None:
            if self.strand == STRAND.POS:
                utr.append(Interval(exons[0].start, self.convert_cdna_to_genomic(self.cds_start)))
            else:
                utr.append(Interval(self.convert_cdna_to_genomic(self.cds_start), exons[-1].end))
        if self.cds_end is not None:
            if self.strand == STRAND.POS:
                utr.append(Interval(self.convert_cdna_to_genomic(self.cds_end), exons[-1].end))
            else:
                utr.append(Interval(exons[0].start, self.convert_cdna_to_genomic(self.cds_end)))
        return utr

    @property
    def genomic_end(self):
        return self.end

    def _genomic_to_cdna_mapping(self):
        mapping = {}

        exons = sorted(self.exons, key=lambda x: x.start)
        if self.strand == STRAND.POS:
            pass
        elif self.strand == STRAND.NEG:
            exons.reverse()
        else:
            raise StrandSpecificityError('cannot convert without strand information')

        l = 1
        for e in exons:
            mapping[Interval(e.start, e.end)] = Interval(l, l + len(e) - 1)
            l += len(e)
        return mapping

    def _cdna_to_genomic_mapping(self):
        mapping = {}
        for k, v in self._genomic_to_cdna_mapping().items():
            mapping[v] = k
        return mapping

    def convert_genomic_to_cdna(self, pos):
        mapping = self._genomic_to_cdna_mapping()
        return Interval.convert_pos(mapping, pos)

    def convert_cdna_to_genomic(self, pos):
        mapping = self._cdna_to_genomic_mapping()
        return Interval.convert_pos(mapping, pos)

    def convert_aa_to_cdna(self, pos):
        return Interval(self.cds_start - 1 + (pos - 1) * 3 + 1, self.cds_start - 1 + pos * 3)

    @property
    def key(self):
        return (self.gene, self.name, self.start, self.end, self.cds_start, self.cds_end)

    def exon_number(self, exon):
        for i, e in enumerate(sorted(exons)):
            if e is exon:
                return i
        raise AttributeError('can only calculate phase on associated exons')


class Exon(BioInterval):
    """
    """
    def __init__(self, start, end, transcript=None, name=None, intact_start_splice=True, intact_end_splice=True):
        """
        Args:
            start (int): the genomic start position
            end (int): the genomic end position
            name (str, optional): the name of the exon
            transcript (Transcript, optional): the 'parent' transcript this exon belongs to
        Raises:
            AttributeError: if the exon start > the exon end
        Example:
            >>> Exon(15, 78)
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end)
        if self.transcript is not None:
            if self not in self.transcript.exons:
                self.transcript.exons.append(self)
        self.intact_start_splice = intact_start_splice
        self.intact_end_splice = intact_end_splice
        if end - start + 1 < SPLICE_SITE_RADIUS * sum([intact_start_splice, intact_end_splice]):
            raise AttributeError('exons must be greater than double the length of a splice site')

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.Transcript`): the transcript this exon belongs to"""
        return self.reference_object

    def __hash__(self):
        return hash((self.transcript, self.start, self.end, self.name))

    @property
    def start_phase(self):
        if self.transcript is None:
            raise AttributeError('cannot calculate exon phase unless there is an associated transcript')
        cdna = self.transcript.convert_genomic_to_cdna(self.start)
        cdna = abs(self.transcript.cds_start - 1 - cdna)
        return cdna % CODON_SIZE

    @property
    def end_phase(self):
        if self.transcript is None:
            raise AttributeError('cannot calculate exon phase unless there is an associated transcript')
        cdna = self.transcript.convert_genomic_to_cdna(self.end)
        cdna = abs(self.transcript.cds_start - 1 - cdna)
        return cdna % CODON_SIZE

    @property
    def start_splice_site(self):
        return Interval(self.start - SPLICE_SITE_RADIUS, self.start + SPLICE_SITE_RADIUS - 1)

    @property
    def end_splice_site(self):
        return Interval(self.end - SPLICE_SITE_RADIUS + 1, self.end + SPLICE_SITE_RADIUS)

class Domain:
    """
    """
    def __init__(self, name, regions, transcript=None):
        """
        Args:
            name (str): the name of the domain i.e. PF00876
            regions (List[Interval]): the amino acid ranges that are part of the domain
            transcript (Transcript, optional): the 'parent' transcript this domain belongs to
        Raises:
            AttributeError: if the end of any region is less than the start
        Example:
            >>> Domain('DNA binding domain', [(1, 4), (10, 24)], transcript)
        """
        self.reference_object = transcript
        self.name = name
        self.regions = sorted(list(set(regions)))  # remove duplicates

        for i, region in enumerate(self.regions):
            if region[0] > region[1]:
                raise AttributeError('domain region start must be <= end')
        self.regions = Interval.min_nonoverlapping(*self.regions)
        if self.reference_object is not None:
            self.reference_object.domains.add(self)

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.Transcript`): the transcript this domain belongs to"""
        return self.reference_object

    @property
    def key(self):
        return tuple([self.name, self.transcript] + self.regions)


def load_masking_regions(filepath):
    """
    reads a file of regions. The expect input format for the file is tab-delimited and

    +---------------+---------------+-----------------------+
    | column name   | value type    | description           |
    +===============+===============+=======================+
    | chr           | string        | the chromosome name   |
    +---------------+---------------+-----------------------+
    | start         | int           | the start position    |
    +---------------+---------------+-----------------------+
    | end           | int           | the end position      |
    +---------------+---------------+-----------------------+
    | name          | string        | label for the region  |
    +---------------+---------------+-----------------------+

    Args:
        filepath (string): path to the input tab-delimited file
    Returns:
        Dict[str,List[BioInterval]]:
            a dictionary keyed by chromosome name with values of lists of regions on the chromosome

    """
    header, rows = TSV.read_file(
        filepath,
        require=['chr', 'start', 'end', 'name'],
        cast={'start': int, 'end': int, 'chr': lambda x: re.sub('^chr', '', x)}
    )
    regions = {}
    for row in rows:
        r = BioInterval(reference_object=row['chr'], start=row['start'], end=row['end'], name=row['name'])
        regions.setdefault(r.reference_object, []).append(r)
    return regions


def load_reference_genes(filepath):
    """
    given a file in the std input format (see below) reads and return a list of genes (and sub-objects)

    +-----------------------+---------------------------+-----------------------------------------------------------+
    | column name           | example                   | description                                               |
    +=======================+===========================+===========================================================+
    | ensembl_transcript_id | ENST000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | ensembl_gene_id       | ENSG000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | strand                | -1                        | positive or negative 1                                    |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_start     | 44                        | where translation begins relative to the start of the cdna|
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_end       | 150                       | where translation terminates                              |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | genomic_exon_ranges   | 100-201;334-412;779-830   | semi-colon demitited exon start/ends                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | AA_domain_ranges      | DBD:220-251,260-271       | semi-colon delimited list of domains                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | hugo_names            | KRAS                      | hugo gene name                                            |
    +-----------------------+---------------------------+-----------------------------------------------------------+

    Args:
        filepath (string): path to the input tab-delimited file

    Returns:
        Dict[str,List[Gene]]: a dictionary keyed by chromosome name with values of list of genes on the chromosome
    """
    def parse_exon_list(row):
        if not row:
            return []
        exons = []
        for temp in row.split(';'):
            try:
                s, t = temp.split('-')
                exons.append(Exon(int(s), int(t)))
            except:
                print('exon error:', temp, row)
        return exons

    def parse_domain_list(row):
        if not row:
            return []
        domains = []
        for d in row.split(';'):
            try:
                name, temp = d.rsplit(':')
                temp = temp.split(',')
                temp = [x.split('-') for x in temp]
                temp = [(int(x), int(y)) for x, y in temp]
                d = Domain(name, temp)
            except:
                print('error in domain:', d, row)
        return domains

    def nullable_int(row):
        try:
            row = int(row)
        except ValueError:
            row = TSV.null(row)
        return row

    def parse_strand(row):
        if row in ['-1', '-']:
            return STRAND.NEG
        elif row in ['1', '+', '+1']:
            return STRAND.POS
        raise ValueError('cast to strand failed')

    header, rows = TSV.read_file(
        filepath,
        require=[
            'ensembl_gene_id',
            'chr',
            'ensembl_transcript_id'
        ],
        add={
            'cdna_coding_start': 'null',
            'cdna_coding_end': 'null',
            'AA_domain_ranges': '',
            'genomic_exon_ranges': '',
            'hugo_names': '',
            'transcript_genomic_start': 'null',
            'transcript_genomic_end': 'null'
        },
        cast={
            'genomic_exon_ranges': parse_exon_list,
            'AA_domain_ranges': parse_domain_list,
            'cdna_coding_end': nullable_int,
            'cdna_coding_start': nullable_int,
            'transcript_genomic_end': nullable_int,
            'transcript_genomic_start': nullable_int,
            'strand': parse_strand,
            'gene_start': int,
            'gene_end': int
        }
    )
    genes = {}
    for row in rows:
        g = Gene(
            row['chr'],
            row['gene_start'],
            row['gene_end'],
            name=row['ensembl_gene_id'],
            strand=row['strand'],
            aliases=row['hugo_names'].split(';')
        )
        if g.name in genes:
            g = genes[g.name]
        else:
            genes[g.name] = g

        t = Transcript(
            name=row['ensembl_transcript_id'],
            gene=g,
            genomic_start=row['transcript_genomic_start'],
            genomic_end=row['transcript_genomic_end'],
            exons=row['genomic_exon_ranges'],
            domains=row['AA_domain_ranges'],
            cds_start=row['cdna_coding_start'],
            cds_end=row['cdna_coding_end']
        )

    ref = {}
    for g in genes.values():
        if g.chr not in ref:
            ref[g.chr] = []
        ref[g.chr].append(g)
    return ref


def overlapping_transcripts(ref_ann, breakpoint):
    """
    Args:
        ref_ann (Dict[str,List[Gene]]): the reference list of genes split by chromosome
        breakpoint (Breakpoint): the breakpoint in question
    Returns:
        List[Transcript]: a list of possible transcripts
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

    return sorted(pos_overlapping_transcripts, key=lambda x: x.position), sorted(neg_overlapping_transcripts, key=lambda x: x.position)


def gather_annotations(ref, bp):  # TODO
    """
    Args:
        ref (Dict[str,List[Gene]]): the list of reference genes hashed by chromosomes
        breakpoint_pairs (List[BreakpointPair]): breakpoint pairs we wish to annotate as events

    each annotation is defined by the annotations selected at the breakpoints
    the other annotations are given relative to this
    the annotation at the breakpoint can be a transcript or an intergenic region

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

        a = Annotation(bpp, a1, a2)

        for gene in ref[bp.break1.chr]:
            a.add_gene(gene)
        if bp.interchromosomal:
            for gene in ref[bp.break2.chr]:
                a.add_gene(gene)
        annotations[(a1, a2)] = a
    return list(annotations.values())


def load_reference_genome(filename):
    HUMAN_REFERENCE_GENOME = None
    with open(filename, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    return HUMAN_REFERENCE_GENOME
