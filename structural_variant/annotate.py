import TSV
from .interval import Interval
from .constants import *
import itertools
from .error import *
from .breakpoint import BreakpointPair, Breakpoint
from Bio import SeqIO
import re
import warnings


def calculate_ORF(spliced_cdna_sequence, min_orf_size=None):
    """
    calculate all possible open reading frames given a spliced cdna sequence (no introns)

    Args:
        spliced_cdna_sequence (str): the sequence

    Returns:
        :any:`list` of :any:`Interval`: list of open reading frame positions on the input sequence
    """
    # do not revcomp
    assert(START_AA != STOP_AA)
    cds_orfs = []  # (cds_start, cds_end)
    for offset in range(0, CODON_SIZE):
        aa_sequence = translate(spliced_cdna_sequence, offset)
        # now calc the open reading frames
        starts = []
        stops = []
        for i, aa in enumerate(aa_sequence):
            if aa == START_AA:
                starts.append(i * CODON_SIZE + 1)
            elif aa == STOP_AA:
                stops.append((i + 1) * CODON_SIZE)

        orfs = []

        for s in starts:
            for t in sorted(stops):
                if t > s:
                    i = Interval(s, t)
                    if not min_orf_size or len(i) >= min_orf_size:
                        orfs.append(Interval(s, t))
                    break

        temp = {}
        for orf in orfs:
            if orf.end not in temp:
                temp[orf.end] = orf
            elif len(orf) > len(temp[orf.end]):
                temp[orf.end] = orf

        cds_orfs.extend(temp.values())
    return cds_orfs


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


class BioInterval:
    def __init__(self, reference_object, start, end, name=None):
        """
        Args:
            reference_object: the object this interval is on
            start (int) start of the interval (inclusive)
            end (int): end of the interval (inclusive)
            name: optional

        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> b[0]
            12572784
            >>> b[1]
            12578898
        """
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
        """
        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> len(b)
            6115
        """
        return len(self.position)

    @property
    def key(self):
        return self.reference_object, self.position, self.name

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        else:
            return self.key == other.key

    def __lt__(self, other):
        if other.reference_object != self.reference_object:
            if self.reference_object < other.reference_object:
                return True
            else:
                return False
        elif self.position < other.position:
            return True
        else:
            return False

    def __hash__(self):
        return hash(self.key)

    def get_sequence(self, REFERENCE_GENOME=None):
        raise NotImplementedError('abstract method must be overidden')


class IntergenicRegion(BioInterval):
    def __init__(self, chr, start, end, strand):
        """
        Args:
            chr (str): the reference object/chromosome for this region
            start (int): the start of the IntergenicRegion
            end (int): the end of the IntergenicRegion
            strand (STRAND): the strand the region is defined on

        Example:
            >>> IntergenicRegion('1', 1, 100, '+')
        """
        BioInterval.__init__(self, chr, start, end)
        self.strand = STRAND.enforce(strand)

    @property
    def key(self):
        return self.reference_object, self.start, self.end, self.strand


class Gene(BioInterval):
    """
    """
    def __init__(self, chr, start, end, name=None, strand=STRAND.NS, aliases=None, sequence=None):
        """
        Args:
            chr (str): the chromosome
            name (str): the gene name/id i.e. ENSG0001
            strand (STRAND): the genomic strand '+' or '-'
            aliases (:class:`list` of :class:`str`): a list of aliases. For example the hugo name could go here
            sequence (str): genomic sequence of the gene
        Example:
            >>> Gene('X', 1, 1000, 'ENG0001', '+', ['KRAS'])
        """
        aliases = [] if aliases is None else aliases
        BioInterval.__init__(self, name=name, reference_object=chr, start=start, end=end)
        self.transcripts = set()
        self.strand = STRAND.enforce(strand)
        self.aliases = aliases
        self.sequence = sequence

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    @property
    def key(self):
        return (self.name, self.strand, self.chr, self.start, self.end)

    def get_sequence(self, REFERENCE_GENOME):
        """
        gene sequence is always given wrt to the positive forward strand regardless of gene strand

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            str: the sequence of the gene
        """
        if self.sequence:
            return self.sequence
        elif REFERENCE_GENOME is None:
            raise AttributeError('reference genome is required to retrieve the gene sequence')
        else:
            return str(REFERENCE_GENOME[self.chr].seq[self.start - 1:self.end]).upper()


class Translation(BioInterval):
    def __init__(self, start, end, transcript, splicing_pattern, domains=None, sequence=None):
        """
        describes the splicing pattern and cds start and end with reference to a particular transcript

        Args:
            start (int): start of the coding sequence (cds) relative to the start of the first exon in the transcript
            end (int): end of the coding sequence (cds) relative to the start of the first exon in the transcript
            transcript (Transcript): the transcript this is a Translation of
            splicing_pattern (:class:`list` of :any:`int`): a list of splicing positions to be used
            domains (:class:`list` of :any:`Domain`): a list of the domains on this translation
            sequence (str): the cds sequence
        """
        domains = [] if domains is None else domains
        BioInterval.__init__(self, reference_object=transcript, name=None, start=start, end=end)
        self.splicing_pattern = sorted(splicing_pattern)
        self.domains = []
        self.sequence = sequence
        for d in domains:
            self.add_domain(d)
        if self.transcript:
            self.transcript.add_translation(self)

    @property
    def transcript(self):
        return self.reference_object

    def convert_aa_to_cdna(self, pos):
        """
        Args:
            pos (int): the amino acid position

        Returns:
            int: the cdna equivalent
        """
        return Interval(self.start - 1 + (pos - 1) * 3 + 1, self.start - 1 + pos * 3)

    def add_domain(self, domain):
        """
        Args:
            domain (Domain): the domain to be added
        """
        domain.reference_object = self
        if domain not in self.domains:
            self.domains.append(domain)

    def genomic_utr_regions(self):
        utr = []
        if self.transcript.strand not in [STRAND.POS, STRAND.NEG]:
            raise AttributeError('strand must be positive or negative to calculate regions')

        exons = sorted(self.transcript.exons, key=lambda x: x.start)

        if self.transcript.strand == STRAND.POS:
            utr.append(Interval(exons[0].start, self.transcript.convert_cdna_to_genomic(self.start)))
        else:
            utr.append(Interval(self.transcript.convert_cdna_to_genomic(self.start), exons[-1].end))
        if self.transcript.strand == STRAND.POS:
            utr.append(Interval(self.transcript.convert_cdna_to_genomic(self.end), exons[-1].end))
        else:
            utr.append(Interval(exons[0].start, self.transcript.convert_cdna_to_genomic(self.end)))
        return utr

    def get_sequence(self, REFERENCE_GENOME=None):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            str: the cds sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        if self.sequence:
            return self.sequence
        elif self.transcript and self.transcript.strand:
            seq = self.transcript.get_spliced_cdna_sequence(self.splicing_pattern, REFERENCE_GENOME)
            return seq[self.start - 1:self.end]
        raise AttributeError('insufficient sequence information')

    def get_cds_sequence(self, REFERENCE_GENOME=None):
        """
        wrapper for the sequence method

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name
        """
        return self.get_sequence(REFERENCE_GENOME)

    def get_AA_sequence(self, REFERENCE_GENOME=None):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            str: the amino acid sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        cds = self.get_cds_sequence(REFERENCE_GENOME)
        return translate(cds)

    @property
    def key(self):
        return (self.reference_object, self.start, self.end, self.splicing_pattern)


class Transcript(BioInterval):
    """
    """
    def __init__(
        self,
        cds_start=None,
        cds_end=None,
        exons=None,
        genomic_start=None,
        genomic_end=None,
        gene=None,
        name=None,
        strand=None,
        domains=None,
        translations=None,
        sequence=None
    ):
        """ creates a new transcript object

        Args:
            cds_start (int): the start of the coding sequence relative to the start of the first exon
            cds_end (int): the end of the coding sequence relative to the start of the first exon
            exons (:class:`list` of :any:`Exon`): list of Exon that make up the transcript
            genomic_start (int): genomic start position of the transcript
            genomic_end (int): genomic end position of the transcript
            gene (Gene): the gene this transcript belongs to
            name (str): name of the transcript
            strand (STRAND): strand the transcript is on, defaults to the strand of the Gene if not specified
            domains (:class:`list` of :any:`Domain`): list of domains to add to translations
            translations (:class:`list` of :any:`Translation`): Translation associated with this transcript
            sequence (str): unspliced cDNA sequence
        """
        # cannot use mutable default args in the function decl
        exons = [] if exons is None else exons
        domains = [] if domains is None else domains
        translations = [] if translations is None else translations

        if genomic_start is None and len(exons) > 0:
            genomic_start = min([e[0] for e in exons])
        if genomic_end is None and len(exons) > 0:
            genomic_end = max([e[1] for e in exons])

        BioInterval.__init__(self, reference_object=gene, name=name, start=genomic_start, end=genomic_end)

        self.exons = []
        self._strand = strand
        self.translations = translations
        self.sequence = sequence

        if self._strand and self.gene and hasattr(self.gene, 'strand') and self.gene.strand != self._strand:
            raise AttributeError('strand does not match reference object')

        for e in exons:
            self.add_exon(e)

        if self.gene is not None and hasattr(self.gene, 'transcripts'):
            self.gene.transcripts.add(self)

        if cds_start is not None or cds_end is not None:
            for sp in self.splicing_patterns():
                tl = Translation(cds_start, cds_end, self, splicing_pattern=sp, domains=domains)

    def splicing_patterns(self):
        """
        returns a list of splice sites to be connected as a splicing pattern

        Returns:
            :class:`list` of :any:`int`: List of positions to be spliced together
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
        """(:any:`Gene`): the gene this transcript belongs to"""
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
        """
        Args:
            pos (int): the genomic position to be converted

        Returns:
            int: the cdna equivalent

        Raises:
            :class:`~structural_variant.error.DiscontinuousMappingError`: when a genomic position not present in the cdna is attempted to be converted
        """
        mapping = self._genomic_to_cdna_mapping()
        return Interval.convert_pos(mapping, pos)

    def convert_cdna_to_genomic(self, pos):
        """
        Args:
            pos (int): cdna position

        Returns:
            int: the genomic equivalent
        """
        mapping = self._cdna_to_genomic_mapping()
        return Interval.convert_pos(mapping, pos)

    def add_exon(self, exon):
        """
        adds a given exon to the transcript

        Args:
            exon (Exon): the exon to be added

        Raises:
            AttributeError: if the input exon overlaps any of the existing exons
        """
        if not isinstance(exon, Exon):
            exon = Exon(exon[0], exon[1])

        for e in self.exons:
            if Interval.overlaps(e, exon):
                raise AttributeError('exons cannot overlap', (e[0], e[1]), (exon[0], exon[1]))

        self.exons.append(exon)
        exon.reference_object = self

        self.exons = sorted(self.exons, key=lambda x: x.start)

    def add_translation(self, translation):
        """
        adds a translation to the transcript

        Args:
            translation (Translation): the translation to be added
        """
        translation.reference_object = self
        if translation not in self.translations:
            self.translations.append(translation)

    @property
    def key(self):
        return (self.gene, self.name, self.start, self.end)

    def exon_number(self, exon):
        """
        exon numbering is based on the direction of translation

        Args:
            exon (Exon): the exon to be numbered

        Returns:
            int: the exon number (1 based)

        Raises:
            AttributeError: if the strand is not given or the exon does not belong to the transcript
        """
        for i, e in enumerate(self.exons):
            if exon != e:
                continue
            if self.strand == STRAND.POS:
                return i + 1
            elif self.strand == STRAND.NEG:
                return len(self.exons) - i
            else:
                raise AttributeError('strand must be pos or neg to calculate the exon number')
        raise AttributeError('can only calculate phase on associated exons')

    def get_sequence(self, REFERENCE_GENOME=None):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            str: the DNA sequence of the transcript including introns
        """
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
        elif REFERENCE_GENOME is None:
            raise AttributeError('reference genome is required to retrieve the gene sequence')
        else:
            if self.strand == STRAND.NEG:
                return reverse_complement(REFERENCE_GENOME[self.gene.chr].seq[self.start - 1:self.end])
            else:
                return str(REFERENCE_GENOME[self.gene.chr].seq[self.start - 1:self.end])

    def get_spliced_cdna_sequence(self, splicing_pattern, REFERENCE_GENOME=None):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): the list of splicing positions
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            str: the spliced cDNA sequence
        """
        temp = sorted([self.start] + splicing_pattern + [self.end])
        m = min(temp)
        conti = []
        for i in range(0, len(temp) - 1, 2):
            conti.append(Interval(temp[i] - m, temp[i + 1] - m))
        seq = self.get_sequence(REFERENCE_GENOME)
        if self.strand == STRAND.NEG:
            # adjust the continuous intervals for the min and flip if revcomp
            seq = reverse_complement(seq)
        spliced_seq = ''.join([str(seq[i.start:i.end + 1]) for i in conti])
        return spliced_seq if self.strand == STRAND.POS else reverse_complement(spliced_seq)


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
    def build(cls, ann, REFERENCE_GENOME, min_orf_size=None, max_orf_cap=None):
        """
        Args:
            ann (Annotation): the annotation object we want to build a FusionTranscript for
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

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
            spliced_cdna_seq = ft.get_spliced_cdna_sequence(spl_patt, REFERENCE_GENOME)
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
                tl = Translation(orf.start, orf.end, ft, spl_patt)

        # remap the domains from the original translations to the current translations
        for new_tl in ft.translations:
            aa_seq = new_tl.get_AA_sequence(REFERENCE_GENOME)
            for tl in ann.transcript1.translations + ann.transcript2.translations:
                for dom in tl.domains:
                    try:
                        regions = dom.align_seq(aa_seq, REFERENCE_GENOME)
                        Domain(dom.name, regions, new_tl)
                    except UserWarning:
                        pass
        return ft

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


class Exon(BioInterval):
    """
    """
    def __init__(self, start, end, transcript=None, name=None, intact_start_splice=True, intact_end_splice=True):
        """
        Args:
            start (int): the genomic start position
            end (int): the genomic end position
            name (str): the name of the exon
            transcript (Transcript): the 'parent' transcript this exon belongs to
        Raises:
            AttributeError: if the exon start > the exon end
        Example:
            >>> Exon(15, 78)
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end)
        if self.transcript is not None:
            self.transcript.add_exon(self)
        self.intact_start_splice = intact_start_splice
        self.intact_end_splice = intact_end_splice
        if end - start + 1 < SPLICE_SITE_RADIUS * sum([intact_start_splice, intact_end_splice]):
            warnings.warn('exons must be greater than double the length of a splice site')

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.Transcript`): the transcript this exon belongs to"""
        return self.reference_object

    def __hash__(self):
        return hash((self.transcript, self.start, self.end, self.name))

    @property
    def start_splice_site(self):
        """(:class:`~structural_variant.interval.Interval`): the genomic range describing the splice site"""
        return Interval(self.start - SPLICE_SITE_RADIUS, self.start + SPLICE_SITE_RADIUS - 1)

    @property
    def end_splice_site(self):
        """(:class:`~structural_variant.interval.Interval`): the genomic range describing the splice site"""
        return Interval(self.end - SPLICE_SITE_RADIUS + 1, self.end + SPLICE_SITE_RADIUS)


class DomainRegion(Interval):
    def __init__(self, start, end, seq=None):
        Interval.__init__(self, start, end)
        self.sequence = seq
        if seq and len(seq) != len(self):
            raise AttributeError('domain region sequence must be of equal length', self, seq)


class Domain:
    """
    """
    def __init__(self, name, regions, translation=None):
        """
        Args:
            name (str): the name of the domain i.e. PF00876
            regions (:any:`list` of :any:`DomainRegion`): the amino acid ranges that are part of the domain
            transcript (Transcript): the 'parent' transcript this domain belongs to
        Raises:
            AttributeError: if the end of any region is less than the start
        Example:
            >>> Domain('DNA binding domain', [(1, 4), (10, 24)], transcript)
        """
        self.reference_object = translation
        self.name = name
        self.regions = sorted(list(set(regions)))  # remove duplicates
        if len(regions) == 0:
            raise AttributeError('at least one region must be given')
        for r1, r2 in itertools.combinations(self.regions, 2):
            if Interval.overlaps(r1, r2):
                raise AttributeError('regions cannot overlap')

        for i in range(0, len(self.regions)):
            curr = self.regions[i]
            if not hasattr(curr, 'sequence'):
                self.regions[i] = DomainRegion(curr[0], curr[1])

        if self.translation is not None:
            self.translation.add_domain(self)

    @property
    def translation(self):
        """(:class:`~structural_variant.annotate.Translation`): the Translation this domain belongs to"""
        return self.reference_object

    @property
    def key(self):
        return tuple([self.name, self.translation])

    def score_region_mapping(self, REFERENCE_GENOME=None):
        """
        compares the sequence in each DomainRegion to the sequence collected for that domain region from the translation object

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            tuple of int and int: tuple contains

                - int: the number of matching amino acids
                - int: the total number of amino acids
        """
        if self.translation:
            aa = self.translation.get_AA_sequence(REFERENCE_GENOME)
            total = 0
            matches = 0
            for region in self.regions:
                if not region.sequence:
                    raise AttributeError('insufficient sequence information')
                ref = aa[region.start - 1:region.end]
                for c1, c2 in zip(ref, region.sequence):
                    if c1 == c2:
                        matches += 1
                    total += 1
            return matches, total
        else:
            raise AttributeError('insufficient sequence information')

    def get_sequences(self, REFERENCE_GENOME=None):
        """
        returns the amino acid sequences for each of the domain regions associated with
        this domain in the order of the regions (sorted by start)

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            :class:`list` of :class:`str`: list of amino acid sequences for each DomainRegion

        Raises:
            AttributeError: if there is not enough sequence information given to determine this
        """
        sequences = {}
        for region in self.regions:
            if region.sequence:
                sequences[region] = region.sequence
        if any([r not in sequences for r in self.regions]):
            if self.translation:
                aa = self.translation.get_AA_sequence(REFERENCE_GENOME)
                for region in self.regions:
                    s = aa[region.start - 1:region.end]
                    if region not in sequences:
                        sequences[region] = s
            else:
                raise AttributeError(
                    'Insufficient information to gather sequence data. reference_object must have sequence information or '
                    'REFERENCE_GENOME must be provided')
        return [sequences[r] for r in self.regions]

    def align_seq(self, input_sequence, REFERENCE_GENOME=None):
        """
        align each region to the input sequence starting with the last one.
        then take the subset of sequence that remains to align the second last and so on
        return a list of intervals for the alignment. If multiple alignments are found,
        then raise an error

        Args:
            input_sequence (str): the sequence to be aligned to
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by template/chr name

        Returns:
            :any:`list` of :any:`DomainRegion`: the list of domain regions on the new input sequence

        Raises:
            AttributeError: if sequence information is not available
            UserWarning: if a valid alignment could not be found or no best alignment was found
        """
        seq_list = self.get_sequences(REFERENCE_GENOME)

        results = {}
        for seq in set(seq_list):
            # align the current sequence to find the best matches
            scores = []
            for p in range(0, len(input_sequence) - len(seq)):
                score = 0
                for i in range(0, len(seq)):
                    if input_sequence[p + i].upper() == seq[i].upper():
                        score += 1
                if score > 0:
                    scores.append((Interval(p, p + len(seq) - 1), score))
            if len(scores) == 0:
                raise UserWarning('could not align a given region')
            results[seq] = scores

        # take the best score for each region and see if they work in sequence
        best = []
        for seq in seq_list:
            temp = max([t[1] for t in results[seq]])
            curr = [t for t in results[seq] if t[1] == temp]
            best.append(curr)

        combinations = [best[0]]
        for options in best[1:]:
            new_combinations = []
            for current_set in combinations:
                last, last_score = current_set[-1]
                for pos, score in options:
                    if pos.start > last.end:
                        new_combinations.append(current_set + [(pos, score)])
            combinations = new_combinations

        # now go through the list for the highest score
        scored_combinations = []
        for pos_list in combinations:
            score = sum([t[1] for t in pos_list])
            scored_combinations.append((score, [t[0] for t in pos_list]))

        if len(scored_combinations) == 0:
            raise UserWarning('could not map the sequences to the input')
        else:
            high = max([s for s, p in scored_combinations])
            temp = []
            for score, pl in scored_combinations:
                if score == high:
                    temp.append(pl)
            if len(temp) > 1:
                raise UserWarning('multiple mappings of equal score')
            else:
                regions = []
                for pos, seq in zip(temp[0], seq_list):
                    regions.append(DomainRegion(pos.start + 1, pos.end + 1, seq))
                return regions


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
        filepath (str): path to the input tab-delimited file
    Returns:
        :class:`dict` of :class:`str` and :class:`list` of :class:`BioInterval`:
            a dictionary keyed by chromosome name with values of lists of regions on the chromosome

    Example:
        >>> m = load_masking_regions('filename')
        >>> m['1']
        [BioInterval(), BioInterval(), ...]
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


def load_reference_genes(filepath, verbose=True):
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
        filepath (str): path to the input tab-delimited file

    Returns:
        :class:`dict` of :class:`str` and :class:`list` of :any:`Gene`: a dictionary keyed by chromosome name with values of list of genes on the chromosome

    Example:
        >>> ref = load_reference_genes('filename')
        >>> ref['1']
        [Gene(), Gene(), ....]
    """
    def parse_exon_list(row):
        if not row:
            return []
        exons = []
        for temp in re.split('[; ]', row):
            try:
                s, t = temp.split('-')
                exons.append((int(s), int(t)))
            except Exception as err:
                if verbose:
                    print('exon error:', repr(temp), repr(err))
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
                temp = Interval.min_nonoverlapping(*temp)
                d = Domain(name, temp)
                domains.append(d)
            except Exception as err:
                if verbose:
                    print('error in domain:', d, row, repr(err))
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
        try:
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
        except AttributeError as err:
            print('failed loading', row['ensembl_transcript_id'], row['hugo_names'].split(';'), repr(err))

    ref = {}
    for g in genes.values():
        if g.chr not in ref:
            ref[g.chr] = []
        ref[g.chr].append(g)
    return ref


def overlapping_transcripts(ref_ann, breakpoint):
    """
    Args:
        ref_ann (:class:`dict` of :class:`str` and :class:`list` of :any:`Gene`): the reference list of genes split by chromosome
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
        ref_ann (:class:`dict` of :class:`str` and :class:`list` of :class:`Gene`): the reference annotations split into lists of genes by chromosome
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
        ref (:class:`dict` of :class:`str` and :class:`list` of :any:`Gene`): the list of reference genes hashed by chromosomes
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


def load_reference_genome(filename):
    """
    Args:
        filename (str): the path to the file containing the input fasta genome

    Returns:
        :class:`dict` of :class:`str` and :class:`Bio.SeqRecord`: a dictionary representing the sequences in the fasta file
    """
    HUMAN_REFERENCE_GENOME = None
    with open(filename, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    return HUMAN_REFERENCE_GENOME
