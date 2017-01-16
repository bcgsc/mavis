from ..constants import STRAND, SPLICE_SITE_RADIUS, reverse_complement, CODON_SIZE
from ..interval import Interval
from ..error import StrandSpecificityError
from .protein import Translation
from .base import BioInterval
import warnings


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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by
                template/chr name

        Returns:
            str: the sequence of the gene
        """
        if self.sequence:
            return self.sequence
        elif REFERENCE_GENOME is None:
            raise AttributeError('reference genome is required to retrieve the gene sequence')
        else:
            return str(REFERENCE_GENOME[self.chr].seq[self.start - 1:self.end]).upper()


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
            :class:`~structural_variant.error.DiscontinuousMappingError`: when a genomic position not present in the
                cdna is attempted to be converted
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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

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
        seq = self.get_sequence(REFERENCE_GENOME)
        if self.strand == STRAND.NEG:
            # adjust the continuous intervals for the min and flip if revcomp
            seq = reverse_complement(seq)
        spliced_seq = ''.join([str(seq[i.start:i.end + 1]) for i in conti])
        return spliced_seq if self.strand == STRAND.POS else reverse_complement(spliced_seq)
