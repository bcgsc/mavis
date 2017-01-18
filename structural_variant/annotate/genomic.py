from ..constants import STRAND, SPLICE_SITE_RADIUS, reverse_complement
from ..interval import Interval
from ..error import StrandSpecificityError
from .base import BioInterval
import warnings
import itertools


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

    def key(self):
        return BioInterval.key(self), self.strand
    
    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object


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
        BioInterval.__init__(self, name=name, reference_object=chr, start=start, end=end, sequence=sequence)
        self.unspliced_transcripts = []
        self.strand = STRAND.enforce(strand)
        self.aliases = aliases
    
    @property
    def transcripts(self):
        return self.unspliced_transcripts

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    def key(self):
        return BioInterval.key(self), self.strand

    def get_sequence(self, REFERENCE_GENOME, ignore_cache=False):
        """
        gene sequence is always given wrt to the positive forward strand regardless of gene strand

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence by
                template/chr name

        Returns:
            str: the sequence of the gene
        """
        if self.sequence and not ignore_cache:
            return self.sequence
        elif REFERENCE_GENOME is None:
            raise AttributeError('reference genome is required to retrieve the gene sequence')
        else:
            return str(REFERENCE_GENOME[self.chr].seq[self.start - 1:self.end]).upper()
    
    @property
    def spliced_transcripts(self):
        spl = []
        for t in self.unspliced_transcripts:
            spl.extend(t.spliced_transcripts)
        return spl


class Exon(BioInterval):
    """
    """
    def __init__(
            self, start, end,
            transcript=None,
            name=None,
            intact_start_splice=True,
            intact_end_splice=True,
            sequence=None):
        """
        Args:
            start (int): the genomic start position
            end (int): the genomic end position
            name (str): the name of the exon
            transcript (usTranscript): the 'parent' transcript this exon belongs to
            intact_start_splice (bool): if the starting splice site has been abrogated
            intact_end_splice (bool): if the end splice site has been abrogated
        Raises:
            AttributeError: if the exon start > the exon end
        Example:
            >>> Exon(15, 78)
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end, sequence=sequence)
        self.intact_start_splice = intact_start_splice
        self.intact_end_splice = intact_end_splice
        if end - start + 1 < SPLICE_SITE_RADIUS * sum([intact_start_splice, intact_end_splice]):
            warnings.warn('exons must be greater than double the length of a splice site')

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.usTranscript`): the transcript this exon belongs to"""
        return self.reference_object

    @property
    def start_splice_site(self):
        """(:class:`~structural_variant.interval.Interval`): the genomic range describing the splice site"""
        return Interval(self.start - SPLICE_SITE_RADIUS, self.start + SPLICE_SITE_RADIUS - 1)

    @property
    def end_splice_site(self):
        """(:class:`~structural_variant.interval.Interval`): the genomic range describing the splice site"""
        return Interval(self.end - SPLICE_SITE_RADIUS + 1, self.end + SPLICE_SITE_RADIUS)


class usTranscript(BioInterval):
    """
    """
    def __init__(
        self,
        exons=None,
        start=None,
        end=None,
        gene=None,
        name=None,
        strand=None,
        spliced_transcripts=None,
        sequence=None
    ):
        """ creates a new transcript object

        Args:
            exons (:class:`list` of :any:`Exon`): list of Exon that make up the transcript
            genomic_start (int): genomic start position of the transcript
            genomic_end (int): genomic end position of the transcript
            gene (Gene): the gene this transcript belongs to
            name (str): name of the transcript
            strand (STRAND): strand the transcript is on, defaults to the strand of the Gene if not specified
            sequence (str): unspliced cDNA sequence
        """
        # cannot use mutable default args in the function decl
        self.exons = [] if exons is None else exons
        self.spliced_transcripts = [] if spliced_transcripts is None else spliced_transcripts
        self.strand = strand

        if start is None and len(self.exons) > 0:
            start = min([e[0] for e in self.exons])
        if end is None and len(self.exons) > 0:
            end = max([e[1] for e in self.exons])

        BioInterval.__init__(self, gene, start, end, name=name, sequence=sequence)
        
        for i in range(0, len(self.exons)):
            curr = self.exons[i]
            try:
                curr.reference_object = self
            except AttributeError:
                self.exons[i] = Exon(curr[0], curr[1], self)

        for e1, e2 in itertools.combinations(self.exons, 2):
            if Interval.overlaps(e1, e2):
                raise AttributeError('exons cannot overlap')

        for s in self.spliced_transcripts:
            s.reference_object = self
        
        try:
            if self.get_strand() != self.gene.get_strand():
                raise AssertionError('gene strand and transcript strand conflict')
        except AttributeError:
            pass

    def generate_splicing_patterns(self):
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

    @property
    def gene(self):
        """(:any:`Gene`): the gene this transcript belongs to"""
        return self.reference_object

    def _genomic_to_cdna_mapping(self, splicing_pattern):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating
        """
        mapping = {}

        exons = sorted(self.exons, key=lambda x: x.start)
        if self.get_strand() == STRAND.POS:
            pass
        elif self.get_strand() == STRAND.NEG:
            exons.reverse()
        else:
            raise StrandSpecificityError('cannot convert without strand information')

        l = 1
        for e in exons:
            mapping[Interval(e.start, e.end)] = Interval(l, l + len(e) - 1)
            l += len(e)
        return mapping

    def _cdna_to_genomic_mapping(self, splicing_pattern):
        """
        Args:
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating
        """
        mapping = {}
        for k, v in self._genomic_to_cdna_mapping(splicing_pattern).items():
            mapping[v] = k
        return mapping

    def convert_genomic_to_cdna(self, pos, splicing_pattern):
        """
        Args:
            pos (int): the genomic position to be converted
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating

        Returns:
            int: the cdna equivalent

        Raises:
            :class:`~structural_variant.error.DiscontinuousMappingError`: when a genomic position not present in the
                cdna is attempted to be converted
        """
        mapping = self._genomic_to_cdna_mapping(splicing_pattern)
        return Interval.convert_pos(mapping, pos)

    def convert_cdna_to_genomic(self, pos, splicing_pattern):
        """
        Args:
            pos (int): cdna position
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating

        Returns:
            int: the genomic equivalent
        """
        mapping = self._cdna_to_genomic_mapping(splicing_pattern)
        return Interval.convert_pos(mapping, pos)

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
            if self.get_strand() == STRAND.POS:
                return i + 1
            elif self.get_strand() == STRAND.NEG:
                return len(self.exons) - i
            else:
                raise AttributeError('strand must be pos or neg to calculate the exon number')
        raise AttributeError('can only calculate phase on associated exons')

    def get_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

        Returns:
            str: the sequence of the transcript including introns (but relative to strand)
        """
        if self.sequence and not ignore_cache:
            return self.sequence
        elif self.gene and self.gene.sequence and not ignore_cache:
            # gene has a sequence set
            start = self.start - self.gene.start
            end = self.end - self.gene.end + len(self.gene.sequence)
            if self.get_strand() == STRAND.NEG:
                return reverse_complement(self.gene.sequence[start:end])
            else:
                return self.gene.sequence[start:end]
        elif REFERENCE_GENOME is None:
            raise AttributeError('reference genome is required to retrieve the gene sequence')
        else:
            if self.get_strand() == STRAND.NEG:
                return reverse_complement(REFERENCE_GENOME[self.gene.chr].seq[self.start - 1:self.end]).upper()
            else:
                return str(REFERENCE_GENOME[self.gene.chr].seq[self.start - 1:self.end]).upper()

    def get_cdna_sequence(self, splicing_pattern, REFERENCE_GENOME=None, ignore_cache=False):
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
        seq = self.get_sequence(REFERENCE_GENOME, ignore_cache)
        if self.get_strand() == STRAND.NEG:
            # adjust the continuous intervals for the min and flip if revcomp
            seq = reverse_complement(seq)
        spliced_seq = ''.join([str(seq[i.start:i.end + 1]) for i in conti])
        spliced_seq = spliced_seq.upper()
        return spliced_seq if self.get_strand() == STRAND.POS else reverse_complement(spliced_seq)
    
    @property
    def translations(self):
        tx = []
        for t in self.spliced_transcripts:
            for tl in t.translations:
                tx.append(tx)
        return tx
    
    @property
    def transcripts(self):
        return self.spliced_transcripts


class Transcript(BioInterval):
    def __init__(self, ust, splicing_patt, sequence=None, translations=None):
        """
        splicing pattern is given in genomic coordinates

        Args:
            us_transcript (usTranscript): the unspliced transcript
            splicing_patt (:class:`list` of :class:`int`): the list of splicing positions
            sequence (str): the cdna sequence
            translations (:class:`list` of :class:`Translation`): the list of translations of this transcript
        """
        pos = sorted([ust.start, ust.end] + splicing_patt)
        self.splicing_pattern = sorted(splicing_patt)
        self.exons = [Exon(s, t, self) for s, t in zip(pos[::2], pos[1::2])]
        
        BioInterval.__init__(self, ust, 1, sum([len(e) for e in self.exons]), sequence=None)
        self.translations = [] if translations is None else [tx for tx in translations]
        
        for tx in self.translations:
            tx.reference_object = self
        if len(splicing_patt) > 0 and (min(splicing_patt) < ust.start or max(splicing_patt) > ust.end):
            raise AssertionError('splicing pattern must be contained within the unspliced transcript')
        elif len(splicing_patt) % 2 != 0:
            raise AssertionError('splicing pattern must be a list of 3\'5\' splicing positions')
    
    def convert_genomic_to_cdna(self, pos):
        """
        Args:
            pos (int): the genomic position to be converted
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating

        Returns:
            int: the cdna equivalent

        Raises:
            :class:`~structural_variant.error.DiscontinuousMappingError`: when a genomic position not present in the
                cdna is attempted to be converted
        """
        mapping = self.unspliced_transcript._genomic_to_cdna_mapping(self.splicing_pattern)
        return Interval.convert_pos(mapping, pos)

    def convert_cdna_to_genomic(self, pos):
        """
        Args:
            pos (int): cdna position
            splicing_pattern (:class:`list` of :class:`int`): list of genomic splice sites 3'5' repeating

        Returns:
            int: the genomic equivalent
        """
        mapping = self.unspliced_transcript._cdna_to_genomic_mapping(self.splicing_pattern)
        return Interval.convert_pos(mapping, pos)
    
    def get_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        if self.sequence and not ignore_cache:
            return self.sequence
        else:
            s = self.unspliced_transcript.get_cdna_sequence(self.splicing_pattern, REFERENCE_GENOME, ignore_cache)
            return s[self.start - 1:self.end]
    
    @property
    def unspliced_transcript(self):
        return self.reference_object

