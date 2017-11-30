from copy import copy
import itertools

from .base import BioInterval, ReferenceName
from .constants import SPLICE_SITE_TYPE
from .splicing import SpliceSite, SplicingPattern
from ..constants import ORIENT, reverse_complement, STRAND
from ..error import NotSpecifiedError
from ..interval import Interval


class Template(BioInterval):

    def __init__(self, name, start, end, seq=None, bands=None):
        bands = [] if bands is None else bands
        name = ReferenceName(name)
        BioInterval.__init__(self, None, start, end, name=name, seq=seq)
        self.bands = bands
        for i in range(0, len(bands)):
            bands[i].reference_object = self
        self.bands.sort()

    def __str__(self):
        return str(self.name)

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(self.name)


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
        """see :func:`structural_variant.annotate.base.BioInterval.key`"""
        return BioInterval.key(self), self.strand

    @property
    def chr(self):
        """returns the name of the chromosome that this region resides on"""
        return self.reference_object

    def __repr__(self):
        return 'IntergenicRegion({}:{}_{}{})'.format(self.chr, self.start, self.end, self.strand)

    def to_dict(self):
        """see :func:`structural_variant.annotate.base.BioInterval.to_dict`"""
        d = BioInterval.to_dict(self)
        d['strand'] = self.strand
        return d


class Gene(BioInterval):
    """
    """

    def __init__(self, chr, start, end, name=None, strand=STRAND.NS, aliases=None, seq=None):
        """
        Args:
            chr (str): the chromosome
            name (str): the gene name/id i.e. ENSG0001
            strand (STRAND): the genomic strand '+' or '-'
            aliases (:class:`list` of :class:`str`): a list of aliases. For example the hugo name could go here
            seq (str): genomic seq of the gene
        Example:
            >>> Gene('X', 1, 1000, 'ENG0001', '+', ['KRAS'])
        """
        aliases = [] if aliases is None else aliases
        chr = ReferenceName(chr)
        BioInterval.__init__(self, name=name, reference_object=chr, start=start, end=end, seq=seq)
        self.unspliced_transcripts = []
        self.strand = STRAND.enforce(strand)
        self.aliases = aliases

    def transcript_priority(self, transcript):
        """
        prioritizes transcripts from 0 to n-1 based on best transcript flag
        and then alphanumeric name sort

        Warning:
            Lower number means higher priority. This is to make sort work by default
        """
        def sort_key(t):
            return (
                0 if t.is_best_transcript else 1,
                t.name, t.start - t.end, t.start, t.end
            )
        priority = sorted(self.transcripts, key=sort_key)
        for i, curr_transcript in enumerate(priority):
            if curr_transcript == transcript:
                return i
        raise ValueError('input transcript is not associated with this gene', transcript)

    @property
    def transcripts(self):
        """:any:`list` of :class:`PreTranscript`: list of unspliced transcripts"""
        return self.unspliced_transcripts

    @property
    def translations(self):
        """:any:`list` of :class:`~mavis.annotate.protein.Translation`: list of translations"""
        translations = []
        for pre_transcript in self.unspliced_transcripts:
            for tx in pre_transcript.transcripts:
                for tl in tx.translations:
                    translations.append(tl)
        return translations

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    def key(self):
        """see :func:`structural_variant.annotate.base.BioInterval.key`"""
        return BioInterval.key(self), self.strand

    def get_seq(self, reference_genome, ignore_cache=False):
        """
        gene sequence is always given wrt to the positive forward strand regardless of gene strand

        Args:
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by
                template/chr name
            ignore_cache (bool): if True then stored sequences will be ignored and the function will attempt to retrieve the sequence using the positions and the input reference_genome

        Returns:
            str: the sequence of the gene
        """
        if self.seq and not ignore_cache:
            return self.seq
        elif reference_genome is None:
            raise NotSpecifiedError('reference genome is required to retrieve the gene sequence')
        else:
            return str(reference_genome[self.chr].seq[self.start - 1:self.end]).upper()

    @property
    def spliced_transcripts(self):
        """:any:`list` of :class:`Transcript`: list of transcripts"""
        spl = []
        for t in self.unspliced_transcripts:
            spl.extend(t.spliced_transcripts)
        return spl

    def to_dict(self):
        """see :func:`structural_variant.annotate.base.BioInterval.to_dict`"""
        d = BioInterval.to_dict(self)
        d['strand'] = self.strand
        return d


class Exon(BioInterval):
    """
    """

    def __init__(
            self, start, end,
            transcript=None,
            name=None,
            intact_start_splice=True,
            intact_end_splice=True,
            seq=None,
            strand=None):
        """
        Args:
            start (int): the genomic start position
            end (int): the genomic end position
            name (str): the name of the exon
            transcript (PreTranscript): the 'parent' transcript this exon belongs to
            intact_start_splice (bool): if the starting splice site has been abrogated
            intact_end_splice (bool): if the end splice site has been abrogated
        Raises:
            AttributeError: if the exon start > the exon end
        Example:
            >>> Exon(15, 78)
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end, seq=seq, strand=strand)

        if self.is_reverse:
            self.start_splice_site = SpliceSite(
                self.transcript, self.start, site_type=SPLICE_SITE_TYPE.DONOR, strand=STRAND.NEG, intact=intact_start_splice)
            self.end_splice_site = SpliceSite(
                self.transcript, self.end, site_type=SPLICE_SITE_TYPE.ACCEPTOR, strand=STRAND.NEG, intact=intact_end_splice)
        else:
            self.start_splice_site = SpliceSite(
                self.transcript, self.start, site_type=SPLICE_SITE_TYPE.ACCEPTOR, strand=STRAND.POS, intact=intact_start_splice)
            self.end_splice_site = SpliceSite(
                self.transcript, self.end, site_type=SPLICE_SITE_TYPE.DONOR, strand=STRAND.POS, intact=intact_end_splice)

    @property
    def transcript(self):
        """:class:`PreTranscript`: the transcript this exon belongs to"""
        return self.reference_object

    @property
    def donor_splice_site(self):
        """:class:`~mavis.interval.Interval`: the genomic range describing the splice site"""
        if self.is_reverse:
            return self.start_splice_site
        else:
            return self.end_splice_site

    @property
    def acceptor_splice_site(self):
        """:class:`~mavis.interval.Interval`: the genomic range describing the splice site"""
        if self.is_reverse:
            return self.end_splice_site
        else:
            return self.start_splice_site

    @property
    def donor(self):
        """`int`: returns the genomic exonic position of the donor splice site"""
        if self.is_reverse:
            return self.start
        else:
            return self.end

    @property
    def acceptor(self):
        """`int`: returns the genomic exonic position of the acceptor splice site"""
        if self.is_reverse:
            return self.end
        else:
            return self.start

    def __repr__(self):
        return 'Exon({}{}, {}{})'.format(
            self.start, '' if self.start_splice_site.intact else '*',
            self.end, '' if self.end_splice_site.intact else '*')


class PreTranscript(BioInterval):
    """
    """

    def __init__(
        self,
        exons,
        gene=None,
        name=None,
        strand=None,
        spliced_transcripts=None,
        seq=None,
        is_best_transcript=False
    ):
        """ creates a new transcript object

        Args:
            exons (:class:`list` of :any:`Exon`): list of Exon that make up the transcript
            genomic_start (int): genomic start position of the transcript
            genomic_end (int): genomic end position of the transcript
            gene (Gene): the gene this transcript belongs to
            name (str): name of the transcript
            strand (STRAND): strand the transcript is on, defaults to the strand of the Gene if not specified
            seq (str): unspliced cDNA seq
        """
        # cannot use mutable default args in the function decl
        self.exons = exons
        self.spliced_transcripts = [] if spliced_transcripts is None else spliced_transcripts
        self.is_best_transcript = is_best_transcript

        if len(exons) == 0:
            raise AttributeError('exons must be given')

        start = min([e[0] for e in self.exons])
        end = max([e[1] for e in self.exons])

        BioInterval.__init__(self, gene, start, end, name=name, seq=seq, strand=strand)

        for i, curr in enumerate(self.exons):
            if isinstance(curr, Exon):
                curr.reference_object = self
            else:
                self.exons[i] = Exon(curr[0], curr[1], self)
        self.exons = sorted(self.exons, key=lambda x: x.start)
        for ex in self.exons:
            if ex.end > self.end or ex.start < self.start:
                raise AssertionError('exon is outside transcript', self, ex)
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
            :class:`list` of :class:`SplicingPattern`: List of positions to be spliced together

        see :ref:`theory - predicting splicing patterns <theory-predicting-splicing-patterns>`
        """
        if len(self.exons) < 2:
            return [SplicingPattern()]
        sites = [self.exons[0].end_splice_site]
        for exon in self.exons[1:-1]:
            sites.extend([exon.start_splice_site, exon.end_splice_site])
        sites.append(self.exons[-1].start_splice_site)
        return SplicingPattern.generate_patterns(sites, is_reverse=self.is_reverse)

    @property
    def gene(self):
        """:any:`Gene`: the gene this transcript belongs to"""
        return self.reference_object

    def _genomic_to_cdna_mapping(self, splicing_pattern):
        """
        Args:
            splicing_pattern (SplicingPattern): list of genomic splice sites 3'5' repeating
        """
        mapping = {}
        length = 1
        pos = sorted([s.pos for s in splicing_pattern] + [self.start, self.end])
        genome_intervals = [Interval(s, t) for s, t in zip(pos[::2], pos[1::2])]

        if self.get_strand() == STRAND.POS:
            pass
        elif self.get_strand() == STRAND.NEG:
            genome_intervals.reverse()
        else:
            raise NotSpecifiedError('cannot convert without strand information')

        for exon in genome_intervals:
            mapping[exon] = Interval(length, length + len(exon) - 1)
            length += len(exon)
        return mapping

    def _cdna_to_genomic_mapping(self, splicing_pattern):
        """
        Args:
            splicing_pattern (SplicingPattern): list of genomic splice sites 3'5' repeating
        """
        mapping = {v: k for k, v in self._genomic_to_cdna_mapping(splicing_pattern).items()}
        return mapping

    def convert_genomic_to_cdna(self, pos, splicing_pattern):
        """
        Args:
            pos (int): the genomic position to be converted
            splicing_pattern (SplicingPattern): list of genomic splice sites 3'5' repeating

        Returns:
            int: the cdna equivalent

        Raises:
            :class:`~mavis.error.IndexError`: when a genomic position not present in the
                cdna is attempted to be converted
        """
        cdna_pos, shift = self.convert_genomic_to_nearest_cdna(pos, splicing_pattern)
        if shift != 0:
            raise IndexError('outside of exonic regions', pos, splicing_pattern, cdna_pos, shift)
        return cdna_pos

    def convert_genomic_to_nearest_cdna(self, pos, splicing_pattern, stick_direction=None, allow_outside=True):
        """
        converts a genomic position to its cdna equivalent or (if intronic) the nearest cdna and shift

        Args:
            pos (int): the genomic position
            splicing_pattern (SplicingPattern): the splicing pattern

        Returns:
            tuple of int and int:
                * *int* - the exonic cdna position
                * *int* - the intronic shift

        """
        mapping = self._genomic_to_cdna_mapping(splicing_pattern)
        exons = sorted(list(mapping.keys()))
        # exonic
        for ex in exons:
            if pos <= ex.end and pos >= ex.start:
                # in the current exon
                cdna_pos = Interval.convert_pos(mapping, pos, True if self.get_strand() == STRAND.NEG else False)
                return cdna_pos, 0
        # intronic
        for ex1, ex2 in zip(exons, exons[1::]):
            if pos > ex1.end and pos < ex2.start:
                # in the current intron
                if (abs(pos - ex1.end) <= abs(pos - ex2.start) or stick_direction == ORIENT.LEFT) and stick_direction != ORIENT.RIGHT:
                    # closest to the first exon
                    cdna_pos = Interval.convert_pos(mapping, ex1.end, True if self.get_strand() == STRAND.NEG else False)
                    return cdna_pos, pos - ex1.end if self.get_strand() == STRAND.POS else ex1.end - pos
                cdna_pos = Interval.convert_pos(mapping, ex2.start, True if self.get_strand() == STRAND.NEG else False)
                return cdna_pos, pos - ex2.start if self.get_strand() == STRAND.POS else ex2.start - pos
        if allow_outside:
            cdna_length = sum([len(e) for e in exons])
            if pos < exons[0].start:  # before the first exon
                return cdna_length if self.is_reverse else 1, pos - exons[0].start
            elif pos > exons[-1].end:
                return 1 if self.is_reverse else cdna_length, pos - exons[-1].end
            else:
                raise NotImplementedError('Unexpected error', self.exons, pos)
        raise IndexError('position does not fall within the current transcript', pos, mapping)

    def convert_cdna_to_genomic(self, pos, splicing_pattern):
        """
        Args:
            pos (int): cdna position
            splicing_pattern (SplicingPattern): list of genomic splice sites 3'5' repeating

        Returns:
            int: the genomic equivalent
        """
        mapping = self._cdna_to_genomic_mapping(splicing_pattern)
        exons = sorted(mapping.values())
        length = sum([len(e) for e in mapping])
        if pos < 0:
            if self.is_reverse:
                return exons[-1].end + abs(pos)
            return exons[0].start + pos
        if pos > length:
            pos -= length
            if self.is_reverse:
                return exons[0].start - pos
            return exons[-1].end + pos
        return Interval.convert_pos(mapping, pos, True if self.get_strand() == STRAND.NEG else False)

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
        for i, current_exon in enumerate(self.exons):
            if exon != current_exon:
                continue
            if self.get_strand() == STRAND.POS:
                return i + 1
            elif self.get_strand() == STRAND.NEG:
                return len(self.exons) - i
            else:
                raise NotSpecifiedError('strand must be pos or neg to calculate the exon number')
        raise AttributeError('can only calculate phase on associated exons')

    def get_seq(self, reference_genome=None, ignore_cache=False):
        """
        Args:
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name
            ignore_cache (bool): if True then stored sequences will be ignored and the function will attempt to retrieve the sequence using the positions and the input reference_genome

        Returns:
            str: the sequence of the transcript including introns (but relative to strand)
        """
        if self.seq and not ignore_cache:
            return self.seq
        elif self.gene and self.gene.seq and not ignore_cache:
            # gene has a seq set
            start = self.start - self.gene.start
            end = self.end - self.gene.end + len(self.gene.seq)
            if self.get_strand() == STRAND.NEG:
                return reverse_complement(self.gene.seq[start:end])
            return self.gene.seq[start:end]
        elif reference_genome is None:
            raise NotSpecifiedError('reference genome is required to retrieve the gene sequence')
        if self.get_strand() == STRAND.NEG:
            return reverse_complement(reference_genome[self.gene.chr].seq[self.start - 1:self.end]).upper()
        return str(reference_genome[self.gene.chr].seq[self.start - 1:self.end]).upper()

    def get_cdna_seq(self, splicing_pattern, reference_genome=None, ignore_cache=False):
        """
        Args:
            splicing_pattern (SplicingPattern): the list of splicing positions
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name
            ignore_cache (bool): if True then stored sequences will be ignored and the function will attempt to retrieve the sequence using the positions and the input reference_genome

        Returns:
            str: the spliced cDNA sequence
        """
        temp = sorted([self.start] + [s.pos for s in splicing_pattern] + [self.end])
        cdna_start = min(temp)
        conti = []
        for i in range(0, len(temp) - 1, 2):
            conti.append(Interval(temp[i] - cdna_start, temp[i + 1] - cdna_start))
        seq = self.get_seq(reference_genome, ignore_cache)
        if self.get_strand() == STRAND.NEG:
            # adjust the continuous intervals for the min and flip if revcomp
            seq = reverse_complement(seq)
        spliced_seq = ''.join([str(seq[i.start:i.end + 1]) for i in conti])
        spliced_seq = spliced_seq.upper()
        return spliced_seq if self.get_strand() == STRAND.POS else reverse_complement(spliced_seq)

    @property
    def translations(self):
        """:class:`list` of :class:`~mavis.annotate.protein.Translation`: list of translations associated with this transcript"""
        translations = []
        for spl_tx in self.spliced_transcripts:
            for translation in spl_tx.translations:
                translations.append(translation)
        return translations

    @property
    def transcripts(self):
        """:class:`list` of :class:`Transcript`: list of spliced transcripts"""
        return self.spliced_transcripts


class Transcript(BioInterval):

    def __init__(self, pre_transcript, splicing_patt, seq=None, translations=None):
        """
        splicing pattern is given in genomic coordinates

        Args:
            pre_transcript (PreTranscript): the unspliced transcript
            splicing_patt (:class:`list` of :class:`int`): the list of splicing positions
            seq (str): the cdna sequence
            translations (:class:`list` of :class:`~mavis.annotate.protein.Translation`):
             the list of translations of this transcript
        """
        pos = sorted([pre_transcript.start, pre_transcript.end] + [s.pos for s in splicing_patt])
        splicing_patt.sort()
        self.splicing_pattern = splicing_patt
        length = sum([t - s + 1 for s, t in zip(pos[::2], pos[1::2])])
        BioInterval.__init__(self, pre_transcript, 1, length, seq=None)
        self.exons = [Exon(s, t, self) for s, t in zip(pos[::2], pos[1::2])]
        self.translations = [] if translations is None else [tx for tx in translations]

        for translation in self.translations:
            translation.reference_object = self
        if splicing_patt and (min(splicing_patt).pos < pre_transcript.start or max(splicing_patt).pos > pre_transcript.end):
            raise AssertionError('splicing pattern must be contained within the unspliced transcript')
        elif len(splicing_patt) % 2 != 0:
            raise AssertionError('splicing pattern must be a list of 3\'5\' splicing positions')

    def convert_genomic_to_cdna(self, pos):
        """
        Args:
            pos (int): the genomic position to be converted

        Returns:
            int: the cdna equivalent

        Raises:
            IndexError: when a genomic position not present in the cdna is attempted to be converted
        """
        return self.unspliced_transcript.convert_genomic_to_cdna(pos, self.splicing_pattern)

    def convert_genomic_to_nearest_cdna(self, pos, **kwargs):
        return self.reference_object.convert_genomic_to_nearest_cdna(pos, self.splicing_pattern, **kwargs)

    def convert_cdna_to_genomic(self, pos):
        """
        Args:
            pos (int): cdna position

        Returns:
            int: the genomic equivalent
        """
        return self.unspliced_transcript.convert_cdna_to_genomic(pos, self.splicing_pattern)

    def get_seq(self, reference_genome=None, ignore_cache=False):
        """
        Args:
            reference_genome (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence by
                template/chr name
            ignore_cache (bool): if True then stored sequences will be ignored and the function will attempt to retrieve the sequence using the positions and the input reference_genome

        Returns:
            str: the sequence corresponding to the spliced cdna
        """
        if self.seq and not ignore_cache:
            return self.seq
        seq = self.unspliced_transcript.get_cdna_seq(self.splicing_pattern, reference_genome, ignore_cache)
        return seq[self.start - 1:self.end]

    @property
    def unspliced_transcript(self):
        """:class:`PreTranscript`: the unspliced transcript this splice variant belongs to"""
        return self.reference_object
