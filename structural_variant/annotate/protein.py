from ..interval import Interval
from .base import BioInterval
from ..constants import translate, START_AA, STOP_AA, CODON_SIZE
import itertools
from ..error import NotSpecifiedError, DiscontinuousMappingError


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
                starts.append(i * CODON_SIZE + 1 + offset)
            elif aa == STOP_AA:
                stops.append((i + 1) * CODON_SIZE + offset)

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


class DomainRegion(BioInterval):
    def __init__(self, start, end, sequence=None, domain=None, name=None):
        BioInterval.__init__(self, domain, start, end, sequence=sequence, name=name)
        if sequence and len(sequence) != len(self):
            raise AssertionError('domain region sequence must be of equal length', self, sequence)


class Domain:
    """
    """
    def __init__(self, name, regions, translation=None, data=None):
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
        self.data = dict()
        if data is not None :
            self.data.update(data)
        if len(regions) == 0:
            raise AttributeError('at least one region must be given')
        for r1, r2 in itertools.combinations(self.regions, 2):
            if Interval.overlaps(r1, r2):
                raise AttributeError('regions cannot overlap')

        for i in range(0, len(self.regions)):
            curr = self.regions[i]
            if not hasattr(curr, 'sequence'):
                self.regions[i] = DomainRegion(curr[0], curr[1])

    @property
    def translation(self):
        """:class:`~structural_variant.annotate.Translation`: the Translation this domain belongs to"""
        return self.reference_object

    def key(self):
        """:class:`tuple`: a tuple representing the items expected to be unique. for hashing and comparing"""
        return tuple([self.name, self.translation])

    def score_region_mapping(self, REFERENCE_GENOME=None):
        """
        compares the sequence in each DomainRegion to the sequence collected for that domain region from the
        translation object

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

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
                    raise NotSpecifiedError('insufficient sequence information')
                ref = aa[region.start - 1:region.end]
                for c1, c2 in zip(ref, region.sequence):
                    if c1 == c2:
                        matches += 1
                    total += 1
            return matches, total
        else:
            raise NotSpecifiedError('insufficient sequence information')

    def get_sequences(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        returns the amino acid sequences for each of the domain regions associated with
        this domain in the order of the regions (sorted by start)

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            :class:`list` of :class:`str`: list of amino acid sequences for each DomainRegion

        Raises:
            AttributeError: if there is not enough sequence information given to determine this
        """
        sequences = {}
        if not ignore_cache:
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
                raise NotSpecifiedError('insufficient sequence information')
        return [sequences[r] for r in self.regions]

    def align_seq(self, input_sequence, REFERENCE_GENOME=None):
        """
        align each region to the input sequence starting with the last one.
        then take the subset of sequence that remains to align the second last and so on
        return a list of intervals for the alignment. If multiple alignments are found,
        then raise an error

        Args:
            input_sequence (str): the sequence to be aligned to
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            tuple: tuple contains

                - int: the number of matches
                - int: the total number of amino acids to be aligned
                - :any:`list` of :any:`DomainRegion`: the list of domain regions on the new input sequence

        Raises:
            AttributeError: if sequence information is not available
            UserWarning: if a valid alignment could not be found or no best alignment was found
        """
        seq_list = self.get_sequences(REFERENCE_GENOME)
        total = sum([len(s) for s in seq_list])

        results = {}
        for seq in set(seq_list):
            # align the current sequence to find the best matches
            scores = []
            for p in range(0, len(input_sequence) - len(seq) + 1):
                score = 0
                for i in range(0, len(seq)):
                    if input_sequence[p + i].upper() == seq[i].upper():
                        score += 1
                if score > 0:
                    scores.append((Interval(p + 1, p + len(seq)), score))
            if len(scores) == 0:
                raise UserWarning('could not align a given region')
            results.setdefault(seq, []).extend(scores)

        # take the best score for each region and see if they work in sequence
        best = []
        for seq in seq_list:
            temp = max([s for i, s in results[seq]])
            curr = [(i, s) for i, s in results[seq] if s == temp]
            best.append(curr)

        combinations = []
        for combo in itertools.product(*best):
            total_score = sum([s for i, s in combo])
            valid = True
            for i in range(1, len(combo)):
                if combo[i][0].start <= combo[i - 1][0].end:
                    valid = False
                    break
            if valid:
                combinations.append((total_score, [i for i, s in combo]))

        if len(combinations) == 0:
            raise UserWarning('could not map the sequences to the input')
        else:
            high = max([s for s, pl in combinations])
            temp = []
            for score, pl in combinations:
                if score == high:
                    temp.append(pl)
            if len(temp) > 1:
                raise UserWarning('multiple mappings of equal score')
            else:
                regions = []
                for pos, seq in zip(temp[0], seq_list):
                    regions.append(DomainRegion(pos.start, pos.end, seq))
                return high, total, regions


class Translation(BioInterval):
    def __init__(self, start, end, transcript=None, domains=None, sequence=None, name=None):
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
        BioInterval.__init__(self, reference_object=transcript, name=name, start=start, end=end, sequence=sequence)
        self.domains = [d for d in domains]

        if start <= 0:
            raise AttributeError('start must be a positive integer')
        if transcript and end > len(transcript):
            raise AttributeError('translation cannot be outside of related transcript range', end, len(transcript))

        for d in domains:
            d.reference_object = self

    @property
    def transcript(self):
        """:class:`~structural_variant.annotate.genomic.Transcript`: the spliced transcript this translation belongs to"""
        return self.reference_object

    def convert_aa_to_cdna(self, pos):
        """
        Args:
            pos (int): the amino acid position

        Returns:
            Interval: the cdna equivalent (with CODON_SIZE uncertainty)
        """
        return Interval(self.start - 1 + (pos - 1) * 3 + 1, self.start - 1 + pos * 3)

    def convert_cdna_to_aa(self, pos):
        """
        Args:
            pos (int): the cdna position

        Returns:
            Interval: the cdna equivalent (with CODON_SIZE uncertainty)

        Raises:
            AttributeError: the cdna position is not translated
        """
        if pos < self.start or pos > self.end:
            raise IndexError('position is out of bounds')
        pos = pos - self.start + 1
        aa = pos // CODON_SIZE
        if pos % CODON_SIZE != 0:
            aa += 1
        return aa

    def convert_genomic_to_cds(self, pos):
        cdna_pos = self.transcript.convert_genomic_to_cdna(pos)
        if cdna_pos < self.start:
            return cdna_pos - self.start
        return cdna_pos - self.start + 1

    def convert_genomic_to_cds_notation(self, pos):
        try:
            cds_pos = self.convert_genomic_to_cds(pos)
            if cds_pos > len(self):
                return '*{}'.format(cds_pos - len(self))
            print(cds_pos, self)
            return '{}'.format(cds_pos)
        except DiscontinuousMappingError as err:  # should give you the nearest positions
            # between two exons?
            exon_list = self.transcript.exons
            for ex1, ex2 in zip(self.transcript.exons[0::], self.transcript.exons[1::]):
                if abs(Interval.dist(ex1, ex2)) > 0:
                    intron = Interval(ex1.end + 1, ex2.start - 1)
                    if pos >= intron.start and pos <= intron.end:
                        # inside this intron
                        if abs(pos - intron.end) > abs(pos - intron.start):  # prefer +
                            ref_pos = self.convert_genomic_to_cds(ex1.end)
                            shift = pos - intron.start + 1
                            return '{}+{}'.format(ref_pos, shift)
                        else:
                            ref_pos = self.convert_genomic_to_cds(ex2.start)
                            shift = intron.end - pos + 1
                            return '{}-{}'.format(ref_pos, shift)
            raise err


    def get_cds_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            str: the cds sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        if self.sequence and not ignore_cache:
            return self.sequence
        elif self.transcript and self.transcript.get_strand():
            seq = self.transcript.get_sequence(REFERENCE_GENOME, ignore_cache)
            return seq[self.start - 1:self.end]
        raise NotSpecifiedError('insufficient sequence information')

    def get_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        wrapper for the sequence method

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name
        """
        return self.get_cds_sequence(REFERENCE_GENOME, ignore_cache)

    def get_AA_sequence(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            str: the amino acid sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        cds = self.get_cds_sequence(REFERENCE_GENOME, ignore_cache)
        return translate(cds)

    def key(self):
        """see :func:`structural_variant.annotate.base.BioInterval.key`"""
        return BioInterval.key(self), self.splicing_pattern
