from ..interval import Interval
from .base import BioInterval
from ..constants import STRAND, translate, START_AA, STOP_AA, CODON_SIZE
import itertools


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
        compares the sequence in each DomainRegion to the sequence collected for that domain region from the
        translation object

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

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
                    'Insufficient information to gather sequence data. reference_object must have sequence'
                    ' information or REFERENCE_GENOME must be provided')
        return [sequences[r] for r in self.regions]

    def align_seq(self, input_sequence, REFERENCE_GENOME=None):
        """
        align each region to the input sequence starting with the last one.
        then take the subset of sequence that remains to align the second last and so on
        return a list of intervals for the alignment. If multiple alignments are found,
        then raise an error

        Args:
            input_sequence (str): the sequence to be aligned to
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
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
            for p in range(0, len(input_sequence) - len(seq)):
                score = 0
                for i in range(0, len(seq)):
                    if input_sequence[p + i].upper() == seq[i].upper():
                        score += 1
                if score > 0:
                    scores.append((Interval(p, p + len(seq) - 1), score))
            if len(scores) == 0:
                raise UserWarning('could not align a given region')
            results.setdefault(seq, []).extend(scores)

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
                return high, total, regions


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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

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
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name
        """
        return self.get_sequence(REFERENCE_GENOME)

    def get_AA_sequence(self, REFERENCE_GENOME=None):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`str` and :class:`Bio.SeqRecord`): dict of reference sequence
                by template/chr name

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
