from ..interval import Interval
from .base import BioInterval
from ..constants import translate, START_AA, STOP_AA, CODON_SIZE
import itertools
from ..error import NotSpecifiedError


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
        orf_intervals = []
        current_start = None
        for i, aa in enumerate(aa_sequence):
            if aa == START_AA:
                p = i * CODON_SIZE + 1 + offset
                if current_start is None:
                    current_start = p
            elif aa == STOP_AA:
                p = (i + 1) * CODON_SIZE + offset
                if current_start is not None:  # close the current interval
                    itvl = Interval(current_start, p)
                    if min_orf_size is None or len(itvl) >= min_orf_size:
                        orf_intervals.append(itvl)
                    current_start = None
        cds_orfs.extend(orf_intervals)
    return cds_orfs


class DomainRegion(BioInterval):
    def __init__(self, start, end, seq=None, domain=None, name=None):
        BioInterval.__init__(self, domain, start, end, seq=seq, name=name)
        if seq and len(seq) != len(self):
            raise AssertionError('domain region seq must be of equal length', self, seq)


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
        if data is not None:
            self.data.update(data)
        if len(regions) == 0:
            raise AttributeError('at least one region must be given')
        for r1, r2 in itertools.combinations(self.regions, 2):
            if Interval.overlaps(r1, r2):
                raise AttributeError('regions cannot overlap')

        for i in range(0, len(self.regions)):
            curr = self.regions[i]
            if not hasattr(curr, 'seq'):
                self.regions[i] = DomainRegion(curr[0], curr[1])

    @property
    def translation(self):
        """:class:`~mavis.annotate.Translation`: the Translation this domain belongs to"""
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
            aa = self.translation.get_AA_seq(REFERENCE_GENOME)
            total = 0
            matches = 0
            for region in self.regions:
                if not region.seq:
                    raise NotSpecifiedError('insufficient sequence information')
                ref = aa[region.start - 1:region.end]
                for c1, c2 in zip(ref, region.seq):
                    if c1 == c2:
                        matches += 1
                    total += 1
            return matches, total
        else:
            raise NotSpecifiedError('insufficient sequence information')

    def get_seqs(self, REFERENCE_GENOME=None, ignore_cache=False):
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
                if region.seq:
                    sequences[region] = region.seq
        if any([r not in sequences for r in self.regions]):
            if self.translation:
                aa = self.translation.get_AA_seq(REFERENCE_GENOME)
                for region in self.regions:
                    s = aa[region.start - 1:region.end]
                    if region not in sequences:
                        sequences[region] = s
            else:
                raise NotSpecifiedError('insufficient sequence information')
        return [sequences[r] for r in self.regions]

    def align_seq(self, input_sequence, REFERENCE_GENOME=None, min_region_match=0.5):
        """
        align each region to the input sequence starting with the last one.
        then take the subset of sequence that remains to align the second last and so on
        return a list of intervals for the alignment. If multiple alignments are found,
        then raise an error

        Args:
            input_sequence (str): the sequence to be aligned to
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name
            min_region_match (float): percent between 0 and 1. Each region must have a score len(seq) * min_region_match

        Returns:
            tuple: tuple contains

                - int: the number of matches
                - int: the total number of amino acids to be aligned
                - :any:`list` of :any:`DomainRegion`: the list of domain regions on the new input sequence

        Raises:
            AttributeError: if sequence information is not available
            UserWarning: if a valid alignment could not be found or no best alignment was found
        """
        seq_list = self.get_seqs(REFERENCE_GENOME)

        dr_by_seq = {s: d for s, d in zip(seq_list, self.regions)}
        seq_list = sorted(seq_list, key=lambda x: dr_by_seq[x].start)
        total = sum([len(s) for s in seq_list])

        if total > len(input_sequence):
            raise UserWarning('could not map the sequences to the input')

        results = []
        last_min_end = 0
        for seq in seq_list:
            # align the current sequence to find the best matches
            scores = []
            min_match = max(1, int(round(len(seq) * min_region_match, 0)))
            for p in range(last_min_end, len(input_sequence) - len(seq) + 1):
                score = 0
                for i in range(0, len(seq)):
                    if input_sequence[p + i].upper() == seq[i].upper():
                        score += 1
                if score > min_match:
                    scores.append((Interval(p + 1, p + len(seq)), score))
            if len(scores) == 0:
                raise UserWarning('could not align a given region', seq)
            best_score = max([s[1] for s in scores])
            best = [s for s in scores if s[1] == best_score]
            results.append(best)
            last_min_end = min([s[0].end for s in best])
        # take the best score for each region and see if they work in sequence
        combinations = [[p] for p in results[0]]
        for scores in results[1:]:
            temp_combos = []
            for pos, score in scores:
                for curr in combinations:
                    if pos.start > curr[-1][0].end:
                        temp_combos.append(curr[:] + [(pos, score)])
            combinations = temp_combos
            if len(combinations) == 0:
                break

        # compute the cumulative scores for the putative alignments
        best_score = None
        best_scoring_alignments = []
        for alignment in combinations:
            cumu_score = sum([s[1] for s in alignment])
            if len(best_scoring_alignments) == 0 or cumu_score > best_score:
                best_scoring_alignments = [alignment]
                best_score = cumu_score
            elif cumu_score == best_score:
                best_scoring_alignments.append(alignment)

        if len(best_scoring_alignments) == 0:
            raise UserWarning('could not map the sequences to the input')
        elif len(best_scoring_alignments) > 1:
            raise UserWarning('multiple mappings of equal score')
        else:
            alignment = [s[0] for s in best_scoring_alignments[0]]
            regions = []
            for itvl, seq in zip(alignment, seq_list):
                regions.append(DomainRegion(itvl.start, itvl.end, seq))
            return best_score, total, regions


class Translation(BioInterval):
    def __init__(self, start, end, transcript=None, domains=None, seq=None, name=None):
        """
        describes the splicing pattern and cds start and end with reference to a particular transcript

        Args:
            start (int): start of the coding sequence (cds) relative to the start of the first exon in the transcript
            end (int): end of the coding sequence (cds) relative to the start of the first exon in the transcript
            transcript (Transcript): the transcript this is a Translation of
            domains (:class:`list` of :any:`Domain`): a list of the domains on this translation
            sequence (str): the cds sequence
        """
        domains = [] if domains is None else domains
        BioInterval.__init__(self, reference_object=transcript, name=name, start=start, end=end, seq=seq)
        self.domains = [d for d in domains]

        if start <= 0:
            raise AttributeError('start must be a positive integer')
        if transcript and end > len(transcript):
            raise AttributeError('translation cannot be outside of related transcript range', end, len(transcript))

        for d in domains:
            d.reference_object = self

    @property
    def transcript(self):
        """:class:`~mavis.annotate.genomic.Transcript`: the spliced transcript this translation belongs to
        """
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
            int: the protein/amino-acid position

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
        """
        converts a genomic position to its cds (coding sequence) equivalent

        Args:
            pos (int): the genomic position

        Returns:
            int: the cds position (negative if before the initiation start site)
        """
        cds, shift = self.convert_genomic_to_nearest_cds(pos)
        if shift != 0:
            raise IndexError('conversion failed. position is outside the exonic region')
        return cds

    def convert_genomic_to_nearest_cds(self, pos):
        """
        converts a genomic position to its cds equivalent or (if intronic) the nearest cds and shift

        Args:
            pos (int): the genomic position

        Returns:
            tuple of int and int:
                * *int* - the cds position
                * *int* - the intronic shift

        """
        c, shift = self.transcript.convert_genomic_to_nearest_cdna(pos)
        if c >= self.start:
            c = c - self.start + 1
        else:
            c -= self.start
        return c, shift

    def convert_genomic_to_cds_notation(self, pos):
        """
        converts a genomic position to its cds (coding sequence) equivalent using
        `hgvs <http://www.hgvs.org/mutnomen/recs-DNA.html>`_ cds notation

        Args:
            pos (int): the genomic position

        Returns:
            str: the cds position notation

        Example:
            >>> tl = Translation(...)
            # a position before the translation start
            >>> tl.convert_genomic_to_cds_notation(1010)
            '-50'
            # a position after the translation end
            >>> tl.convert_genomic_to_cds_notation(2031)
            '*72'
            # an intronic position
            >>> tl.convert_genomic_to_cds_notation(1542)
            '50+10'
            >>> tl.convert_genomic_to_cds_notation(1589)
            '51-14'
        """
        c, shift = self.convert_genomic_to_nearest_cds(pos)
        sc = ''
        if shift > 0:
            sc = '+{}'.format(shift)
        elif shift < 0:
            sc = str(shift)

        if c > len(self):
            return '*{}{}'.format(c - len(self), sc)
        return '{}{}'.format(c, sc)

    def get_cds_seq(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            str: the cds sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        if self.seq and not ignore_cache:
            return self.seq
        elif self.transcript and self.transcript.get_strand():
            seq = self.transcript.get_seq(REFERENCE_GENOME, ignore_cache)
            return seq[self.start - 1:self.end]
        raise NotSpecifiedError('insufficient seq information')

    def get_seq(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        wrapper for the sequence method

        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name
        """
        return self.get_cds_seq(REFERENCE_GENOME, ignore_cache)

    def get_AA_seq(self, REFERENCE_GENOME=None, ignore_cache=False):
        """
        Args:
            REFERENCE_GENOME (:class:`dict` of :class:`Bio.SeqRecord` by :class:`str`): dict of reference sequence
                by template/chr name

        Returns:
            str: the amino acid sequence

        Raises:
            AttributeError: if the reference sequence has not been given and is not set
        """
        cds = self.get_cds_seq(REFERENCE_GENOME, ignore_cache)
        return translate(cds)

    def key(self):
        """see :func:`structural_variant.annotate.base.BioInterval.key`"""
        return BioInterval.key(self), self.splicing_pattern
