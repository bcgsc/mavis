import itertools

from .base import BioInterval
from .constants import ACCEPTOR_SEQ, DONOR_SEQ, SPLICE_SITE_RADIUS, SPLICE_SITE_TYPE, SPLICE_TYPE
from ..constants import reverse_complement, STRAND
from ..interval import Interval


class SplicingPattern(list):

    def __init__(self, *args, splice_type=SPLICE_TYPE.NORMAL):
        list.__init__(self, *args)
        self.splice_type = splice_type

    def __str__(self):
        temp = []
        for site in self:
            temp.append('{}{}{}'.format('D' if site.type == SPLICE_SITE_TYPE.DONOR else 'A', site.pos, '' if site.intact else '*'))
        return '[{}]'.format(', '.join(temp))

    @classmethod
    def classify(cls, pattern, original_sites):
        # now need to decide the type for each set
        pattern = sorted(pattern)
        r_introns = []
        s_exons = []
        assert len(pattern) % 2 == 0

        for donor, acceptor in zip(pattern[0::2], pattern[1::2]):
            # check if any original splice positions are between this donor and acceptor
            temp = []
            for site in original_sites:
                if site > donor and site < acceptor:
                    temp.append(site)
            assert len(temp) % 2 == 0
            s_exons.extend(temp)

        for acceptor, donor in zip(pattern[1::2], pattern[2::2]):
            temp = []
            for site in original_sites:
                if site > acceptor and site < donor:
                    temp.append(site)
            assert len(temp) % 2 == 0
            r_introns.extend(temp)
        if pattern:
            # any skipped positions before the first donor or after the last acceptor
            temp = []
            for site in original_sites:
                if site < pattern[0]:
                    temp.append(site)
            assert len(temp) % 2 == 0
            r_introns.extend(temp)
            temp = []
            for site in original_sites:
                if site > pattern[-1]:
                    temp.append(site)
            r_introns.extend(temp)
            assert len(temp) % 2 == 0
        rintron_count = 0
        for i in range(0, len(r_introns) - 1):
            if abs(r_introns[i].pos - r_introns[i + 1].pos) > 1:
                rintron_count += 1
        sexon_count = len(s_exons) // 2
        # now classifying the pattern
        if rintron_count + sexon_count == 0:
            return SPLICE_TYPE.NORMAL
        elif rintron_count == 0:
            if sexon_count > 1:
                return SPLICE_TYPE.MULTI_SKIP
            return SPLICE_TYPE.SKIP
        elif sexon_count == 0:
            if rintron_count > 1:
                return SPLICE_TYPE.MULTI_RETAIN
            return SPLICE_TYPE.RETAIN
        return SPLICE_TYPE.COMPLEX

    @classmethod
    def generate_patterns(cls, sites, is_reverse=False):
        """
        returns a list of splice sites to be connected as a splicing pattern

        Returns:
            :class:`list` of :class:`SplicingPattern`: List of positions to be spliced together

        see :ref:`theory - predicting splicing patterns <theory-predicting-splicing-patterns>`
        """
        if not sites:
            return [SplicingPattern()]
        sites = sorted(sites, reverse=is_reverse)

        patterns = []
        for site in sites:
            if site.intact:
                if patterns and patterns[-1][0].type == site.type:
                    patterns[-1].append(site)
                else:
                    patterns.append([site])
        if patterns and patterns[0][0].type == SPLICE_SITE_TYPE.ACCEPTOR:
            patterns = patterns[1:]
        if patterns and patterns[-1][0].type == SPLICE_SITE_TYPE.DONOR:
            patterns = patterns[:-1]
        if not patterns:
            return [SplicingPattern()]
        patterns = list(itertools.product(*patterns))
        for i, patt in enumerate(patterns):
            patterns[i] = SplicingPattern(patt, splice_type=cls.classify(patt, sites))
        return patterns


class SpliceSite(BioInterval):

    def __init__(self, ref, pos, site_type, intact=True, start=None, end=None, strand=None, seq=None):
        if start is None or end is None:
            self.strand = strand if strand else ref.get_strand()
            if self.strand == STRAND.NEG:
                if site_type == SPLICE_SITE_TYPE.DONOR:
                    if start is None:
                        start = pos - SPLICE_SITE_RADIUS
                    if end is None:
                        end = pos + SPLICE_SITE_RADIUS - 1
                else:
                    if start is None:
                        start = pos - SPLICE_SITE_RADIUS + 1
                    if end is None:
                        end = pos + SPLICE_SITE_RADIUS
            else:
                if site_type == SPLICE_SITE_TYPE.ACCEPTOR:
                    if start is None:
                        start = pos - SPLICE_SITE_RADIUS
                    if end is None:
                        end = pos + SPLICE_SITE_RADIUS - 1
                else:
                    if start is None:
                        start = pos - SPLICE_SITE_RADIUS + 1
                    if end is None:
                        end = pos + SPLICE_SITE_RADIUS
        BioInterval.__init__(self, ref, start, end, seq=seq, strand=strand)
        assert pos <= self.end and pos >= self.start
        self.pos = pos
        self.intact = intact
        self.type = SPLICE_SITE_TYPE.enforce(site_type)

    def __or__(self, other):
        return Interval.__or__(self, other)

    def __repr__(self):
        cls = self.__class__.__name__
        refname = self.reference_object
        try:
            refname = self.reference_object.name
        except AttributeError:
            pass
        seq = '' if not self.seq else ', seq=' + self.seq
        return '{}(type={}, {}:{}({}-{}){}, strand={})'.format(
            cls, SPLICE_SITE_TYPE.reverse(self.type),
            refname, self.pos, self.start, self.end, seq, self.get_strand())


def predict_splice_sites(input_sequence, is_reverse=False):
    """
    looks for the expected splice site sequence patterns in the
    input strings and returns a list of putative splice sites

    Args:
        input_sequence (str): input sequence with respect to the positive/forward strand
        is_reverse (bool): True when the sequences is transcribed on the reverse strand

    Return:
        list of SpliceSite: list of putative splice sites
    """
    if is_reverse:
        sequence = reverse_complement(input_sequence)
    else:
        sequence = input_sequence

    def convert_match_to_ss(match, splice_type):
        prefix = match.group(1)
        suffix = match.group(2)
        return SpliceSite(
            None, start=match.start() + 1, end=match.end(),
            pos=match.start() + len(prefix),
            seq=prefix + suffix,
            site_type=splice_type,
            strand=STRAND.POS)

    sites = []
    positions = set()
    for regex in DONOR_SEQ:
        for match in regex.finditer(sequence):
            donor_site = convert_match_to_ss(match, SPLICE_SITE_TYPE.DONOR)
            if donor_site.pos not in positions:
                sites.append(donor_site)
                positions.add(donor_site.pos)
    positions = set()
    for regex in ACCEPTOR_SEQ:
        for match in regex.finditer(sequence):
            acceptor_site = convert_match_to_ss(match, SPLICE_SITE_TYPE.ACCEPTOR)
            if acceptor_site.pos not in positions:
                sites.append(acceptor_site)
                positions.add(acceptor_site.pos)
    if is_reverse:
        temp = []
        # flip all the sites
        for site in sites:
            offset = site.end - site.pos
            start = len(sequence) - site.end + 1
            new_site = SpliceSite(
                None, start=start, end=len(sequence) - site.start + 1,
                seq=reverse_complement(site.seq),
                strand=STRAND.NEG,
                pos=start + offset,
                site_type=site.type)
            temp.append(new_site)
        sites = temp
    return sites
