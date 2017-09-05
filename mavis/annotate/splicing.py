from ..constants import STRAND, reverse_complement
from .constants import SPLICE_SITE_RADIUS, SPLICE_TYPE, SPLICE_SITE_TYPE, DONOR_SEQ, ACCEPTOR_SEQ
from ..interval import Interval
from .base import BioInterval


class SplicingPattern(list):

    def __init__(self, *args, splice_type=SPLICE_TYPE.NORMAL):
        list.__init__(self, *args)
        self.splice_type = splice_type

    @staticmethod
    def classify(pattern, original_sites):
        # now need to decide the type for each set
        pattern = sorted(pattern)
        r_introns = 0
        s_exons = 0
        assert(len(pattern) % 2 == 0)

        for d, a in zip(pattern[0::2], pattern[1::2]):
            # check if any original splice positions are between this donor and acceptor
            temp = 0
            for s in original_sites:
                if s > d and s < a:
                    temp += 1
            assert(temp % 2 == 0)
            s_exons += temp // 2

        for a, d in zip(pattern[1::2], pattern[2::2]):
            temp = 0
            for s in original_sites:
                if s > a and s < d:
                    temp += 1
            assert(temp % 2 == 0)
            r_introns += temp // 2

        if len(pattern) > 0:
            # any skipped positions before the first donor or after the last acceptor
            temp = 0
            for s in original_sites:
                if s < pattern[0]:
                    temp += 1
            assert(temp % 2 == 0)
            r_introns += temp // 2
            temp = 0
            for s in original_sites:
                if s > pattern[-1]:
                    temp += 1
            r_introns += temp // 2
            assert(temp % 2 == 0)

        # now classifying the pattern
        if r_introns + s_exons == 0:
            return SPLICE_TYPE.NORMAL
        elif r_introns == 0:
            if s_exons > 1:
                return SPLICE_TYPE.MULTI_SKIP
            else:
                return SPLICE_TYPE.SKIP
        elif s_exons == 0:
            if r_introns > 1:
                return SPLICE_TYPE.MULTI_RETAIN
            else:
                return SPLICE_TYPE.RETAIN
        else:
            return SPLICE_TYPE.COMPLEX


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
        assert(pos <= self.end and pos >= self.start)
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
            d = convert_match_to_ss(match, SPLICE_SITE_TYPE.DONOR)
            if d.pos not in positions:
                sites.append(d)
                positions.add(d.pos)
    positions = set()
    for regex in ACCEPTOR_SEQ:
        for match in regex.finditer(sequence):
            d = convert_match_to_ss(match, SPLICE_SITE_TYPE.ACCEPTOR)
            if d.pos not in positions:
                sites.append(d)
                positions.add(d.pos)
    if is_reverse:
        temp = []
        l = len(sequence)
        # flip all the sites
        for site in sites:
            offset = site.end - site.pos
            start = l - site.end + 1
            new_site = SpliceSite(
                None, start=start, end=l - site.start + 1,
                seq=reverse_complement(site.seq),
                strand=STRAND.NEG,
                pos=start + offset,
                site_type=site.type)
            temp.append(new_site)
        sites = temp
    return sites
