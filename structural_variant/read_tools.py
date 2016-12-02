import pysam
import warnings
from structural_variant.constants import *


class BamCache:
    """
    caches reads by name to facilitate getting read mates without jumping around
    the file if we've already read that section
    """
    def __init__(self, bamfile):
        """
        Args:
            bamfile (str): path to the input bam file
        """
        self.cache = {}
        self.fh = bamfile
        if not hasattr(bamfile, 'fetch'):
            self.fh = pysam.AlignmentFile(bamfile, 'rb')

    def add_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read to add to the cache
        """
        self.cache.setdefault(read.query_name, set()).add(read)

    def reference_id(self, chrom):
        """
        Args:
            chrom (str): the chromosome/reference name
        Returns:
            int: the reference id corresponding to input chromosome name
        """
        tid = self.fh.get_tid(chrom)
        if tid == -1:
            raise KeyError('invalid reference name not present in bam file')
        return tid

    def chr(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read we want the chromosome name for
        Returns:
            str: the name of the chromosome
        """
        return self.fh.get_reference_name(read.reference_id)

    @classmethod
    def _generate_fetch_bins(cls, start, stop, sample_bins, bin_gap_size):
        """
        Args:
            start (int): the start if the area to fetch reads from
            stop (int): the end of the region
            sample_bins (int): the number of bins to split the region into
            bin_gap_size (int): the space to skip between bins
        """
        bin_size = int(((stop - start + 1) - bin_gap_size * (sample_bins - 1)) / sample_bins)

        fetch_regions = [(start, start + bin_size - 1)]  # exclusive ranges for fetch
        for i in range(0, sample_bins - 1):
            st = fetch_regions[-1][1] + bin_gap_size + 1
            end = st + bin_size
            fetch_regions.append((st, end))
        fetch_regions[-1] = fetch_regions[-1][0], stop
        return fetch_regions

    def fetch(self, chrom, start, stop, read_limit=10000, cache=False, sample_bins=3,
              cache_if=lambda x: True, bin_gap_size=0):
        """
        wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
        position from within the loop. Also caches reads if requested and can return a limited read number
        """
        # try using the cache to avoid fetching regions more than once
        result = []
        bin_limit = int(read_limit / sample_bins) if read_limit else None
        # split into multiple fetches based on the 'sample_bins'
        for fstart, fend in self.__class__._generate_fetch_bins(start, stop, sample_bins, bin_gap_size):
            count = 0
            for read in self.fh.fetch(chrom, fstart, fend):
                if bin_limit is not None and count >= bin_limit:
                    warnings.warn(
                        'hit read limit. Fetch will not cache all reads for this bin: {}-{}'.format(fstart, fend))
                    break
                result.append(read)
                if cache:
                    if cache_if(read):
                        self.add_read(read)
                count += 1
        return set(result)

    def get_mate(self, read, primary_only=True, allow_file_access=True):
        """
        Args:
            read (pysam.AlignedSegment): the read
            primary_only (boolean, default=True): ignore secondary alignments
            allow_file_access (boolean, default=True): determines if the bam can be accessed to try to find the mate
        Returns:
            List[pysam.AlignedSegment]: list of mates of the input read
        """
        # NOTE: will return all mate alignments that have been cached
        putative_mates = self.cache.get(read.query_name, set())
        # if SUFFIX_DELIM in read.query_name:
        #     prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
        #     putative_mates.update(self.cache.get(temp, set()))
        mates = []
        for mate in putative_mates:
            if any([read.is_read1 == mate.is_read1,
                    read.is_read2 == mate.is_read2,
                    read.next_reference_start != mate.reference_start,
                    read.next_reference_id != mate.reference_id,
                    primary_only and mate.is_secondary]):
                continue
            mates.append(mate)
        if len(mates) == 0:
            if not allow_file_access or read.mate_is_unmapped:
                raise KeyError('mate is not found in the cache')
            else:
                warnings.warn(
                    'looking for uncached mate of {0}. This requires file access and'.format(read.query_name) +
                    ' requests may be slow. This should also not be using in a loop iterating using the file pointer ' +
                    ' as it will change the file pointer position')
                m = self.fh.mate(read)
                self.add_read(m)
                return [m]
        return mates

    def close(self):
        """
        close the bam file handle
        """
        self.fh.close()


class CigarTools:
    """
    holds methods related to processing cigar tuples. Cigar tuples are generally
    an iterable list of tuples where the first element in each tuple is the
    CIGAR value (i.e. 1 for an insertion), and the second value is the frequency
    """
    @classmethod
    def recompute_cigar_mismatch(cls, read, ref):
        """
        for cigar tuples where M is used, recompute to replace with X/= for increased
        utility and specificity
        """
        temp = []
        offset = 0

        ref_pos = read.reference_start
        seq_pos = 0

        for cigar_value, freq in read.cigar:
            if cigar_value in [CIGAR.S, CIGAR.I]:
                temp.append((cigar_value, freq))
                seq_pos += freq
            elif cigar_value in [CIGAR.D, CIGAR.N]:
                temp.append((cigar_value, freq))
                ref_pos += freq
            elif cigar_value in [CIGAR.M, CIGAR.X, CIGAR.EQ]:
                for offset in range(0, freq):
                    if DNA_ALPHABET.match(ref[ref_pos], read.query_sequence[seq_pos]):
                        if len(temp) == 0 or temp[-1][0] != CIGAR.EQ:
                            temp.append((CIGAR.EQ, 1))
                        else:
                            temp[-1] = (CIGAR.EQ, temp[-1][1] + 1)
                    else:
                        if len(temp) == 0 or temp[-1][0] != CIGAR.X:
                            temp.append((CIGAR.X, 1))
                        else:
                            temp[-1] = (CIGAR.X, temp[-1][1] + 1)
                    ref_pos += 1
                    seq_pos += 1
            else:
                raise NotImplementedError('unexpected CIGAR value {0} is not supported currently'.format(cigar_value))
        assert(sum([x[1] for x in temp]) == sum(x[1] for x in read.cigar))
        return temp

    @classmethod
    def longest_fuzzy_match(cls, cigar, max_fuzzy_interupt=1):
        """
        computes the longest sequence of exact matches allowing for 'x' event interrupts

        """
        temp = CigarTools.join(cigar)
        longest_fuzzy_match = 0
        for pos, c in enumerate(temp):
            if c[0] != CIGAR.EQ:
                continue
            current_fuzzy_match = c[1]
            pos += 1
            fuzzy_count = 0
            while pos < len(temp) and fuzzy_count <= max_fuzzy_interupt:
                v, f = temp[pos]
                if v != CIGAR.EQ:
                    fuzzy_count += 1
                else:
                    current_fuzzy_match += f
                pos += 1
            if current_fuzzy_match > longest_fuzzy_match:
                longest_fuzzy_match = current_fuzzy_match
        return longest_fuzzy_match

    @classmethod
    def longest_exact_match(cls, cigar):
        """
        returns the longest consecutive exact match
        """
        return CigarTools.longest_fuzzy_match(cigar, 0)

    @classmethod
    def score(cls, cigar, **kwargs):
        """scoring based on sw alignment properties with gap extension penalties

        Args:
            cigar (List<(CIGAR,int)>): list of cigar tuple values
            MISMATCH (int, default=-1): mismatch penalty
            MATCH (int, default=2): match penalty
            GAP (int, default=-4): initial gap penalty
            GAP_EXTEND (int, default=-1): gap extension penalty

        Returns:
            int: the score value
        """

        MISMATCH = kwargs.pop('MISMATCH', -1)
        MATCH = kwargs.pop('MATCH', 2)
        GAP = kwargs.pop('GAP', -4)
        GAP_EXTEND = kwargs.pop('GAP_EXTEND', -1)

        score = 0
        for v, freq in cigar:
            if v == CIGAR.EQ:
                score += MATCH * freq
            elif v == CIGAR.X:
                score += MISMATCH * freq
            elif v in [CIGAR.I, CIGAR.D]:
                score += GAP + GAP_EXTEND * (freq - 1)
            elif v in [CIGAR.S, CIGAR.N]:
                pass
            else:
                raise AssertionError('unexpected cigar value', v)
        return score

    @classmethod
    def match_percent(cls, cigar):
        """
        calculates the percent of aligned bases (matches or mismatches) that are matches
        """
        matches = 0
        mismatches = 0
        for v, f in cigar:
            if v == CIGAR.EQ:
                matches += f
            elif v == CIGAR.X:
                mismatches += f
            elif v == CIGAR.M:
                raise AttributeError('cannot calculate match percent with non-specific alignments', cigar)
        if matches + mismatches == 0:
            raise AttributeError('input cigar string does not have any aligned sections (X or =)', cigar)
        else:
            return matches / (matches + mismatches)

    @classmethod
    def join(cls, *pos):
        """
        given a number of cigar lists, joins them and merges any consecutive tuples
        with the same cigar value

        Example:
            >>> CigarTools.join([(1, 1), (4, 7)], [(4, 3), (2, 4)])
            [(1, 1), (4, 10), (2, 4)]
        """
        result = []
        for cigar in pos:
            for v, f in cigar:
                if len(result) > 0 and result[-1][0] == v:
                    result[-1] = (v, f + result[-1][1])
                else:
                    result.append((v, f))
        return result

    @classmethod
    def extend_softclipping(cls, cigar, min_exact_to_stop_softclipping):
        """
        given some input cigar, extends softclipping if there are mismatches/insertions/deletions
        close to the end of the aligned portion. The stopping point is defined by the
        min_exact_to_stop_softclipping parameter. this function will throw an error if there is no
        exact match aligned portion to signal stop

        Args:
            original_cigar (List[CIGAR,int]): the input cigar
            min_exact_to_stop_softclipping (int): number of exact matches to terminate extension

        Returns:
            (List[CIGAR,int], int): the new cigar string and a number representing the shift from the original
                start position
        """
        ref_start_shift = 0
        # determine how far to scoop for the front softclipping
        new_cigar = []
        temp = [(v, f) for v, f in cigar if (v in [CIGAR.EQ, CIGAR.M] and f >= min_exact_to_stop_softclipping)]
        if len(temp) == 0:
            raise AttributeError('cannot compute on this cigar as there is no stop point')

        match_satisfied = False
        for v, f in cigar:
            if v in [CIGAR.M, CIGAR.EQ] and f >= min_exact_to_stop_softclipping:
                match_satisfied = True
            if match_satisfied and (len(new_cigar) == 0 or new_cigar[-1][0] != CIGAR.S):
                new_cigar.append((v, f))
            elif match_satisfied:  # first after SC
                if v in [CIGAR.D, CIGAR.X, CIGAR.N]:
                    ref_start_shift += f
                    new_cigar.append((CIGAR.S, f))
                elif v == CIGAR.I:
                    new_cigar.append((CIGAR.S, f))
                else:
                    new_cigar.append((v, f))
            else:
                if v in [CIGAR.D, CIGAR.N]:
                    ref_start_shift += f
                    pass
                elif v in [CIGAR.I, CIGAR.S]:
                    new_cigar.append((CIGAR.S, f))
                else:
                    new_cigar.append((CIGAR.S, f))
                    ref_start_shift += f
        cigar = new_cigar[::-1]

        new_cigar = []

        match_satisfied = False
        for v, f in cigar:
            if v in [CIGAR.M, CIGAR.EQ] and f >= min_exact_to_stop_softclipping:
                match_satisfied = True
            if match_satisfied and (len(new_cigar) == 0 or new_cigar[-1][0] != CIGAR.S):
                new_cigar.append((v, f))
            elif match_satisfied:  # first after SC
                if v in [CIGAR.D, CIGAR.X, CIGAR.N]:
                    new_cigar.append((CIGAR.S, f))
                elif v == CIGAR.I:
                    new_cigar.append((CIGAR.S, f))
                else:
                    new_cigar.append((v, f))
            else:
                if v in [CIGAR.D, CIGAR.N]:
                    pass
                elif v in [CIGAR.I, CIGAR.S]:
                    new_cigar.append((CIGAR.S, f))
                else:
                    new_cigar.append((CIGAR.S, f))
        new_cigar.reverse()
        return new_cigar, ref_start_shift

    @classmethod
    def compute(cls, ref, alt, force_softclipping=True, min_exact_to_stop_softclipping=6):
        """
        given a ref and alt sequence compute the cigar string representing the alt

        returns the cigar tuples along with the start position of the alt relative to the ref
        """
        if not force_softclipping:
            min_exact_to_stop_softclipping = 1

        if len(ref) != len(alt):
            raise AttributeError('ref and alt must be the same length')
        cigar = []
        for r, a in zip(ref, alt):
            if r == a and r == GAP:
                pass
            elif r == GAP:
                cigar.append((CIGAR.I, 1))
            elif a == GAP:
                cigar.append((CIGAR.D, 1))
            elif DNA_ALPHABET.match(r, a):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))
        cigar = cls.join(cigar)

        try:
            c, rs = CigarTools.extend_softclipping(cigar, min_exact_to_stop_softclipping)
            return c, rs
        except AttributeError:
            return cigar, 0

    @classmethod
    def convert_for_igv(cls, cigar):
        """
        igv does not support the extended CIGAR values for match v mismatch

        Example:
            >>> CigarTools.convert_for_igv([(7, 4), (8, 1), (7, 5)])
            [(0, 10)]
        """
        result = []
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ]:
                v = CIGAR.M
            result.append((v, f))
        return cls.join(result)

    @classmethod
    def alignment_matches(cls, cigar):
        """
        counts the number of aligned bases irrespective of match/mismatch
        this is equivalent to counting all CIGAR.M
        """
        result = 0
        for v, f in cigar:
            if v in [CIGAR.X, CIGAR.EQ, CIGAR.M]:
                result += f
        return result


def breakpoint_pos(read, orient=ORIENT.NS):
    """
    assumes the breakpoint is the position following softclipping on the side with more
    softclipping (unless and orientation has been specified)

    Args:
        read (psyam.AlignedSegment): the read object
        orient (ORIENT): the orientation

    Returns:
        int: the position of the breakpoint in the input read
    """
    typ, freq = read.cigar[0]
    end_typ, end_freq = read.cigar[-1]
    ORIENT.enforce(orient)

    if typ != CIGAR.S and end_typ != CIGAR.S:
        raise AttributeError('cannot compute breakpoint for a read without soft-clipping')

    if orient == ORIENT.NS:
        if (typ == CIGAR.S and end_typ == CIGAR.S and freq > end_freq) \
                or typ == CIGAR.S and end_typ != CIGAR.S:
            orient = ORIENT.RIGHT
            # soft clipped to the left
        else:
            # soft clipped to the right
            orient = ORIENT.LEFT

    if orient == ORIENT.RIGHT:
        if typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint')
        return read.reference_start
    else:
        if end_typ != CIGAR.S:
            raise AttributeError('soft clipping doesn\'t support input orientation for a breakpoint')
        return read.reference_end - 1


def nsb_align(ref, seq, weight_of_score=0.5, min_overlap_percent=100):
    """
    given some reference string and a smaller sequence string computes the best non-space-breaking alignment
    i.e. an alignment that does not allow for indels (straight-match)
    """
    if len(ref) < 1 or len(seq) < 1:
        raise AttributeError('cannot overlap on an empty sequence')

    if min_overlap_percent <= 0 or min_overlap_percent > 100:
        raise AttributeError('percent must be greater than 0 and up to 100', min_overlap_percent)

    min_overlap = int(round(min_overlap_percent * len(seq) / 100, 0))

    # store to improve speed and space (don't need to store all alignments)
    best_score = 0
    results = []

    for ref_start in range(min_overlap - len(seq), len(ref) + len(seq) - min_overlap):
        score = 0
        cigar = []
        for i in range(0, len(seq)):
            r = ref_start + i
            if r < 0 or r >= len(ref):
                cigar.append((CIGAR.S, 1))
                continue
            if DNA_ALPHABET.match(ref[r], seq[i]):
                cigar.append((CIGAR.EQ, 1))
            else:
                cigar.append((CIGAR.X, 1))

        cigar = CigarTools.join(cigar)
        # end mismatches we set as soft-clipped
        if cigar[0][0] == CIGAR.X:
            cigar[0] = (CIGAR.S, cigar[0][1])
        if cigar[-1][0] == CIGAR.X:
            cigar[-1] = (CIGAR.S, cigar[-1][1])

        qstart = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]

        score = CigarTools.score(cigar) * weight_of_score + \
            CigarTools.longest_exact_match(cigar) * (1 - weight_of_score)
        score = CigarTools.score(cigar)
        a = pysam.AlignedSegment()
        a.query_sequence = str(seq)
        a.reference_start = ref_start + qstart
        a.cigar = cigar

        if score >= best_score:
            best_score = score
            results.append((a, score))

    filtered = [x for x, y in results if y == best_score]

    return filtered


def median_insert_size(reads):
    isizes = sorted([abs(read.template_length) for read in reads])
    if len(isizes) % 2 == 0:
        i = len(isizes) // 2
        return (isizes[i - 1] + isizes[i]) / 2
    else:
        return isizes[len(isizes) // 2]


def stderr_insert_size(reads, centre):
    err = 0
    for read in reads:
        err += math.pow(abs(read.template_length) - centre, 2)
    return err


def stdev_insert_size(reads, centre):
    return math.sqrt(stderr_insert_size(reads, centre))
