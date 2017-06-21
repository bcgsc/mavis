import pysam
import warnings
import atexit
import re
from ..annotate.base import ReferenceName
from ..interval import Interval


class BamCache:
    """
    caches reads by name to facilitate getting read mates without jumping around
    the file if we've already read that section
    """

    def __init__(self, bamfile, stranded=False):
        """
        Args:
            bamfile (str): path to the input bam file
        """
        self.cache = {}
        self.stranded = stranded
        self.fh = bamfile
        if not hasattr(bamfile, 'fetch'):
            self.fh = pysam.AlignmentFile(bamfile, 'rb')
        else:
            try:
                self.fh = bamfile.fh
            except AttributeError:
                pass
        atexit.register(self.close)  # makes the file 'auto close' on normal python exit

    def add_read(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read to add to the cache
        """
        self.cache.setdefault(read.query_name, set()).add(read)

    def has_read(self, read):
        if read.query_name in self.cache:
            return True
        else:
            return False

    def reference_id(self, chrom):
        """
        Args:
            chrom (str): the chromosome/reference name
        Returns:
            int: the reference id corresponding to input chromosome name
        """
        tid = self.fh.get_tid(chrom)
        if tid == -1:
            tid = self.fh.get_tid(re.sub('^chr', '', chrom))
        if tid == -1:
            tid = self.fh.get_tid('chr' + chrom)
        if tid == -1:
            raise KeyError('invalid reference name not present in bam file', chrom)
        return tid

    def get_read_reference_name(self, read):
        """
        Args:
            read (pysam.AlignedSegment): the read we want the chromosome name for
        Returns:
            str: the name of the chromosome
        """
        return ReferenceName(self.fh.get_reference_name(read.reference_id))

    @classmethod
    def _generate_fetch_bins(cls, start, stop, sample_bins, min_bin_size):
        """
        Args:
            start (int): the start if the area to fetch reads from
            stop (int): the end of the region
            sample_bins (int): the number of bins to split the region into
            min_bin_size (int): the minimum bin size
        """
        assert(min_bin_size > 0)
        length = stop - start + 1
        bin_sizes = []
        if sample_bins * min_bin_size > length:
            sample_bins = max([1, length // min_bin_size])
        s = length // sample_bins  # base size
        bin_sizes = [s for i in range(0, sample_bins)]
        remainder = (length - sum(bin_sizes)) // sample_bins
        bin_sizes = [b + remainder for b in bin_sizes]
        unsplittable_remainder = length - sum(bin_sizes)
        for i in range(0, unsplittable_remainder):
            bin_sizes[i] += 1
        fetch_regions = [(start, start + bin_sizes[0] - 1)]
        for b in bin_sizes[1:]:
            last_end = fetch_regions[-1][1]
            fetch_regions.append((last_end + 1, last_end + b))
        return [Interval(s, t) for s, t in fetch_regions]

    def fetch(
        self, input_chrom, start, stop, limit=10000, cache_if=lambda x: True, filter_if=lambda x: False
    ):
        # try using the cache to avoid fetching regions more than once
        result = []
        chrom = input_chrom
        # split into multiple fetches based on the 'sample_bins'
        if str(chrom) not in self.fh.references:
            chrom = re.sub('^chr', '', chrom)
            if chrom not in self.fh.references:
                chrom = 'chr' + chrom
            if chrom not in self.fh.references:
                raise KeyError('bam file does not contain the expected reference', input_chrom)
        temp_cache = set()
        count = 0

        for read in self.fh.fetch(chrom, start, stop):
            if limit is not None and count >= limit:
                break
            if not filter_if(read):
                result.append(read)
            if cache_if(read):
                self.add_read(read)
            if read.query_name not in temp_cache:
                count += 1
                temp_cache.add(read.query_name)
        return set(result)

    def fetch_from_bins(
        self, input_chrom, start, stop, read_limit=10000, cache=False, sample_bins=3,
        cache_if=lambda x: True, min_bin_size=10, filter_if=lambda x: False
    ):
        """
        wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
        position from within the loop. Also caches reads if requested and can return a limited read number

        Args:
            chrom (str): the chromosome
            start (int): the start position
            stop (int): the end position
            read_limit (int): the maximum number of reads to parse
            cache (bool): flag to store reads
            sample_bins (int): number of bins to split the region into
            cache_if (callable): function to check to against a read to determine if it should be cached
            bin_gap_size (int): gap between the bins for the fetch area

        Returns:
            :class:`set` of :class:`pysam.AlignedSegment`: set of reads gathered from the region
        """
        # try using the cache to avoid fetching regions more than once
        result = []
        bin_limit = int(read_limit / sample_bins) if read_limit else None
        chrom = input_chrom
        # split into multiple fetches based on the 'sample_bins'
        if str(chrom) not in self.fh.references:
            chrom = re.sub('^chr', '', chrom)
            if chrom not in self.fh.references:
                chrom = 'chr' + chrom
            if chrom not in self.fh.references:
                raise KeyError('bam file does not contain the expected reference', input_chrom)
        bins = self.__class__._generate_fetch_bins(start, stop, sample_bins, min_bin_size)
        running_surplus = 0
        temp_cache = set()
        for fstart, fend in bins:
            count = 0
            running_surplus += bin_limit

            for read in self.fh.fetch(chrom, fstart, fend):
                if bin_limit is not None and count >= running_surplus:
                    break
                if not filter_if(read):
                    result.append(read)
                if read.query_name not in temp_cache:
                    count += 1
                    temp_cache.add(read.query_name)
                if cache and cache_if(read):
                    self.add_read(read)
            running_surplus -= count
        return set(result)

    def get_mate(self, read, primary_only=True, allow_file_access=False):
        """
        Args:
            read (pysam.AlignedSegment): the read
            primary_only (bool): ignore secondary alignments
            allow_file_access (bool): determines if the bam can be accessed to try to find the mate
        Returns:
            :class:`list` of :class:`pysam.AlignedSegment`: list of mates of the input read
        """
        # NOTE: will return all mate alignments that have been cached
        putative_mates = self.cache.get(read.query_name, set())
        mates = []
        for mate in putative_mates:
            if not mate.is_unmapped:
                if any([
                    read.is_read1 == mate.is_read1,
                    read.next_reference_start != mate.reference_start,
                    read.next_reference_id != mate.reference_id,
                    primary_only and (mate.is_secondary or mate.is_supplementary),
                    abs(read.template_length) != abs(mate.template_length)
                ]):
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
        try:
            self.fh.close()
        except AttributeError:
            pass
