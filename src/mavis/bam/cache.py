import atexit
import re
from typing import Callable, Dict, List, Set, Union

import pysam

from .. import util as _util
from ..annotate.base import ReferenceName
from ..interval import Interval
from .read import SamRead


class BamCache:
    """
    caches reads by name to facilitate getting read mates without jumping around
    the file if we've already read that section
    """

    fh: pysam.AlignmentFile
    stranded: bool
    cache: Dict

    def __init__(self, bamfile: Union[pysam.AlignmentFile, str], stranded: bool = False):
        """
        Args:
            bamfile: path to the input bam file
        """
        self.cache = {}
        self.stranded = stranded

        if not hasattr(bamfile, 'fetch'):
            self.fh = pysam.AlignmentFile(bamfile, 'rb')
        else:
            try:
                self.fh = bamfile.fh
            except AttributeError:
                self.fh = bamfile
        atexit.register(self.close)  # makes the file 'auto close' on normal python exit

    def valid_chr(self, chrom: str) -> bool:
        """
        checks if a reference name exists in the bam file header
        """
        try:
            self.reference_id(chrom)
            return True
        except KeyError:
            return False

    def add_read(self, read: pysam.AlignedSegment):
        """
        Args:
            read: the read to add to the cache
        """
        if not read.is_unmapped and read.reference_start == read.reference_end:
            _util.logger.debug(f'ignoring invalid read: {read.query_name}')
            return
        if not isinstance(read, SamRead):
            read = SamRead.copy(read)
        self.cache.setdefault(read.query_name, set())
        if read not in self.cache[read.query_name]:
            self.cache[read.query_name].add(read)

    def has_read(self, read: pysam.AlignedSegment) -> bool:
        """
        checks if a read query name exists in the current cache
        """
        if read.query_name not in self.cache:
            return False
        elif read in self.cache[read.query_name]:
            return True
        return False

    def reference_id(self, chrom: str) -> int:
        """
        Args:
            chrom: the chromosome/reference name
        Returns:
            the reference id corresponding to input chromosome name
        """
        tid = self.fh.get_tid(chrom)
        if tid == -1:
            tid = self.fh.get_tid(re.sub('^chr', '', chrom))
        if tid == -1:
            tid = self.fh.get_tid('chr' + chrom)
        if tid == -1:
            raise KeyError('invalid reference name not present in bam file', chrom)
        return tid

    def get_read_reference_name(self, read: pysam.AlignedSegment) -> str:
        """
        Args:
            read: the read we want the chromosome name for
        Returns:
            the name of the chromosome
        """
        return ReferenceName(self.fh.get_reference_name(read.reference_id))

    @classmethod
    def _generate_fetch_bins(
        cls, start: int, stop: int, sample_bins: int, min_bin_size: int
    ) -> List[Interval]:
        """
        Args:
            start: the start if the area to fetch reads from
            stop: the end of the region
            sample_bins: the number of bins to split the region into
            min_bin_size: the minimum bin size
        """
        assert min_bin_size > 0
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
        self,
        input_chrom: str,
        start: int,
        stop: int,
        limit: int = 10000,
        cache_if: Callable = lambda x: True,
        filter_if: Callable = lambda x: False,
        stop_on_cached_read: bool = False,
    ) -> Set[pysam.AlignedSegment]:
        """
        Args:
            input_chrom: chromosome name
            start: start position
            end: end position
            limit: maximum number of reads to fetch
            cache_if:  if returns True then the read is added to the cache
            filter_if: if returns True then the read is not returned as part of the result
            stop_on_cached_read: stop reading at the first read found that is already in the cache
        Note:
            the cache_if and filter_if functions must be any function that takes a read as input and returns a boolean

        Returns:
            a set of reads which overlap the input region
        """
        # try using the cache to avoid fetching regions more than once
        result = []
        chrom = input_chrom
        # split into multiple fetches based on the 'sample_bins'
        if str(chrom) not in self.fh.references:
            chrom = re.sub('^chr', '', chrom)
            if chrom not in self.fh.references:
                chrom = 'chr' + chrom
            if chrom not in self.fh.references:
                raise KeyError(
                    'bam file does not contain the expected reference',
                    input_chrom,
                    self.fh.references,
                )
        temp_cache = set()
        count = 0

        for read in self.fh.fetch(chrom, start, stop):
            if limit is not None and count >= limit:
                break
            if stop_on_cached_read and self.has_read(read):
                break
            if not read.is_unmapped and read.reference_start == read.reference_end:
                _util.logger.debug(f'ignoring invalid read {read.query_name}')
                continue
            read = SamRead.copy(read)
            if not filter_if(read):
                result.append(read)
            if cache_if(read):
                self.add_read(read)
            if read.query_name not in temp_cache:
                count += 1
                temp_cache.add(read.query_name)
        return set(result)

    def fetch_from_bins(
        self,
        input_chrom: str,
        start: int,
        stop: int,
        read_limit: int = 10000,
        cache: bool = False,
        sample_bins: int = 3,
        cache_if: Callable = lambda x: True,
        min_bin_size: int = 10,
        filter_if: Callable = lambda x: False,
    ) -> Set[pysam.AlignedSegment]:
        """
        wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
        position from within the loop. Also caches reads if requested and can return a limited read number

        Args:
            input_chrom: the chromosome
            start: the start position
            stop: the end position
            read_limit: the maximum number of reads to parse
            cache: flag to store reads
            sample_bins: number of bins to split the region into
            cache_if: function to check to against a read to determine if it should be cached
            bin_gap_size: gap between the bins for the fetch area

        Returns:
            set of reads gathered from the region
        """
        # try using the cache to make grabbing mate pairs easier
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
        bins: List[Interval] = self.__class__._generate_fetch_bins(
            start, stop, sample_bins, min_bin_size
        )
        running_surplus = 0
        temp_cache = set()
        for fstart, fend in bins:
            count = 0
            running_surplus += bin_limit

            for read in self.fh.fetch(chrom, fstart, fend):
                if bin_limit is not None and count >= running_surplus:
                    break
                if not read.is_unmapped and read.reference_start == read.reference_end:
                    _util.logger.debug(f'ignoring invalid read {read.query_name}')
                    continue
                read = SamRead.copy(read)
                if not filter_if(read):
                    result.append(read)
                if read.query_name not in temp_cache:
                    count += 1
                    temp_cache.add(read.query_name)
                if cache and cache_if(read):
                    self.add_read(read)
            running_surplus -= count
        return set(result)

    def get_mate(
        self, read: pysam.AlignedSegment, primary_only: bool = True, allow_file_access: bool = False
    ) -> List[pysam.AlignedSegment]:
        """
        Args:
            read: the read
            primary_only: ignore secondary alignments
            allow_file_access: determines if the bam can be accessed to try to find the mate
        Returns:
            list of mates of the input read
        """
        # NOTE: will return all mate alignments that have been cached
        putative_mates = self.cache.get(read.query_name, set())
        mates = []
        for mate in putative_mates:
            if not mate.is_unmapped:
                if any(
                    [
                        read.is_read1 == mate.is_read1,
                        read.next_reference_start != mate.reference_start,
                        read.next_reference_id != mate.reference_id,
                        primary_only and (mate.is_secondary or mate.is_supplementary),
                        abs(read.template_length) != abs(mate.template_length),
                    ]
                ):
                    continue
            mates.append(mate)
        if len(mates) == 0:
            if not allow_file_access or read.mate_is_unmapped:
                raise KeyError('mate is not found in the cache')
            else:
                _util.logger.warning(
                    f'looking for uncached mate of {read.query_name}. This requires file access and'
                    ' requests may be slow. This should also not be using in a loop iterating using the file pointer '
                    ' as it will change the file pointer position'
                )
                m = self.fh.mate(read)
                m = SamRead.copy(m)
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
