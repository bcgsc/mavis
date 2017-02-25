import pysam
import warnings


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

    def fetch(
        self, chrom, start, stop, read_limit=10000, cache=False, sample_bins=3,
        cache_if=lambda x: True, bin_gap_size=0, filter_if=lambda x: False
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
        # split into multiple fetches based on the 'sample_bins'
        for fstart, fend in self.__class__._generate_fetch_bins(start, stop, sample_bins, bin_gap_size):
            count = 0
            for read in self.fh.fetch(chrom, fstart, fend):
                if bin_limit is not None and count >= bin_limit:
                    warnings.warn(
                        'hit read limit. Fetch will not cache all reads for this bin: {}-{}'.format(fstart, fend))
                    break
                if not filter_if(read):
                    result.append(read)
                if cache and cache_if(read):
                    self.add_read(read)
                count += 1
        return set(result)

    def get_mate(self, read, primary_only=True, allow_file_access=True):
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
        # if SUFFIX_DELIM in read.query_name:
        #     prefix, temp = read.query_name.split(SUFFIX_DELIM, 1)
        #     putative_mates.update(self.cache.get(temp, set()))
        mates = []
        for mate in putative_mates:
            if not mate.is_unmapped:
                if any([
                    read.is_read1 == mate.is_read1,
                    read.next_reference_start != mate.reference_start,
                    read.next_reference_id != mate.reference_id,
                    primary_only and mate.is_secondary,
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
        self.fh.close()
