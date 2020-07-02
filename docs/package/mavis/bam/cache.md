# mavis.bam.cache

## class BamCache

caches reads by name to facilitate getting read mates without jumping around
the file if we've already read that section

### BamCache.\_\_init\_\_()

```python
def __init__(self, bamfile, stranded=False):
```

**Args**

- bamfile (`str`): path to the input bam file
- stranded

### BamCache.valid\_chr()

checks if a reference name exists in the bam file header

```python
def valid_chr(self, chrom):
```

**Args**

- chrom


### BamCache.has\_read()

checks if a read query name exists in the current cache

```python
def has_read(self, read):
```

**Args**

- read





### BamCache.fetch\_from\_bins()

wrapper around the fetch method, returns a list to avoid errors with changing the file pointer
position from within the loop. Also caches reads if requested and can return a limited read number

```python
def fetch_from_bins(
    self,
    input_chrom,
    start,
    stop,
    read_limit=10000,
    cache=False,
    sample_bins=3,
    cache_if=lambda x: True,
    min_bin_size=10,
    filter_if=lambda x: False,
):
```

**Args**

- input_chrom
- start (`int`): the start position
- stop (`int`): the end position
- read_limit (`int`): the maximum number of reads to parse
- cache (`bool`): flag to store reads
- sample_bins (`int`): number of bins to split the region into
- cache_if (`callable`): function to check to against a read to determine if it should be cached
- min_bin_size
- filter_if

**Returns**

: :class:`set` of :class:`pysam.AlignedSegment`: set of reads gathered from the region


### BamCache.close()

close the bam file handle

```python
def close(self):
```

