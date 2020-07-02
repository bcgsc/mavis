# mavis.bam.stats

## os.environ[OMP_NUM_THREADS]

```python
os.environ["OMP_NUM_THREADS"] = "4"  # export OMP_NUM_THREADS=4
```

## os.environ[OPENBLAS_NUM_THREADS]

```python
os.environ["OPENBLAS_NUM_THREADS"] = "4"  # export OPENBLAS_NUM_THREADS=4
```

## os.environ[MKL_NUM_THREADS]

```python
os.environ["MKL_NUM_THREADS"] = "4"  # export MKL_NUM_THREADS=6
```

## os.environ[VECLIB_MAXIMUM_THREADS]

```python
os.environ["VECLIB_MAXIMUM_THREADS"] = "4"  # export VECLIB_MAXIMUM_THREADS=4
```

## os.environ[NUMEXPR_NUM_THREADS]

```python
os.environ["NUMEXPR_NUM_THREADS"] = "6"  # export NUMEXPR_NUM_THREADS=6
```

## class BamStats





## class Histogram

**inherits** `dict`

### Histogram.add()

add a key to the histogram with a default frequency of 1

```python
def add(self, item, freq=1):
```

**Args**

- item
- freq

### Histogram.median()

flattens the histogram to compute the median value

```python
def median(self):
```




## compute\_transcriptome\_bam\_stats()

computes various statistical measures relating the input bam file

```python
def compute_transcriptome_bam_stats(
    bam_cache,
    annotations,
    sample_size,
    min_mapping_quality=1,
    stranded=True,
    sample_cap=10000,
    distribution_fraction=0.97,
):
```

**Args**

- bam_cache
- annotations (`object`): see :func:`~mavis.annotate.load_reference_genes`
- sample_size (`int`): the number of genes to compute stats over
- min_mapping_quality (`int`): the minimum mapping quality for a read to be used
- stranded (`bool`): if True then reads must match the gene strand
- sample_cap (`int`): maximum number of reads to collect for any given sample region
- distribution_fraction (`float`): the proportion of the distribution to use in computing stdev

**Returns**

- `BamStats`: the fragment size median, stdev and the read length in a object

## compute\_genome\_bam\_stats()

computes various statistical measures relating the input bam file

```python
def compute_genome_bam_stats(
    bam_file_handle,
    sample_bin_size,
    sample_size,
    min_mapping_quality=1,
    sample_cap=10000,
    distribution_fraction=0.99,
):
```

**Args**

- bam_file_handle (`pysam.AlignmentFile`): the input bam file handle
- sample_bin_size (`int`): how large to make the sample bin (in bp)
- sample_size (`int`): the number of genes to compute stats over
- min_mapping_quality (`int`): the minimum mapping quality for a read to be used
- sample_cap (`int`): maximum number of reads to collect for any given sample region
- distribution_fraction (`float`): the proportion of the distribution to use in computing stdev

**Returns**

- `BamStats`: the fragment size median, stdev and the read length in a object
