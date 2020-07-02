# mavis.blat

::

    In general the coordinates in psl files are “zero based half open.” The first base in a sequence is numbered
    zero rather than one. When representing a range the end coordinate is not included in the range. Thus the first
    100 bases of a sequence are represented as 0-100, and the second 100 bases are represented as 100-200. There is
    another little unusual feature in the .psl format. It has to do with how coordinates are handled on the
    negative strand. In the qStart/qEnd fields the coordinates are where it matches from the point of view of the forward
    strand (even when the match is on the reverse strand). However on the qStarts[] list, the coordinates are reversed.
##  http://wiki.bits.vib.be/index.php/Blat

## class Blat

### Blat.millibad()

this function is used in calculating percent identity
direct translation of the perl code
# https://genome.ucsc.edu/FAQ/FAQblat.html#blat4

```python
@staticmethod
def millibad(row, is_protein=False, is_mrna=True):
```

**Args**

- row
- is_protein
- is_mrna

### Blat.score()

direct translation from ucsc guidelines on replicating the web blat score
https://genome.ucsc.edu/FAQ/FAQblat.html#blat4

below are lines from the perl code i've re-written in python

::

my $sizeMul = pslIsProtein($blockCount, $strand, $tStart, $tEnd, $tSize, $tStarts, $blockSizes);
sizmul = 1 for DNA
my $pslScore = $sizeMul * ($matches + ($repMatches >> 1) ) - $sizeMul * $misMatches - $qNumInsert - $tNumIns
ert)

```python
@staticmethod
def score(row, is_protein=False):
```

**Args**

- row
- is_protein



### Blat.pslx\_row\_to\_pysam()

given a 'row' from reading a pslx file. converts the row to a BlatAlignedSegment object

```python
@staticmethod
def pslx_row_to_pysam(row, bam_cache, reference_genome):
```

**Args**

- row (`dict of str`): a row object from the 'read_pslx' method
- bam_cache (`BamCache`): the bam file/cache to use as a template for creating reference_id from chr name



## process\_blat\_output()

converts the blat output pslx (unheadered file) to bam reads

```python
def process_blat_output(
    input_bam_cache,
    query_id_mapping,
    reference_genome,
    aligner_output_file='aligner_out.temp',
    blat_min_percent_of_max_score=0.8,
    blat_min_identity=0.7,
    blat_limit_top_aln=25,
    is_protein=False,
):
```

**Args**

- input_bam_cache
- query_id_mapping
- reference_genome
- aligner_output_file
- blat_min_percent_of_max_score
- blat_min_identity
- blat_limit_top_aln
- is_protein
