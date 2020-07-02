# mavis.annotate.protein

## class DomainRegion

**inherits** `BioInterval`



## class Domain


### Domain.translation()

:class:`~mavis.annotate.Translation`: the Translation this domain belongs to

```python
@property
def translation(self):
```

**Args**

- self

### Domain.key()

:class:`tuple`: a tuple representing the items expected to be unique. for hashing and comparing

```python
def key(self):
```

### Domain.score\_region\_mapping()

compares the sequence in each DomainRegion to the sequence collected for that domain region from the
translation object

```python
def score_region_mapping(self, reference_genome=None):
```

**Args**

- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence

**Returns**

- `tuple of int and int`: tuple contains - int: the number of matching amino acids - int: the total number of amino acids

### Domain.get\_seqs()

returns the amino acid sequences for each of the domain regions associated with
this domain in the order of the regions (sorted by start)

```python
def get_seqs(self, reference_genome=None, ignore_cache=False):
```

**Args**

- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence
- ignore_cache

**Returns**

: :class:`list` of :class:`str`: list of amino acid sequences for each DomainRegion

**Raises**

- `AttributeError`: if there is not enough sequence information given to determine this

### Domain.align\_seq()

align each region to the input sequence starting with the last one.
then take the subset of sequence that remains to align the second last and so on
return a list of intervals for the alignment. If multiple alignments are found,
then raise an error

```python
def align_seq(self, input_sequence, reference_genome=None, min_region_match=0.5):
```

**Args**

- input_sequence (`str`): the sequence to be aligned to
- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence
- min_region_match (`float`): percent between 0 and 1. Each region must have a score len(seq) * min_region_match

**Returns**

- `tuple`: tuple contains - int: the number of matches - int: the total number of amino acids to be aligned - :any:`list` of :any:`DomainRegion`: the list of domain regions on the new input sequence

**Raises**

- `AttributeError`: if sequence information is not available
- `UserWarning`: if a valid alignment could not be found or no best alignment was found


## class Translation

**inherits** `BioInterval`

### Translation.\_\_init\_\_()

describes the splicing pattern and cds start and end with reference to a particular transcript

```python
def __init__(self, start, end, transcript=None, domains=None, seq=None, name=None):
```

**Args**

- start (`int`): start of the coding sequence (cds) relative to the start of the first exon in the transcript
- end (`int`): end of the coding sequence (cds) relative to the start of the first exon in the transcript
- transcript (`Transcript`): the transcript this is a Translation of
- domains (`:class:`list` of :any:`Domain``): a list of the domains on this translation
- seq
- name

### Translation.transcript()

:class:`~mavis.annotate.genomic.Transcript`: the spliced transcript this translation belongs to

```python
@property
def transcript(self):
```

**Args**

- self



### Translation.convert\_genomic\_to\_cds()

converts a genomic position to its cds (coding sequence) equivalent

```python
def convert_genomic_to_cds(self, pos):
```

**Args**

- pos (`int`): the genomic position

**Returns**

- `int`: the cds position (negative if before the initiation start site)

### Translation.convert\_genomic\_to\_nearest\_cds()

converts a genomic position to its cds equivalent or (if intronic) the nearest cds and shift

```python
def convert_genomic_to_nearest_cds(self, pos):
```

**Args**

- pos (`int`): the genomic position

**Returns**

- `tuple of int and int`:  * *int* - the cds position * *int* - the intronic shift

### Translation.convert\_genomic\_to\_cds\_notation()

converts a genomic position to its cds (coding sequence) equivalent using
`hgvs <http://www.hgvs.org/mutnomen/recs-DNA.html>`_ cds notation

```python
def convert_genomic_to_cds_notation(self, pos):
```

**Args**

- pos (`int`): the genomic position

**Returns**

- `str`: the cds position notation

**Examples**

```python
tl = Translation(...)
# a position before the translation start
tl.convert_genomic_to_cds_notation(1010)
'-50'
# a position after the translation end
tl.convert_genomic_to_cds_notation(2031)
'*72'
# an intronic position
tl.convert_genomic_to_cds_notation(1542)
'50+10'
tl.convert_genomic_to_cds_notation(1589)
'51-14'
```



### Translation.get\_seq()

wrapper for the sequence method

```python
def get_seq(self, reference_genome=None, ignore_cache=False):
```

**Args**

- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence
- ignore_cache


### Translation.key()

see :func:`structural_variant.annotate.base.BioInterval.key`

```python
def key(self):
```


## calculate\_orf()

calculate all possible open reading frames given a spliced cdna sequence (no introns)

```python
def calculate_orf(spliced_cdna_sequence, min_orf_size=None):
```

**Args**

- spliced_cdna_sequence (`str`): the sequence
- min_orf_size

**Returns**

: :any:`list` of :any:`Interval`: list of open reading frame positions on the input sequence
