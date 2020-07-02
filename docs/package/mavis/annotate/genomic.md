# mavis.annotate.genomic

## class Template

**inherits** `BioInterval`






## class IntergenicRegion

**inherits** `BioInterval`


### IntergenicRegion.key()

see :func:`structural_variant.annotate.base.BioInterval.key`

```python
def key(self):
```

### IntergenicRegion.chr()

returns the name of the chromosome that this region resides on

```python
@property
def chr(self):
```

**Args**

- self


### IntergenicRegion.to\_dict()

see :func:`structural_variant.annotate.base.BioInterval.to_dict`

```python
def to_dict(self):
```


## class Gene

**inherits** `BioInterval`


### Gene.transcript\_priority()

prioritizes transcripts from 0 to n-1 based on best transcript flag
and then alphanumeric name sort

Warning:
Lower number means higher priority. This is to make sort work by default

```python
def transcript_priority(self, transcript):
```

**Args**

- transcript

### Gene.transcripts()

:any:`list` of :class:`PreTranscript`: list of unspliced transcripts

```python
@property
def transcripts(self):
```

**Args**

- self

### Gene.translations()

:any:`list` of :class:`~mavis.annotate.protein.Translation`: list of translations

```python
@property
def translations(self):
```

**Args**

- self

### Gene.chr()

returns the name of the chromosome that this gene resides on

```python
@property
def chr(self):
```

**Args**

- self

### Gene.key()

see :func:`structural_variant.annotate.base.BioInterval.key`

```python
def key(self):
```

### Gene.get\_seq()

gene sequence is always given wrt to the positive forward strand regardless of gene strand

```python
def get_seq(self, reference_genome, ignore_cache=False):
```

**Args**

- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence by
- ignore_cache (`bool`): if True then stored sequences will be ignored and the function will attempt to retrieve the sequence using the positions and the input reference_genome

**Returns**

- `str`: the sequence of the gene

### Gene.spliced\_transcripts()

:any:`list` of :class:`Transcript`: list of transcripts

```python
@property
def spliced_transcripts(self):
```

**Args**

- self

### Gene.to\_dict()

see :func:`structural_variant.annotate.base.BioInterval.to_dict`

```python
def to_dict(self):
```


## class Exon

**inherits** `BioInterval`


### Exon.transcript()

:class:`PreTranscript`: the transcript this exon belongs to

```python
@property
def transcript(self):
```

**Args**

- self

### Exon.donor\_splice\_site()

:class:`~mavis.interval.Interval`: the genomic range describing the splice site

```python
@property
def donor_splice_site(self):
```

**Args**

- self

### Exon.acceptor\_splice\_site()

:class:`~mavis.interval.Interval`: the genomic range describing the splice site

```python
@property
def acceptor_splice_site(self):
```

**Args**

- self

### Exon.donor()

`int`: returns the genomic exonic position of the donor splice site

```python
@property
def donor(self):
```

**Args**

- self

### Exon.acceptor()

`int`: returns the genomic exonic position of the acceptor splice site

```python
@property
def acceptor(self):
```

**Args**

- self



## class PreTranscript

**inherits** `BioInterval`

### PreTranscript.\_\_init\_\_()

creates a new transcript object

```python
def __init__(
    self,
    exons,
    gene=None,
    name=None,
    strand=None,
    spliced_transcripts=None,
    seq=None,
    is_best_transcript=False,
):
```

**Args**

- exons (`:class:`list` of :any:`Exon``): list of Exon that make up the transcript
- gene (`Gene`): the gene this transcript belongs to
- name (`str`): name of the transcript
- strand (`STRAND`): strand the transcript is on, defaults to the strand of the Gene if not specified
- spliced_transcripts
- seq (`str`): unspliced cDNA seq
- is_best_transcript

### PreTranscript.generate\_splicing\_patterns()

returns a list of splice sites to be connected as a splicing pattern

```python
def generate_splicing_patterns(self):
```

**Returns**

: :class:`list` of :class:`SplicingPattern`: List of positions to be spliced together see :ref:`theory - predicting splicing patterns <theory-predicting-splicing-patterns>`

### PreTranscript.gene()

:any:`Gene`: the gene this transcript belongs to

```python
@property
def gene(self):
```

**Args**

- self




### PreTranscript.convert\_genomic\_to\_nearest\_cdna()

converts a genomic position to its cdna equivalent or (if intronic) the nearest cdna and shift

```python
def convert_genomic_to_nearest_cdna(
    self, pos, splicing_pattern, stick_direction=None, allow_outside=True
):
```

**Args**

- pos (`int`): the genomic position
- splicing_pattern (`SplicingPattern`): the splicing pattern
- stick_direction
- allow_outside

**Returns**

- `tuple of int and int`:  * *int* - the exonic cdna position * *int* - the intronic shift


### PreTranscript.exon\_number()

exon numbering is based on the direction of translation

```python
def exon_number(self, exon):
```

**Args**

- exon (`Exon`): the exon to be numbered

**Returns**

- `int`: the exon number (1 based)

**Raises**

- `AttributeError`: if the strand is not given or the exon does not belong to the transcript



### PreTranscript.translations()

:class:`list` of :class:`~mavis.annotate.protein.Translation`: list of translations associated with this transcript

```python
@property
def translations(self):
```

**Args**

- self

### PreTranscript.transcripts()

:class:`list` of :class:`Transcript`: list of spliced transcripts

```python
@property
def transcripts(self):
```

**Args**

- self


## class Transcript

**inherits** `BioInterval`

### Transcript.\_\_init\_\_()

splicing pattern is given in genomic coordinates

```python
def __init__(self, pre_transcript, splicing_patt, seq=None, translations=None):
```

**Args**

- pre_transcript (`PreTranscript`): the unspliced transcript
- splicing_patt (`:class:`list` of :class:`int``): the list of splicing positions
- seq (`str`): the cdna sequence






### Transcript.unspliced\_transcript()

:class:`PreTranscript`: the unspliced transcript this splice variant belongs to

```python
@property
def unspliced_transcript(self):
```

**Args**

- self

