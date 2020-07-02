# mavis.annotate.variant

## class Annotation

**inherits** `BreakpointPair`

a fusion of two transcripts created by the associated breakpoint_pair
will also hold the other annotations for overlapping and encompassed and nearest genes

### Annotation.\_\_init\_\_()

Holds a breakpoint call and a set of transcripts, other information is gathered relative to these

```python
def __init__(
    self, bpp, transcript1=None, transcript2=None, proximity=5000, data=None, **kwargs
):
```

**Args**

- bpp (`BreakpointPair`): the breakpoint pair call. Will be adjusted and then stored based on the transcripts
- transcript1 (`Transcript`): transcript at the first breakpoint
- transcript2 (`Transcript`): Transcript at the second breakpoint
- proximity
- data (`dict`): optional dictionary to hold related attributes

### Annotation.add\_gene()

adds a input_gene to the current set of annotations. Checks which set it should be added to

```python
def add_gene(self, input_gene):
```

**Args**

- input_gene (`input_gene`): the input_gene being added

### Annotation.flatten()

generates a dictionary of the annotation information as strings

```python
def flatten(self):
```

**Returns**

: :class:`dict` of :class:`str` by :class:`str`: dictionary of attribute names and values



## class IndelCall

### IndelCall.\_\_init\_\_()

Given two sequences, Assuming there exists a single difference between the two
call an indel which accounts for the change

```python
def __init__(self, refseq, mutseq):
```

**Args**

- refseq (`str`): The reference (amino acid) sequence
- mutseq (`str`): The mutated (amino acid) sequence

### IndelCall.hgvs\_protein\_notation()

returns the HGVS protein notation for an indel call

```python
def hgvs_protein_notation(self):
```



## flatten\_fusion\_translation()

for a given fusion product (translation) gather the information to be output to the tabbed files

```python
def flatten_fusion_translation(translation):
```

**Args**

- translation (`Translation`): the translation which is on the fusion transcript

**Returns**

- `dict`: the dictionary of column names to values

## call\_protein\_indel()

compare the fusion protein/aa sequence to the reference protein/aa sequence and
return an hgvs notation indel call

```python
def call_protein_indel(ref_translation, fusion_translation, reference_genome=None):
```

**Args**

- ref_translation (`Translation`): the reference protein/translation
- fusion_translation (`Translation`): the fusion protein/translation
- reference_genome: the reference genome object used to fetch the reference translation AA sequence

**Returns**

- `str`: the :term:`HGVS` protein indel notation





## choose\_more\_annotated()

for a given set of annotations if there are annotations which contain transcripts and
annotations that are simply intergenic regions, discard the intergenic region annotations

similarly if there are annotations where both breakpoints fall in a transcript and
annotations where one or more breakpoints lands in an intergenic region, discard those
that land in the intergenic region

```python
def choose_more_annotated(ann_list):
```

**Args**

- ann_list (`list of :class:`Annotation``): list of input annotations

**Returns**

- `list of `: class:`Annotation`: the filtered list

## choose\_transcripts\_by\_priority()

for each set of annotations with the same combinations of genes, choose the
annotation with the most "best_transcripts" or most "alphanumeric" choices
of transcript. Throw an error if they are identical

```python
def choose_transcripts_by_priority(ann_list):
```

**Args**

- ann_list (`list of :class:`Annotation``): input annotations

**Returns**

- `list of `: class:`Annotation`: the filtered list

