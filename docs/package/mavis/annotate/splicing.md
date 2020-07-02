# mavis.annotate.splicing

## class SplicingPattern

**inherits** `list`




### SplicingPattern.generate\_patterns()

returns a list of splice sites to be connected as a splicing pattern

```python
@classmethod
def generate_patterns(cls, sites, is_reverse=False):
```

**Args**

- sites
- is_reverse

**Returns**

: :class:`list` of :class:`SplicingPattern`: List of positions to be spliced together see :ref:`theory - predicting splicing patterns <theory-predicting-splicing-patterns>`


## class SpliceSite

**inherits** `BioInterval`





## predict\_splice\_sites()

looks for the expected splice site sequence patterns in the
input strings and returns a list of putative splice sites

```python
def predict_splice_sites(input_sequence, is_reverse=False):
```

**Args**

- input_sequence (`str`): input sequence with respect to the positive/forward strand
- is_reverse (`bool`): True when the sequences is transcribed on the reverse strand
