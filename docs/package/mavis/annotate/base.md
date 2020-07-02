# mavis.annotate.base

## class ReferenceName

**inherits** `str`

Class for reference sequence names. Ensures that hg19/hg38 chromosome names match.









## class BioInterval


### BioInterval.start()

*int*: the start position

```python
@property
def start(self):
```

**Args**

- self

### BioInterval.end()

*int*: the end position

```python
@property
def end(self):
```

**Args**

- self



### BioInterval.key()

:class:`tuple`: a tuple representing the items expected to be unique. for hashing and comparing

```python
def key(self):
```




### BioInterval.get\_seq()

get the sequence for the current annotation object

```python
def get_seq(self, reference_genome=None, ignore_cache=False):
```

**Args**

- reference_genome
- ignore_cache

**Raises**

- `NotImplementedError`: abstract method

### BioInterval.get\_strand()

pulls strand information from the current object, or follows reference
objects until the strand is found

```python
def get_strand(self):
```

**Returns**

- `STRAND`: the strand of this or any of its reference objects

**Raises**

- `AttributeError`: raised if the strand is not set on this or any of its reference objects

### BioInterval.is\_reverse()

True if the gene is on the reverse/negative strand.

```python
@property
def is_reverse(self):
```

**Args**

- self

**Raises**

- `AttributeError`: if the strand is not specified

### BioInterval.get\_chr()

pulls chromosome information from the current object, or follows reference
objects until the chromosome is found

```python
def get_chr(self):
```

**Returns**

- `str`: the chromosome of this or any of its reference objects

**Raises**

- `AttributeError`: raised if the chromosome is not set on this or any of its reference objects

### BioInterval.to\_dict()

creates a dictionary representing the current object

```python
def to_dict(self):
```

**Returns**

: :class:`dict` by :class:`str`: the dictionary of attribute values


