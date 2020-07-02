# mavis.bam.read

## class SamRead

**inherits** `pysam.AlignedSegment`

Subclass to extend the pysam.AlignedSegment class adding some utility methods and convenient representations

Allows next_reference_name and reference_name to be set directly so that is does not depend on a bam header

### SamRead.\_\_init\_\_()

```python
def __init__(
    self, reference_name=None, next_reference_name=None, alignment_score=None, **kwargs
):
```

**Args**

- reference_name
- next_reference_name
- alignment_score

### SamRead.set\_key()

Warning:
Using this method sets the _key attribute which is used for comparison and hashing. If you alter
this attribute while items are in a hashed state it may lead to unexpected results such as duplicates
of a single object within a set

```python
def set_key(self):
```

### SamRead.key()

uses a stored _key attribute, if available. This is to avoid the hash changing if the reference start (for example)
is changed but also allow this attribute to be used and calculated for non SamRead objects

This way to change the hash behaviour the user must be explicit and use the set_key method

```python
def key(self):
```







### SamRead.deletion\_sequences()

returns the reference sequences for all deletions

```python
def deletion_sequences(self, reference_genome):
```

**Args**

- reference_genome

### SamRead.insertion\_sequences()

returns the inserted sequence for all insertions

```python
def insertion_sequences(self):
```




## pileup()

For a given set of reads generate a pileup of all reads (excluding those for which the filter_func returns True)

```python
def pileup(reads, filter_func=None):
```

**Args**

- reads (`iterable of pysam.AlignedSegment`): reads to pileup
- filter_func (`callable`): function which takes in a  read and returns True if it should be ignored and False otherwise

**Returns**

- `iterable of tuple of int and int`: tuples of genomic position and read count at that position


## breakpoint\_pos()

assumes the breakpoint is the position following softclipping on the side with more
softclipping (unless and orientation has been specified)

```python
def breakpoint_pos(read, orient=ORIENT.NS):
```

**Args**

- read (`:class:`~pysam.AlignedSegment``): the read object
- orient (`ORIENT`): the orientation

**Returns**

- `int`: the position of the breakpoint in the input read

## calculate\_alignment\_score()

calculates a score for comparing alignments

```python
def calculate_alignment_score(read, consec_bonus=1):
```

**Args**

- read (`pysam.AlignedSegment`): the input read
- consec_bonus

**Returns**

- `float`: the score

## nsb\_align()

given some reference string and a smaller sequence string computes the best non-space-breaking alignment
i.e. an alignment that does not allow for indels (straight-match). Positions in the aligned segments are
given relative to the length of the reference sequence (1-based)

```python
def nsb_align(
    ref,
    seq,
    weight_of_score=0.5,
    min_overlap_percent=1,
    min_match=0,
    min_consecutive_match=1,
    scoring_function=calculate_alignment_score,
):
```

**Args**

- ref (`str`): the reference sequence
- seq (`str`): the sequence being aligned
- weight_of_score (`float`): when scoring alignments this determines the amount
- min_overlap_percent (`float`): the minimum amount of overlap of the input sequence to the reference
- min_match (`float`): the minimum number of matches compared to total
- min_consecutive_match
- scoring_function (`callable`): any function that will take a read as input and return a float

**Returns**

: :class:`list` of :class:`~pysam.AlignedSegment`: list of aligned segments

## sequenced\_strand()

determines the strand that was sequenced

```python
def sequenced_strand(read, strand_determining_read=2):
```

**Args**

- read (`:class:`~pysam.AlignedSegment``): the read being used to determine the strand
- strand_determining_read (`int`): which read in the read pair is the same as the sequenced strand

**Returns**

- `STRAND`: the strand that was sequenced

**Raises**

- `ValueError`: if strand_determining_read is not 1 or 2
- `Warning`
: if the input pair is unstranded the information will not be representative of the
: strand sequenced since the assumed convention is not followed

## read\_pair\_type()

assumptions based on illumina pairs: only 4 possible combinations

```python
def read_pair_type(read):
```

**Args**

- read (`:class:`~pysam.AlignedSegment``): the input read

**Returns**

- `READ_PAIR_TYPE`: the type of input read pair

**Raises**

- `NotImplementedError`: for any read that does not fall into the four expected configurations (see below)
: ::
: ++++> <---- is LR same-strand
: ++++> ++++> is LL opposite
: <---- <---- is RR opposite
: <---- ++++> is RL same-strand

## orientation\_supports\_type()

checks if the orientation is compatible with the type of event

```python
def orientation_supports_type(read, event_type):
```

**Args**

- read (`:class:`~pysam.AlignedSegment``): a read from the pair
- event_type (`SVTYPE`): the type of event to check

**Returns**

- `bool`:  - ``True`` - the read pair is in the correct orientation for this event type - ``False`` - the read is not in the correct orientation

## convert\_events\_to\_softclipping()

given an alignment, simplifies the alignment by grouping everything past the first anchor and including the
first event considered too large and unaligning them turning them into softclipping

```python
def convert_events_to_softclipping(read, orientation, max_event_size, min_anchor_size=None):
```

**Args**

- read
- orientation
- max_event_size
- min_anchor_size

## sequence\_complexity()

basic measure of sequence complexity

```python
def sequence_complexity(seq):
```

**Args**

- seq
