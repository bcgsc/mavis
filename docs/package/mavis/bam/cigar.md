# mavis.bam.cigar

holds methods related to processing cigar tuples. Cigar tuples are generally
an iterable list of tuples where the first element in each tuple is the
CIGAR value (i.e. 1 for an insertion), and the second value is the frequency

## EVENT_STATES

```python
EVENT_STATES = {CIGAR.D, CIGAR.I, CIGAR.X}
```

## ALIGNED_STATES

```python
ALIGNED_STATES = {CIGAR.M, CIGAR.X, CIGAR.EQ}
```

## REFERENCE_ALIGNED_STATES

```python
REFERENCE_ALIGNED_STATES = ALIGNED_STATES | {CIGAR.D, CIGAR.N}
```

## QUERY_ALIGNED_STATES

```python
QUERY_ALIGNED_STATES = ALIGNED_STATES | {CIGAR.I, CIGAR.S}
```

## CLIPPING_STATE

```python
CLIPPING_STATE = {CIGAR.S, CIGAR.H}
```

## recompute\_cigar\_mismatch()

for cigar tuples where M is used, recompute to replace with X/= for increased
utility and specificity

```python
def recompute_cigar_mismatch(read, ref):
```

**Args**

- read (`pysam.AlignedSegment`): the input read
- ref (`str`): the reference sequence

**Returns**

: :class:`list` of :class:`tuple` of :class:`int` and :class:`int`: the cigar tuple

## longest\_fuzzy\_match()

computes the longest sequence of exact matches allowing for 'x' event interrupts

```python
def longest_fuzzy_match(cigar, max_fuzzy_interupt=1):
```

**Args**

- cigar: cigar tuples
- max_fuzzy_interupt (`int`): number of mismatches allowed

## longest\_exact\_match()

returns the longest consecutive exact match

```python
def longest_exact_match(cigar):
```

**Args**

- cigar (`:class:`list` of :class:`tuple` of :class:`int` and :class:`int``): the cigar tuples

## score()

scoring based on sw alignment properties with gap extension penalties

```python
def score(cigar, **kwargs):
```

**Args**



**Returns**

- `int`: the score value

## match\_percent()

calculates the percent of aligned bases (matches or mismatches) that are matches

```python
def match_percent(cigar):
```

**Args**

- cigar

## join()

given a number of cigar lists, joins them and merges any consecutive tuples
with the same cigar value

```python
def join(*pos):
```

**Examples**

```python
join([(1, 1), (4, 7)], [(4, 3), (2, 4)])
[(1, 1), (4, 10), (2, 4)]
```


## extend\_softclipping()

given some input cigar, extends softclipping if there are mismatches/insertions/deletions
close to the end of the aligned portion. The stopping point is defined by the
min_exact_to_stop_softclipping parameter. this function will throw an error if there is no
exact match aligned portion to signal stop

```python
def extend_softclipping(cigar, min_exact_to_stop_softclipping):
```

**Args**

- cigar
- min_exact_to_stop_softclipping (`int`): number of exact matches to terminate extension

**Returns**

- `tuple`:  - :class:`list` of :class:`~mavis.constants.CIGAR` and :class:`int` - new cigar list - :class:`int` - shift from the original start position

## compute()

given a ref and alt sequence compute the cigar string representing the alt

returns the cigar tuples along with the start position of the alt relative to the ref

```python
def compute(ref, alt, force_softclipping=True, min_exact_to_stop_softclipping=6):
```

**Args**

- ref
- alt
- force_softclipping
- min_exact_to_stop_softclipping

## convert\_for\_igv()

igv does not support the extended CIGAR values for match v mismatch

```python
def convert_for_igv(cigar):
```

**Args**

- cigar

**Examples**

```python
convert_for_igv([(7, 4), (8, 1), (7, 5)])
[(0, 10)]
```


## alignment\_matches()

counts the number of aligned bases irrespective of match/mismatch
this is equivalent to counting all CIGAR.M

```python
def alignment_matches(cigar):
```

**Args**

- cigar

## merge\_indels()

For a given cigar tuple, merges adjacent insertions/deletions

```python
def merge_indels(cigar):
```

**Args**

- cigar

**Examples**

```python
merge_indels([(CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.D, 4), (CIGAR.I, 2), (CIGAR.D, 2), (CIGAR.EQ, 10)])
[(CIGAR.EQ, 10), (CIGAR.I, 5), (CIGAR.D, 6), (CIGAR.EQ, 10)]
```


## hgvs\_standardize\_cigar()

extend alignments as long as matches are possible.
call insertions before deletions

```python
def hgvs_standardize_cigar(read, reference_seq):
```

**Args**

- read
- reference_seq

## convert\_string\_to\_cigar()

Given a cigar string, converts it to the appropriate cigar tuple

```python
def convert_string_to_cigar(string):
```

**Args**

- string

**Examples**

```python
convert_string_to_cigar('8M2I1D9X')
[(CIGAR.M, 8), (CIGAR.I, 2), (CIGAR.D, 1), (CIGAR.X, 9)]
```



## merge\_internal\_events()

merges events (insertions, deletions, mismatches) within a cigar if they are
between exact matches on either side (anchors) and separated by less exact
matches than the given parameter

does not merge two mismatches, must contain a deletion/insertion

```python
def merge_internal_events(cigar, inner_anchor=10, outer_anchor=10):
```

**Args**

- cigar (`list`): a list of tuples of cigar states and counts
- inner_anchor (`int`): minimum number of consecutive exact matches separating events
- outer_anchor (`int`): minimum consecutively aligned exact matches to anchor an end for merging

**Returns**

- `list`: new list of cigar tuples with merged events

**Examples**

```python
merge_internal_events([(CIGAR.EQ, 10), (CIGAR.X, 1), (CIGAR.EQ, 2), (CIGAR.D, 1), (CIGAR.EQ, 10)])
[(CIGAR.EQ, 10), (CIGAR.I, 3), (CIGAR.D, 4), (CIGAR.EQ, 10)]
```
