# mavis.breakpoint

## class Breakpoint

**inherits** `Interval`

class for storing information about a SV breakpoint
coordinates are given as 1-indexed


### Breakpoint.\_\_init\_\_()

```python
def __init__(self, chr, start, end=None, orient=ORIENT.NS, strand=STRAND.NS, seq=None):
```

**Args**

- chr (`str`): the chromosome
- start (`int`): the genomic position of the breakpoint
- end (`int`): if the breakpoint is uncertain (a range) then specify the end of the range here
- orient (`ORIENT`): the orientation (which side is retained at the break)
- strand (`STRAND`): the strand
- seq (`str`): the seq






## class BreakpointPair






### BreakpointPair.interchromosomal()

:class:`bool`: True if the breakpoints are on different chromosomes, False otherwise

```python
@property
def interchromosomal(self):
```

**Args**

- self








### BreakpointPair.flatten()

returns the key-value self for the breakpoint self information as
can be written directly as a tab row

```python
def flatten(self):
```

### BreakpointPair.classify()

uses the chr, orientations and strands to determine the
possible structural_variant types that this pair could support

```python
@classmethod
def classify(cls, pair, distance=None):
```

**Args**

- pair (`BreakpointPair`): the pair to classify
- distance (`callable`): if defined, will be passed to net size to use in narrowing the list of putative types (del vs ins)

**Returns**

: :class:`list` of :any:`SVTYPE`: a list of possible SVTYPE

**Examples**

```python
bpp = BreakpointPair(Breakpoint('1', 1), Breakpoint('1', 9999), opposing_strands=True)
BreakpointPair.classify(bpp)
['inversion']
bpp = BreakpointPair(Breakpoint('1', 1, orient='L'), Breakpoint('1', 9999, orient='R'), opposing_strands=False)
BreakpointPair.classify(bpp)
{'deletion', 'insertion'}
see :ref:`related theory documentation <theory-classifying-events>`
```


### BreakpointPair.net\_size()

Returns the size of the event for a given pair. Mainly applicable to indels

```python
def net_size(self, distance=lambda x, y: Interval(abs(x - y))):
```

**Args**

- distance


### BreakpointPair.breakpoint\_sequence\_homology()

for a given set of breakpoints matches the sequence opposite the partner breakpoint
this sequence comparison is done with reference to a reference genome and does not
use novel or untemplated sequence in the comparison. For this reason, insertions
will never return any homologous sequence

::

small duplication event CTT => CTTCTT

GATACATTTCTTCTTGAAAA reference
## # ---------<========== first breakpoint>-------- second breakpoint## CT-CT------ first break homologyTT-TT-------- second break homology

```python
def breakpoint_sequence_homology(self, reference_genome):
```

**Args**

- reference_genome (`:class:`dict` of :class:`Bio.SeqRecord` by :class:`str``): dict of reference sequence by template/chr name

**Returns**

- `tuple`:  - :class:`str` - homologous sequence at the first breakpoint - :class:`str` - homologous sequence at the second breakpoint

**Raises**

- `AttributeError`: for non specific breakpoints


### BreakpointPair.untemplated\_shift()

gives a range for each breakpoint on the possible alignment range in the shifting the untemplated
sequence

```python
def untemplated_shift(self, reference_genome):
```

**Args**

- reference_genome
