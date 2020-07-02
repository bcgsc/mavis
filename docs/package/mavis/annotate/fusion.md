# mavis.annotate.fusion

## class FusionTranscript

**inherits** `PreTranscript`

FusionTranscript is a PreTranscript built from two parent PreTranscripts. It has most of the
same functionality as a regular PreTranscript except that it will not have a parent gene and
retains a mapping of the new exons to the exons in the PreTranscript they originated from

Additionally the FusionTranscript is always constructed on the positive strand.

The preferred way to construct a FusionTranscript is through the build method.

### FusionTranscript.\_\_init\_\_()

```python
def __init__(self):
```










## determine\_prime()

determine the side of the transcript 5' or 3' which is 'kept' given the breakpoint

```python
def determine_prime(transcript, breakpoint):
```

**Args**

- transcript (`Transcript`): the transcript
- breakpoint (`Breakpoint`): the breakpoint

**Returns**

- `PRIME`: 5' or 3'

**Raises**

- `AttributeError`: if the orientation of the breakpoint or the strand of the transcript is not specified
