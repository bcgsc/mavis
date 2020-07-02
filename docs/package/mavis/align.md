# mavis.align

Should take in a sam file from a aligner like bwa aln or bwa mem and convert it into a

## SUPPORTED_ALIGNER

```python
SUPPORTED_ALIGNER = MavisNamespace(
    BWA_MEM='bwa mem', BLAT='blat', __name__='~mavis.align.SUPPORTED_ALIGNER'
)
```

## class SplitAlignment

**inherits** `BreakpointPair`




### SplitAlignment.query\_coverage()

interval representing the total region of the input sequence that is covered by the combination of alignments

```python
def query_coverage(self):
```

### SplitAlignment.query\_consumption()

fraction of the query sequence which is aligned (everything not soft-clipped) in either alignment

```python
def query_consumption(self):
```


### SplitAlignment.score()

scores events between 0 and 1 penalizing events interrupting the alignment. Counts a split
alignment as a single event

```python
def score(self, consec_bonus=10):
```

**Args**

- consec_bonus






## get\_aligner\_version()

executes a subprocess to try and run the aligner without arguments and parse the version number from the output

```python
def get_aligner_version(aligner):
```

**Args**

- aligner

**Examples**

```python
get_aligner_version('blat')
'36x2'
```



## convert\_to\_duplication()

Given a breakpoint call, tests if the untemplated sequences matches the preceding
reference sequence. If it does this is annotated as a duplication and the new
breakpoint pair is returned. If not, then the original breakpoint pair is returned

```python
def convert_to_duplication(alignment, reference_genome):
```

**Args**

- alignment
- reference_genome

## call\_read\_events()

Given a read, return breakpoint pairs representing all putative events

```python
def call_read_events(read, secondary_read=None, is_stranded=False):
```

**Args**

- read
- secondary_read
- is_stranded

## read\_breakpoint()

convert a given read to a single breakpoint

```python
def read_breakpoint(read):
```

**Args**

- read

## call\_paired\_read\_event()

For a given pair of reads call all applicable events. Assume there is a major
event from both reads and then call indels from the individual reads

```python
def call_paired_read_event(read1, read2, is_stranded=False):
```

**Args**

- read1
- read2
- is_stranded

## align\_sequences()

calls the alignment tool and parses the return output for a set of sequences

```python
def align_sequences(
    sequences,
    input_bam_cache,
    reference_genome,
    aligner,
    aligner_reference,
    aligner_output_file='aligner_out.temp',
    aligner_fa_input_file='aligner_in.fa',
    aligner_output_log='aligner_out.log',
    blat_limit_top_aln=25,
    blat_min_identity=0.7,
    clean_files=True,
    log=DEVNULL,
    **kwargs
):
```

**Args**

- sequences (`dict of str to str`): dictionary of sequences by name
- input_bam_cache (`BamCache`): bam cache to be used as a template for reading the alignments
- reference_genome: the reference genome
- aligner (`SUPPORTED_ALIGNER`): the name of the aligner to be used
- aligner_reference (`str`): path to the aligner reference file
- aligner_output_file
- aligner_fa_input_file
- aligner_output_log
- blat_limit_top_aln
- blat_min_identity
- clean_files
- log

## select\_contig\_alignments()

standardize/simplify reads and filter bad/irrelevant alignments
adds the contig alignments to the contigs

```python
def select_contig_alignments(evidence, reads_by_query):
```

**Args**

- evidence
- reads_by_query
