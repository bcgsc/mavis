# Column Names

List of column names and their definitions. The types indicated here are
the expected types in a row for a given column name.


## library

Identifier for the library/source

## cluster\_id

Identifier for the merging/clustering step

## cluster\_size

**type**: `int`

The number of breakpoint
pair calls that were grouped in creating the cluster

## validation\_id

Identifier for the validation step

## annotation\_id

Identifier for the annotation step

## product\_id

Unique identifier of the final fusion including splicing and ORF
decision from the annotation step

## event\_type

**type**: `mavis.constants.SVTYPE`

The
classification of the event

## inferred\_pairing

A semi colon delimited of event identifiers i.e.
`<annotation_id>_<splicing pattern>_<cds start>_<cds end>`
which were paired to the current event based on predicted products

## pairing

A semi colon delimited of event identifiers i.e.
`<annotation_id>_<splicing pattern>_<cds start>_<cds end>`
which were paired to the current event based on breakpoint positions

## gene1

Gene for the current annotation at the first breakpoint

## gene1\_direction

**type**: `mavis.constants.PRIME`

The
direction/prime of the gene

## gene2

Gene for the current annotation at the second breakpoint

## gene2\_direction

**type**: `mavis.constants.PRIME`

The
direction/prime of the gene. Has the following possible values

## gene1\_aliases

Other gene names associated with the current annotation at the first
breakpoint

## gene2\_aliases

Other gene names associated with the current annotation at the
second breakpoint

## gene\_product\_type

**type**: `mavis.constants.GENE_PRODUCT_TYPE`

Describes if the putative fusion product will be
sense or anti-sense

## fusion\_cdna\_coding\_end

Position wrt the 5' end of the fusion transcript where coding ends
last base of the stop codon

## transcript1

Transcript for the current annotation at the first breakpoint

## transcript2

Transcript for the current annotation at the second breakpoint

## fusion\_splicing\_pattern

`mavis.constants.SPLICE_TYPE` -
Type of splicing pattern used to create the fusion cDNA.

## fusion\_cdna\_coding\_start

**type**: `int`

Position wrt the 5' end of
the fusion transcript where coding begins first base of the Met
amino acid.

## fusion\_cdna\_coding\_end

**type**: `int`

Position wrt the 5' end of
the fusion transcript where coding ends last base of the stop codon

## fusion\_mapped\_domains

**type**: `JSON`

List of domains in [JSON](#JSON)
format where each domain start and end positions are given wrt to
the fusion transcript and the mapping quality is the number of
matching amino acid positions over the total number of amino acids.
The sequence is the amino acid sequence of the domain on the
reference/original transcript

## fusion\_sequence\_fasta\_id

The sequence identifier for the cdna sequence output fasta file

## fusion\_sequence\_fasta\_file

**type**: `FILEPATH`

Path to the corresponding fasta output file

## annotation\_figure

**type**: `FILEPATH`

File path to the svg drawing representing the
annotation

## annotation\_figure\_legend

**type**: `JSON`

[JSON](#JSON) data for the figure
legend

## genes\_encompassed

Applies to intrachromosomal events only. List of genes which overlap
any region that occurs between both breakpoints. For example in a
deletion event these would be deleted genes.

## genes\_overlapping\_break1

list of genes which overlap the first breakpoint

## genes\_overlapping\_break2

list of genes which overlap the second breakpoint

## genes\_proximal\_to\_break1

list of genes near the breakpoint and the distance away from the
breakpoint

## genes\_proximal\_to\_break2

list of genes near the breakpoint and the distance away from the
breakpoint

## break1\_chromosome

**type**: `str`

The name of the chromosome
on which breakpoint 1 is situated

## break1\_position\_start

**type**: `int`

Start integer inclusive
1-based of the range representing breakpoint 1

## break1\_position\_end

**type**: `int`

End integer inclusive
1-based of the range representing breakpoint 1

## break1\_orientation

**type**: `mavis.constants.ORIENT`

The side
of the breakpoint wrt the positive/forward strand that is retained.

## break1\_strand

**type**: `mavis.constants.STRAND`

The
strand wrt to the reference positive/forward strand at this
breakpoint.

## break1\_seq

**type**: `str`

The sequence up to and
including the breakpoint. Always given wrt to the positive/forward
strand

## break2\_chromosome

The name of the chromosome on which breakpoint 2 is situated

## break2\_position\_start

**type**: `int`

Start integer inclusive
1-based of the range representing breakpoint 2

## break2\_position\_end

**type**: `int`

End integer inclusive
1-based of the range representing breakpoint 2

## break2\_orientation

**type**: `mavis.constants.ORIENT`

The side
of the breakpoint wrt the positive/forward strand that is retained.

## break2\_strand

**type**: `mavis.constants.STRAND`

The
strand wrt to the reference positive/forward strand at this
breakpoint.

## break2\_seq

**type**: `str`

The sequence up to and
including the breakpoint. Always given wrt to the positive/forward
strand

## opposing\_strands

**type**: `bool`

Specifies if breakpoints
are on opposite strands wrt to the reference. Expects a boolean

## stranded

**type**: `bool`

Specifies if the sequencing
protocol was strand specific or not. Expects a boolean

## protocol

`mavis.constants.PROTOCOL` -
Specifies the type of library

## tools

The tools that called the event originally from the cluster step.
Should be a semi-colon delimited list of `<tool name>_<tool version>`

## contigs\_assembled

**type**: `int`

Number of contigs that were
built from split read sequences

## contigs\_aligned

**type**: `int`

Number of contigs that were
able to align

## contig\_alignment\_query\_name

The query name for the contig alignment. Should match the 'read'
name(s) in the .contigs.bam output file

## contig\_seq

**type**: `str`

Sequence of the current
contig wrt to the positive forward strand if not strand specific

## contig\_remap\_score

**type**: `float`

Score representing the
number of sequences from the set of sequences given to the assembly
algorithm that were aligned to the resulting contig with an
acceptable scoring based on user-set thresholds. For any sequence
its contribution to the score is divided by the number of mappings
to give less weight to multimaps

## call\_sequence\_complexity

**type**: `float`

The minimum amount any two
bases account for of the proportion of call sequence. An average for
non-contig calls

## contig\_remapped\_reads

**type**: `int`

the number of reads from the
input bam that map to the assembled contig

## contig\_remapped\_read\_names

read query names for the reads that were remapped. A -1 or -2 has
been appended to the end of the name to indicate if this is the
first or second read in the pair

## contig\_alignment\_score

**type**: `float`

A rank based on the
alignment tool blat etc. of the alignment being used. An average if
split alignments were used. Lower numbers indicate a better
alignment. If it was the best alignment possible then this would be
zero.

## contig\_alignment\_reference\_start

The reference start(s) `<chr>:<position>` of the contig alignment.
Semi-colon delimited

## contig\_alignment\_cigar

The cigar string(s) representing the contig alignment. Semi-colon
delimited

## contig\_remap\_coverage

**type**: `float`

Fraction of the contig
sequence which is covered by the remapped reads

## contig\_build\_score

**type**: `int`

Score representing the edge
weights of all edges used in building the sequence

## contig\_strand\_specific

**type**: `bool`

A flag to indicate if it
was possible to resolve the strand for this contig

## spanning\_reads

**type**: `int`

the number of spanning reads
which support the event

## spanning\_read\_names

read query names of the spanning reads which support the current
event

## call\_method

**type**: `mavis.constants.CALL_METHOD`

The
method used to call the breakpoints

## flanking\_pairs

**type**: `int`

Number of read-pairs where
one read aligns to the first breakpoint window and the second read
aligns to the other. The count here is based on the number of unique
query names

## flanking\_pairs\_compatible

**type**: `int`

Number of flanking pairs of
a compatible orientation type. This applies to insertions and
duplications. Flanking pairs supporting an insertion will be
compatible to a duplication and flanking pairs supporting a
duplication will be compatible to an insertion (possibly indicating
an internal translocation)

## flanking\_median\_fragment\_size

**type**: `int`

The median fragment size of
the flanking reads being used as evidence

## flanking\_stdev\_fragment\_size

**type**: `float`

The standard deviation in
fragment size of the flanking reads being used as evidence

## break1\_split\_reads

**type**: `int`

Number of split reads that
call the exact breakpoint given

## break1\_split\_reads\_forced

**type**: `int`

Number of split reads which
were aligned to the opposite breakpoint window using a targeted
alignment

## break2\_split\_reads

**type**: `int`

Number of split reads that
call the exact breakpoint given

## break2\_split\_reads\_forced

**type**: `int`

Number of split reads which
were aligned to the opposite breakpoint window using a targeted
alignment

## linking\_split\_reads

**type**: `int`

Number of split reads that
align to both breakpoints

## untemplated\_seq

**type**: `str`

The untemplated/novel
sequence between the breakpoints

## break1\_homologous\_seq

**type**: `str`

Sequence in common at the
first breakpoint and other side of the second breakpoint

## break2\_homologous\_seq

**type**: `str`

Sequence in common at the
second breakpoint and other side of the first breakpoint

## break1\_ewindow

**type**: `int-int`

Window where evidence was gathered for the first
breakpoint

## break1\_ewindow\_count

**type**: `int`

Number of reads
processed/looked-at in the first evidence window

## break1\_ewindow\_practical\_coverage

```python
break1_ewindow_practical_coverage: float = break1_ewindow_count / len(break1_ewindow)
```

Not the actual coverage as bins are sampled within and there is a read limit cutoff

## break2\_ewindow

**type**: `int-int`

Window where evidence was gathered for the second
breakpoint

## break2\_ewindow\_count

**type**: `int`

Number of reads
processed/looked-at in the second evidence window

## break2\_ewindow\_practical\_coverage

```python
break2_ewindow_practical_coverage: float = break2_ewindow_count / len(break2_ewindow)
```

Not the actual coverage as bins are sampled within and there is a read limit cutoff

## raw\_flanking\_pairs

**type**: `int`

Number of flanking reads
before calling the breakpoint. The count here is based on the number
of unique query names

## raw\_spanning\_reads

**type**: `int`

Number of spanning reads
collected during evidence collection before calling the breakpoint

## raw\_break1\_split\_reads

**type**: `int`

Number of split reads before
calling the breakpoint

## raw\_break2\_split\_reads

**type**: `int`

Number of split reads before
calling the breakpoint

## cdna\_synon

semi-colon delimited list of transcript ids which have an identical
cdna sequence to the cdna sequence of the current fusion product

## protein\_synon

semi-colon delimited list of transcript ids which produce a
translation with an identical amino-acid sequence to the current
fusion product

## tracking\_id

column used to store input identifiers from the original SV calls.
Used to track calls from the input files to the final outputs.

## fusion\_protein\_hgvs

**type**: `str`

Describes the fusion protein
in HGVS notation. Will be None if the change is not an indel or is
synonymous

## net\_size

**type**: `int-int`

The net size of an event. For translocations and
inversion this will always be 0. For indels it will be negative for
deletions and positive for insertions. It is a range to accommodate
non-specific events.

## supplementary\_call

**type**: `bool`

Flag to indicate if the
current event was a supplementary call, meaning a call that was
found as a result of validating another event.
