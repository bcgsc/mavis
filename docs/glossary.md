# Glossary

## General Terms

### flanking read pair

A pair of reads where one read maps to one side of a set of
    breakpoints and its mate maps to the other.

### split read

A read which aligns next to a breakpoint and is softclipped at one
    or more sides.

### spanning read

Applies primarily to small structural variants. Reads which span
    both breakpoints.

### half-mapped read

A read whose mate is unaligned. Generally this refers to reads in
    the evidence stage that are mapped next to a breakpoint.

### breakpoint

A breakpoint is a genomic position (interval) on some
    reference/template/chromosome which has a strand and orientation.
    The orientation describes the portion of the reference that is
    retained.

### event

Used interchangeably with [structural variant](./../glossary#structural variant).

### event type

Classification for a structural variant. see
    [event_type](./../glossary#event_type).

### structural variant

A genomic alteration that can be described by a pair of breakpoints
    and an [event type](./../glossary#event type). The two
    breakpoints represent regions in the genome that are broken apart
    and reattached together.

### breakpoint pair

Basic definition of a [structural variant](./../glossary#structural variant). Does not automatically imply a classification/type.

### bed

File format specification. See
    <https://genome.ucsc.edu/FAQ/FAQformat#format1>.

### IGV batch file

This is a file format type defined by [IGV](./../glossary#IGV) see [running IGV with a batch
    file](https://software.broadinstitute.org/software/igv/batch).

### BAM

File format specification. See
    <https://genome.ucsc.edu/FAQ/FAQformat#format5.1>.

### 2bit

File format specification. See
    <https://genome.ucsc.edu/FAQ/FAQformat#format7>.

### fasta

File format specification. See
    <https://genome.ucsc.edu/FAQ/FAQformat#format18>.

### psl

File format specification. See
    <https://genome.ucsc.edu/FAQ/FAQformat#format2>.

### pslx

Extended format of a [psl](./../glossary#psl).

### SVG

SVG (Scalable vector graph) is an image format. see
    <https://www.w3schools.com/graphics/svg_intro.asp>.

### JSON

JSON (JavaScript Object Notation) is a data file format. see
    <https://www.w3schools.com/js/js_json_intro.asp>.

### blat

Alignment tool. see <https://genome.ucsc.edu/FAQ/FAQblat.html#blat3>
    for instructions on download and install.

### IGV

Integrative Genomics Viewer is a visualization tool. see
    <http://software.broadinstitute.org/software/igv>.

### SGE

Sun Grid Engine (SGE) is a job scheduling system for cluster
    management see
    <http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html>.

### SLURM

SLURM is a job scheduling system for cluster management see
    <https://slurm.schedmd.com/archive/slurm-17.02.1>.

### TORQUE

TORQUE is a job scheduling system for cluster management see
    <http://www.adaptivecomputing.com/products/open-source/torque/>.

### BWA

BWA is an alignement tool. See <https://github.com/lh3/bwa> for
    install instructions.

### `HGVS <den-Dunnen-2016>`{.interpreted-text role="ref"}

Community based standard of reccommendations for variant notation.
    See <http://varnomen.hgvs.org/>

### BreakDancer

BreakDancer is an SV caller. Source for BreakDancer can be found
    [here](https://github.com/genome/breakdancer) [Chen-2009](./citations#Chen-2009)

### BreakSeq

BreakSeq is an SV caller. Source for BreakSeq can be found
    [here](https://github.com/bioinform/breakseq2) [Abyzov-2015](./citations#Abyzov-2015)

### Chimerascan

Chimerascan is an SV caller. Source for Chimerascan can be found
    [here](https://code.google.com/archive/p/chimerascan)
    [Iyer-2011](./citations#Iyer-2011)

### CNVnator

CNVnator is an SV caller. Source for CNVnator can be found
    [here](https://github.com/abyzovlab/CNVnator) [Abyzov-2011](./citations#Abyzov-2011)

### DeFuse

DeFuse is an SV caller. Source for DeFuse can be found
    [here](https://bitbucket.org/dranew/defuse) [McPherson-2011](./citations#McPherson-2011)

### DELLY

DELLY is an SV caller. Source for DELLY can be found
    [here](https://github.com/dellytools/delly) [Rausch-2012](./citations#Rausch-2012)

### Manta

Manta is an SV caller. Source for Manta can be found
    [here](https://github.com/Illumina/manta) [Chen-2016](./citations#Chen-2016)

### STAR-Fusion

STAR-Fusion is an SV caller. Source for STAR-Fusion can be found
    [here](https://github.com/STAR-Fusion/STAR-Fusion) [Haas-2017](./citations#Haas-2017)

### Strelka

Strelka is an SNV and small indel caller. Only the small indels can
    be processed, since SNVs are not currently suported. Source for
    Strelka can be found [here](https://github.com/Illumina/strelka)
    [Saunders-2012](./citations#Saunders-2012)

### Pindel

Pindel is an SV caller. Source for Pindel can be found
    [here](https://github.com/genome/pindel) [Ye-2009](./citations#Ye-2009)

### Trans-ABySS

Trans-ABySS is an SV caller. Source for Trans-ABySS can be found
    [here](https://github.com/bcgsc/transabyss) [Robertson-2010](./citations#Robertson-2010)

### SV

Structural Variant

## Column Names

List of column names and their definitions. The types indicated here are
the expected types in a row for a given column name.


### library

Identifier for the library/source

### cluster\_id

Identifier for the merging/clustering step

### cluster\_size

`int` - The number of breakpoint
    pair calls that were grouped in creating the cluster

### validation\_id

Identifier for the validation step

### annotation\_id

Identifier for the annotation step

### product\_id

Unique identifier of the final fusion including splicing and ORF
    decision from the annotation step

### event\_type

`~mavis.constants.SVTYPE` - The
    classification of the event

### inferred\_pairing

A semi colon delimited of event identifiers i.e.
    \<annotation\_id\>\_\<splicing pattern\>\_\<cds start\>\_\<cds end\>
    which were paired to the current event based on predicted products

### pairing

A semi colon delimited of event identifiers i.e.
    \<annotation\_id\>\_\<splicing pattern\>\_\<cds start\>\_\<cds end\>
    which were paired to the current event based on breakpoint positions

### gene1

Gene for the current annotation at the first breakpoint

### gene1\_direction

`~mavis.constants.PRIME` - The
    direction/prime of the gene

### gene2

Gene for the current annotation at the second breakpoint

### gene2\_direction

`~mavis.constants.PRIME` - The
    direction/prime of the gene. Has the following possible values

### gene1\_aliases

Other gene names associated with the current annotation at the first
    breakpoint

### gene2\_aliases

Other gene names associated with the current annotation at the
    second breakpoint

### gene\_product\_type

`~mavis.constants.GENE_PRODUCT_TYPE`{.interpreted-text
    role="class"} - Describes if the putative fusion product will be
    sense or anti-sense

### fusion\_cdna\_coding\_end

Position wrt the 5\' end of the fusion transcript where coding ends
    last base of the stop codon

### transcript1

Transcript for the current annotation at the first breakpoint

### transcript2

Transcript for the current annotation at the second breakpoint

### fusion\_splicing\_pattern

`~mavis.constants.SPLICE_TYPE` -
    Type of splicing pattern used to create the fusion cDNA.

### fusion\_cdna\_coding\_start

`int` - Position wrt the 5\' end of
    the fusion transcript where coding begins first base of the Met
    amino acid.

### fusion\_cdna\_coding\_end

`int` - Position wrt the 5\' end of
    the fusion transcript where coding ends last base of the stop codon

### fusion\_mapped\_domains

`JSON` - List of domains in [JSON](./../glossary#JSON)
    format where each domain start and end positions are given wrt to
    the fusion transcript and the mapping quality is the number of
    matching amino acid positions over the total number of amino acids.
    The sequence is the amino acid sequence of the domain on the
    reference/original transcript

### fusion\_sequence\_fasta\_id

The sequence identifier for the cdna sequence output fasta file

### fusion\_sequence\_fasta\_file

`FILEPATH` - Path to the corresponding fasta output file

### annotation\_figure

`FILEPATH` - File path to the svg drawing representing the
    annotation

### annotation\_figure\_legend

`JSON` - [JSON](./../glossary#JSON) data for the figure
    legend

### genes\_encompassed

Applies to intrachromosomal events only. List of genes which overlap
    any region that occurs between both breakpoints. For example in a
    deletion event these would be deleted genes.

### genes\_overlapping\_break1

list of genes which overlap the first breakpoint

### genes\_overlapping\_break2

list of genes which overlap the second breakpoint

### genes\_proximal\_to\_break1

list of genes near the breakpoint and the distance away from the
    breakpoint

### genes\_proximal\_to\_break2

list of genes near the breakpoint and the distance away from the
    breakpoint

### break1\_chromosome

`str` - The name of the chromosome
    on which breakpoint 1 is situated

### break1\_position\_start

`int` - Start integer inclusive
    1-based of the range representing breakpoint 1

### break1\_position\_end

`int` - End integer inclusive
    1-based of the range representing breakpoint 1

### break1\_orientation

`~mavis.constants.ORIENT` - The side
    of the breakpoint wrt the positive/forward strand that is retained.

### break1\_strand

`~mavis.constants.STRAND` - The
    strand wrt to the reference positive/forward strand at this
    breakpoint.

### break1\_seq

`str` - The sequence up to and
    including the breakpoint. Always given wrt to the positive/forward
    strand

### break2\_chromosome

The name of the chromosome on which breakpoint 2 is situated

### break2\_position\_start

`int` - Start integer inclusive
    1-based of the range representing breakpoint 2

### break2\_position\_end

`int` - End integer inclusive
    1-based of the range representing breakpoint 2

### break2\_orientation

`~mavis.constants.ORIENT` - The side
    of the breakpoint wrt the positive/forward strand that is retained.

### break2\_strand

`~mavis.constants.STRAND` - The
    strand wrt to the reference positive/forward strand at this
    breakpoint.

### break2\_seq

`str` - The sequence up to and
    including the breakpoint. Always given wrt to the positive/forward
    strand

### opposing\_strands

`bool` - Specifies if breakpoints
    are on opposite strands wrt to the reference. Expects a boolean

### stranded

`bool` - Specifies if the sequencing
    protocol was strand specific or not. Expects a boolean

### protocol

`~mavis.constants.PROTOCOL` -
    Specifies the type of library

### tools

The tools that called the event originally from the cluster step.
    Should be a semi-colon delimited list of \<tool name\>\_\<tool
    version\>

### contigs\_assembled

`int` - Number of contigs that were
    built from split read sequences

### contigs\_aligned

`int` - Number of contigs that were
    able to align

### contig\_alignment\_query\_name

The query name for the contig alignment. Should match the \'read\'
    name(s) in the .contigs.bam output file

### contig\_seq

`str` - Sequence of the current
    contig wrt to the positive forward strand if not strand specific

### contig\_remap\_score

`float` - Score representing the
    number of sequences from the set of sequences given to the assembly
    algorithm that were aligned to the resulting contig with an
    acceptable scoring based on user-set thresholds. For any sequence
    its contribution to the score is divided by the number of mappings
    to give less weight to multimaps

### call\_sequence\_complexity

`float` - The minimum amount any two
    bases account for of the proportion of call sequence. An average for
    non-contig calls

### contig\_remapped\_reads

`int` - the number of reads from the
    input bam that map to the assembled contig

### contig\_remapped\_read\_names

read query names for the reads that were remapped. A -1 or -2 has
    been appended to the end of the name to indicate if this is the
    first or second read in the pair

### contig\_alignment\_score

`float` - A rank based on the
    alignment tool blat etc. of the alignment being used. An average if
    split alignments were used. Lower numbers indicate a better
    alignment. If it was the best alignment possible then this would be
    zero.

### contig\_alignment\_reference\_start

The reference start(s) \<chr\>:\<position\> of the contig alignment.
    Semi-colon delimited

### contig\_alignment\_cigar

The cigar string(s) representing the contig alignment. Semi-colon
    delimited

### contig\_remap\_coverage

`float` - Fraction of the contig
    sequence which is covered by the remapped reads

### contig\_build\_score

`int` - Score representing the edge
    weights of all edges used in building the sequence

### contig\_strand\_specific

`bool` - A flag to indicate if it
    was possible to resolve the strand for this contig

### spanning\_reads

`int` - the number of spanning reads
    which support the event

### spanning\_read\_names

read query names of the spanning reads which support the current
    event

### call\_method

`~mavis.constants.CALL_METHOD` - The
    method used to call the breakpoints

### flanking\_pairs

`int` - Number of read-pairs where
    one read aligns to the first breakpoint window and the second read
    aligns to the other. The count here is based on the number of unique
    query names

### flanking\_pairs\_compatible

`int` - Number of flanking pairs of
    a compatible orientation type. This applies to insertions and
    duplications. Flanking pairs supporting an insertion will be
    compatible to a duplication and flanking pairs supporting a
    duplication will be compatible to an insertion (possibly indicating
    an internal translocation)

### flanking\_median\_fragment\_size

`int` - The median fragment size of
    the flanking reads being used as evidence

### flanking\_stdev\_fragment\_size

`float` - The standard deviation in
    fragment size of the flanking reads being used as evidence

### break1\_split\_reads

`int` - Number of split reads that
    call the exact breakpoint given

### break1\_split\_reads\_forced

`int` - Number of split reads which
    were aligned to the opposite breakpoint window using a targeted
    alignment

### break2\_split\_reads

`int` - Number of split reads that
    call the exact breakpoint given

### break2\_split\_reads\_forced

`int` - Number of split reads which
    were aligned to the opposite breakpoint window using a targeted
    alignment

### linking\_split\_reads

`int` - Number of split reads that
    align to both breakpoints

### untemplated\_seq

`str` - The untemplated/novel
    sequence between the breakpoints

### break1\_homologous\_seq

`str` - Sequence in common at the
    first breakpoint and other side of the second breakpoint

### break2\_homologous\_seq

`str` - Sequence in common at the
    second breakpoint and other side of the first breakpoint

### break1\_ewindow

`int-int` - Window where evidence was gathered for the first
    breakpoint

### break1\_ewindow\_count

`int` - Number of reads
    processed/looked-at in the first evidence window

### break1\_ewindow\_practical\_coverage

`float` -
    break2\_ewindow\_practical\_coverage, break1\_ewindow\_count /
    len(break1\_ewindow). Not the actual coverage as bins are sampled
    within and there is a read limit cutoff

### break2\_ewindow

`int-int` - Window where evidence was gathered for the second
    breakpoint

### break2\_ewindow\_count

`int` - Number of reads
    processed/looked-at in the second evidence window

### break2\_ewindow\_practical\_coverage

`float` -
    break2\_ewindow\_practical\_coverage, break2\_ewindow\_count /
    len(break2\_ewindow). Not the actual coverage as bins are sampled
    within and there is a read limit cutoff

### raw\_flanking\_pairs

`int` - Number of flanking reads
    before calling the breakpoint. The count here is based on the number
    of unique query names

### raw\_spanning\_reads

`int` - Number of spanning reads
    collected during evidence collection before calling the breakpoint

### raw\_break1\_split\_reads

`int` - Number of split reads before
    calling the breakpoint

### raw\_break2\_split\_reads

`int` - Number of split reads before
    calling the breakpoint

### cdna\_synon

semi-colon delimited list of transcript ids which have an identical
    cdna sequence to the cdna sequence of the current fusion product

### protein\_synon

semi-colon delimited list of transcript ids which produce a
    translation with an identical amino-acid sequence to the current
    fusion product

### tracking\_id

column used to store input identifiers from the original SV calls.
    Used to track calls from the input files to the final outputs.

### fusion\_protein\_hgvs

`str` - Describes the fusion protein
    in HGVS notation. Will be None if the change is not an indel or is
    synonymous

### net\_size

`int-int` - The net size of an event. For translocations and
    inversion this will always be 0. For indels it will be negative for
    deletions and positive for insertions. It is a range to accommodate
    non-specific events.

### supplementary\_call

`bool` - Flag to indicate if the
    current event was a supplementary call, meaning a call that was
    found as a result of validating another event.
