# General

## About

pipeline to merge and validate input from different structural variant callers into a single report

## Running Unit Tests

to run the tests

    cd /path/to/project/dir
    nosetests --with-coverage --cover-package=structural_variant

## Getting Started

### Input file requirements

The requirements are described in [JIRA](https://www.bcgsc.ca/jira/browse/APA-618) and are listed below.
These pertain to the input files from the various tools you want to merge


    | column name       | expected value                    | description                     |
    |-------------------|-----------------------------------|---------------------------------|
    | start_position    | <int>-<int>                       |                                 |
    | start_strand      | <+,-,?>                           | the reference strand aligned to |
    | start_orientation | <L,R,?>                           |                                 |
    | end_chromosome    | <1-22,X,Y,MT>                     |                                 |
    | end_position      | <int>-<int>                       |                                 |
    | end_strand        | <+,-,?>                           | the reference strand aligned to |
    | end_orientation   | <L,R,?>                           |                                 |
    | protocol          | <genome or transcriptome>         |                                 |
    | library           | library id                        |                                 |
    | tool_version      | <tool name>_<tool version number> |                                 |
    | opposing_strand   | <True,False,?>                    |                                 |

### Available Pre-formatting scripts

- Trans-ABySS
    - [convert_ta.py](https://svn.bcgsc.ca/svn/SVIA/sv_compile/tags/<version>/tools/convert_ta.py)
- DELLY
    - DUSTIN TO ADD
- MANTA
    - CALEB TO ADD

## Design

### general process:

#### in-progress
1. reformat input files
2. group breakpoint pairs by library and protocol
3. cluster breakpoint pairs
4. filter 'lowQual' tagged clusters with below a input threshold (i.e. singlets)
5. 'fake' reciprocal events
6. calculate the evidence window (define region to parse reads from the bam) (TODO: calculating the transcriptome window)
7. gather evidence: flanking, spanning, and split reads
8. assemble split reads into contigs
9. align split-read contigs to the genome 

#### TODO
10. split breakpoint pairs (and evidence) by classification
11. split breakpoint pairs (and evidence) by breakpoint calls (from split read contigs, or from flanking reads [approximate])
12. filter by evidence amount
13. annotate calls
14. pair between libraries based on breakpoints and annotations
15. gather further detailed documentation on resulting protein/transcript
16. draw putative fusion structure

#### proposed output structure

```
.
|-- clustering/
|   `-- <library id>/
|       `-- cluster assignment tracking files
|-- validation/
|   `-- <library id>/
|       |-- breakpoint pairs with associated evidence amounts
|       |-- bam files with reads used as evidence
|       `-- bam files with contigs from assembly
|-- annotation/
|   `-- <library id>/
|       `-- annotated breakpoint pairs
`-- pairing/
    `-- cross library putative pairings
```

### split-read evidence

#### purpose

to ensure breakpoints are called consistently and to aid in eliminating false positives

#### process

- gather split reads that fall in the respective evidence windows
- revcomp the split reads from one of the breakpoints if opposing strands have been specified
- assemble the contigs
- align the contigs to the reference genome
- return a list of pairs of alignments (non-redundant, primarily for duplications and translocations) as well as single/gapped alignments
