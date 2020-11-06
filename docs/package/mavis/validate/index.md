# Sub-package Documentation

The validation sub-package is responsible for pulling supporting reads from the bam file
and re-calling events based on the evidence in a standard notation.

## Types of Output Files

A variety of intermediate output files are given for the user. These can be used to "drill down"
further into events and also for developers debugging when adding new features, etc.

| expected name/suffix        | file type/format                                    | content                            |
| --------------------------- | --------------------------------------------------- | ---------------------------------- |
| ``*.raw_evidence.bam``      | [bam](../../../glossary/#bam)                       | raw evidence                       |
| ``*.contigs.bam``           | [bam](../../../glossary/#bam)                       | aligned contigs                    |
| ``*.evidence.bed``          | [bed](../../../glossary/#bed)                       | evidence collection window regions |
| ``*.validation-passed.bed`` | [bed](../../../glossary/#bed)                       | validated event positions          |
| ``*.validation-failed.tab`` | text/tabbed                                         | failed events                      |
| ``*.validation-passed.tab`` | text/tabbed                                         | validated events                   |
| ``*.contigs.fa``            | [fasta](../../../glossary/#fasta)                   | assembled contigs                  |
| ``*.contigs.blat_out.pslx`` | [pslx](../../../glossary/#pslx)                     | results from blatting contigs      |
| ``*.igv.batch``             | [IGV batch file](../../../glossary/#IGV-batch-file) | igv batch file                     |


## Algorithm Overview

- (For each breakpoint pair)

    - [Calculate the window/region](../../../background/theory/#calculating-the-evidence-window) to read from the bam and collect
      evidence
    - Store evidence ([flanking read pair](../../../glossary/#flanking-read-pair), [half-mapped read](../../../glossary/#half-mapped-read), [spanning read](../../../glossary/#spanning-read), [split read](../../../glossary/#split-read),
      [compatible flanking pairs](../../../glossary/#compatible-flanking-pairs)) which match the expected event type and position
    - Assemble a contig from the collected reads. see [theory - assembling contigs](../../../background/theory/#assembling-contigs)

- Generate a [fasta](../../../glossary/#fasta) file containing all the contig sequences
- Align contigs to the reference genome (currently [blat](../../../glossary/#blat) is used to perform this step)
- Make the final event calls. Each level of calls consumes all supporting reads so they are not re-used in subsequent
  levels of calls.
- (For each breakpoint pair)

    - call by contig
    - call by [spanning read](../../../glossary/#spanning-read)
    - call by [split read](../../../glossary/#split-read)
    - call by [flanking read pair](../../../glossary/#flanking-read-pair). see [theory - calling breakpoints by flanking evidence](../../../background/theory/#calling-breakpoints-by-flanking-evidence)

- Output new calls, evidence, contigs, etc
