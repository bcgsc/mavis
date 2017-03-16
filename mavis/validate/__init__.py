"""
Sub-package Documentation
==========================

The validation sub-package is responsible for pulling supporting reads from the bam file
and re-calling events based on the evidence in a standard notation.

Types of Output Files
----------------------

A variety of intermediate output files are given for the user. These can be used to "drill down"
further into events and also for developers debugging when adding new features, etc.

+-----------------------------+------------------------+------------------------------------+
| expected name/suffix        | file type/format       | content                            |
+=============================+========================+====================================+
| ``*.raw_evidence.bam``      | :term:`bam`            | raw evidence                       |
+-----------------------------+------------------------+------------------------------------+
| ``*.contigs.bam``           | :term:`bam`            | aligned contigs                    |
+-----------------------------+------------------------+------------------------------------+
| ``*.evidence.bed``          | :term:`bed`            | evidence collection window regions |
+-----------------------------+------------------------+------------------------------------+
| ``*.validation-passed.bed`` | :term:`bed`            | validated event positions          |
+-----------------------------+------------------------+------------------------------------+
| ``*.validation-failed.tab`` | text/tabbed            | failed events                      |
+-----------------------------+------------------------+------------------------------------+
| ``*.validation-passed.tab`` | text/tabbed            | validated events                   |
+-----------------------------+------------------------+------------------------------------+
| ``*.contigs.fa``            | :term:`fasta`          | assembled contigs                  |
+-----------------------------+------------------------+------------------------------------+
| ``*.contigs.blat_out.pslx`` | :term:`pslx`           | results from blatting contigs      |
+-----------------------------+------------------------+------------------------------------+
| ``*.igv.batch``             | :term:`IGV batch file` | igv batch file                     |
+-----------------------------+------------------------+------------------------------------+


Algorithm Overview
--------------------

- (For each breakpoint pair)

    - :ref:`Calculate the window/region <theory-calculating-the-evidence-window>` to read from the bam and collect
      evidence
    - Store evidence (:term:`flanking read pair`, :term:`half-mapped read`, :term:`spanning read`, :term:`split read`,
      :term:`compatible flanking pairs`) which match the expected event type and position
    - Assemble a contig from the collected reads. see :ref:`theory - assembling contigs <theory-assembling-contigs>`

- Generate a :term:`fasta` file containing all the contig sequences
- Align contigs to the reference genome (currently :term:`blat` is used to perform this step)
- Make the final event calls
- (For each breakpoint pair)

    - call by contig
    - if fails, then call by :term:`spanning read`
    - if fails, then call by :term:`split read`
    - if fails, then call by mixed :term:`split read` / :term:`flanking read pair`
    - if fails, then call by :term:`flanking read pair`. see 
      :ref:`theory - calling breakpoints by flanking evidence <theory-calling-breakpoints-by-flanking-evidence>`
    - if fails, then the event is failed

- (For each breakpoint pair)

    - determine the amount of support for the more specific call. see 
      :ref:`theory - determining flanking support <theory-determining-flanking-support>`

- Output new calls, evidence, contigs, etc

"""
