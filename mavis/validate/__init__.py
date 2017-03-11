"""
Sub-package Documentation
==========================

The validation sub-package is responsible for pulling supporting reads from the bam file 
and re-calling events based on the evidence in a standard notation.

Types of Output Files
----------------------

A variety of intermediate output files are given for the user. These can be used to "drill down"
further into events and also for developers can be useful for debugging when adding new 
features, etc.

+----------------------------+------------------------+------------------------------------+
| expected name/suffix       | file type/format       | content                            |
+============================+========================+====================================+
| ``.raw_evidence.bam``      | :term:`bam file`       | raw evidence                       |
+----------------------------+------------------------+------------------------------------+
| ``.contigs.bam``           | :term:`bam file`       | aligned contigs                    |
+----------------------------+------------------------+------------------------------------+
| ``.evidence.bed``          | :term:`bed file`       | evidence collection window regions |
+----------------------------+------------------------+------------------------------------+
| ``.validation-passed.bed`` | :term:`bed file`       | validated event positions          |
+----------------------------+------------------------+------------------------------------+
| ``.validation-failed.tab`` | text/tabbed            | failed events                      |
+----------------------------+------------------------+------------------------------------+
| ``.validation-passed.tab`` | text/tabbed            | validated events                   |
+----------------------------+------------------------+------------------------------------+
| ``.contigs.fa``            | :term:`fasta file`     | assembled contigs                  |
+----------------------------+------------------------+------------------------------------+
| ``.contigs.blat_out.pslx`` | :term:`pslx file`      | results from blatting contigs      |
+----------------------------+------------------------+------------------------------------+
| ``.igv.batch``             | :term:`IGV batch file` | igv batch file                     |
+----------------------------+------------------------+------------------------------------+


Algorithm Overview
--------------------

- (For each breakpoint pair)

    - :ref:`Calculate the window/region <theory-calculating-the-evidence-window>` to read from the bam and collect 
      evidence
    - Store evidence (:term:`flanking pairs`, :term:`half-mapped reads`, :term:`spanning reads`, :term:`split reads`, 
      :term:`compatible flanking pairs`) which match the expected event type and position
    - :ref:`Assemble a contig <theory-assembling-contigs>` from the collected reads

- Generate a Fasta File containing all the contig sequences
- Align contigs to the reference genome (currently blat is used to perform this step)
- Make the final event calls
- (For each breakpoint pair)

    - call by contig
    - if fails, then call by spanning reads
    - if fails, then call by split reads
    - if fails, then call by mixed split/flanking reads
    - if fails, then :ref:`call by flanking pairs <theory-calling-breakpoints-by-flanking-evidence>`
    - if fails, then the event is failed

- (For each breakpoint pair)

    - :ref:`determine the amount of support <theory-determining-flanking-support>` for the more specific call

- Output new calls, evidence, contigs, etc

"""
from .main import main
