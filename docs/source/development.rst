Guidelines for Contributors
===================================

.. mdinclude:: ./../../.github/CONTRIBUTING.md


Major Assumptions
------------------

Some assumptions have been made when developing this project. The major ones have been listed here to
facilitate debugging/development if any of these are violated in the future.

- The input bam reads have stored the sequence wrt to the positive/forward strand and have not stored the reverse
  complement.
- The distribution of the fragment sizes in the bam file approximately follows a normal distribution.


Current Limitations
---------------------

- Assembling contigs will always fail for repeat sequences as we do not resolve this. Unlike traditional assemblies
  we cannot assume even input coverage as we are taking a select portion of the reads to assemble.
- Currently no attempt is made to group/pair single events into complex events.
- Transcriptome validation uses a collapsed model of all overlapping transcripts and is not isoform specific. Allowing
  for isoform specific validation would be computationally expensive but may be considered as an optional setting for
  future releases.
