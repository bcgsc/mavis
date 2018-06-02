.. _guidelines-for-contributors:

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

Computing Code coverage
-------------------------

Since MAVIS uses multiple processes, it adds complexity to computing the code coverage. Running coverage normally will undereport.
To ensure that the coverage module captures the information from the subprocesses we need to do the following

In our development python virtual environment put a coverage.pth file (ex. ``venv/lib/python3.6/site-packages/coverage.pth``) containing the following

.. code:: python

    import coverage; coverage.process_startup()

Additionally you will need to set the environment variable

.. code:: bash

    export COVERAGE_PROCESS_START=/path/to/mavis/repo/mavis/.coveragerc
