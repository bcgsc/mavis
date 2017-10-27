
.. figure:: _static/acronym.svg

About
---------

MAVIS is a pipeline to merge and validate input from different structural variant callers into a single report.
The pipeline consists four of main steps

- :ref:`cluster <mavis-cluster>`
- :ref:`validate <mavis-validate>`
- :ref:`annotate <mavis-annotate>`
- :ref:`pairing <mavis-pairing>`


Getting started
...................

There are 3 major steps to setting up and installing MAVIS

1. **Install non-python dependencies**

Before MAVIS can be installed, the :ref:`non-python dependencies <non-python-dependencies>` will need to be installed.
After these have been installed MAVIS itself can be installed through pip

2. **Install MAVIS**

.. code:: bash

    pip install mavis

This will install mavis and its python dependencies.

3. **Build reference files**

After MAVIS is installed the :ref:`reference files <reference-input>` must be generated before it can be run.

Once the above 3 steps are complete MAVIS is ready to be run. See :ref:`running the pipeline <pipeline>`.


.. _non-python-dependencies:

Non-python dependencies
.........................

Aligner (:term:`blat`)
+++++++++++++++++++++++++

In addition to the python package dependencies, MAVIS also requires an aligner to be installed. Currently the only
aligners supported are :term:`blat` and :term:`bwa mem`. For MAVIS to run successfully the aligner must be installed and accessible on the 
path. If you have a non-std install you may find it useful to edit the PATH environment variable. For example

.. code:: bash
    
    export PATH=/path/to/directory/containing/blat/binary:$PATH

:term:`Blat <blat>` is the default aligner. To configure MAVIS to use :term:`bwa mem` as a default instead, use the
:ref:`MAVIS environment variables <config-environment>`. Both the :term:`aligner` and :ref:`aligner reference <reference-files-aligner-reference>` settings
should be specified

.. code:: bash

    export MAVIS_ALIGNER='bwa mem'
    export MAVIS_ALIGNER_REFERENCE=/path/to/mem/fasta/ref/file


Samtools
++++++++++++++++++

Samtools is only used in sorting and indexing the intermediary output bams. Eventually this will hopefully be 
accomplished through :term:`pysam` only.
