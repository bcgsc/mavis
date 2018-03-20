
About
======

.. mdinclude:: ../../README.md

.. _non-python-dependencies:

Aligner dependency
-----------------------

In addition to the python package dependencies, MAVIS also requires an aligner to be installed. Currently the only
aligners supported are :term:`blat` and :term:`bwa mem <bwa>`. For MAVIS to run successfully the aligner must be installed and accessible on the 
path. If you have a non-std install you may find it useful to edit the PATH environment variable. For example

.. code:: bash
    
    export PATH=/path/to/directory/containing/blat/binary:$PATH

:term:`Blat <blat>` is the default aligner. To configure MAVIS to use :term:`bwa mem <bwa>` as a default instead, use the
:ref:`MAVIS environment variables <config-environment>`. Both the :term:`aligner` and :ref:`aligner reference <reference-files-aligner-reference>` settings
should be specified

.. code:: bash

    export MAVIS_ALIGNER='bwa mem'
    export MAVIS_ALIGNER_REFERENCE=/path/to/mem/fasta/ref/file

.. note:: 

    Although MAVIS does attempt to standardize alignments there will still be some difference in the coordinates of the final
    call set dependent on the aligner used to align putatative contigs. Additionally the aligner used on the input bam
    will have a more significant impact as it will affect the reads collected in addition to the coordintates of all non-contig
    calls.


.. _resource-requirements:

Resource Requirements
-----------------------

MAVIS has been tested on both unix and linux systems. For the standard pipeline, the validation stage is
the most computationally expensive. This will vary depending on the size of your input bam file and
the number of events input to be validated. There are a number of settings that can be adjusted to reduce
memory and cpu requirements depending on what the user is trying to analyze.  

Uninformative Filter
++++++++++++++++++++++

For example, if the user is only interested in events in genes, then the :term:`uninformative_filter` can be used. 
This will drop all events that are not within a certain distance (:term:`max_proximity`) to any annotation in the 
annotations reference file. These events will be dropped prior to the validation stage which results in 
significant speed up.

This can be set using the environment variable

.. code::

    export MAVIS_UNINFORMATIVE_FILTER=True

or in the pipeline config file

.. code::

    [cluster]
    uninformative_filter = True

or as a command line argument to the cluster stage

.. code::

    mavis cluster --uninformative_filter True ....

Splitting Validation into Cluster Jobs
+++++++++++++++++++++++++++++++++++++++

MAVIS chooses the number of jobs to split validate/annotate stages into based on
two settings: :term:`max_files` and :term:`min_clusters_per_file`.

For example, in the following situation say you have: 1000 clusters, ``max_files=10``, and ``min_clusters_per_file=10``. Then
MAVIS will set up 10 validation jobs each with 100 events.

However, if ``min_clusters_per_file=500``, then MAVIS would only set up 2 jobs each with 500 events. This is because
:term:`min_clusters_per_file` takes precedence over :term:`max_files`. 

