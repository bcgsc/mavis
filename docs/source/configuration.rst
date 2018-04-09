
Configuration and Settings
=============================

.. _pipeline-config:

Pipeline Configuration File
-------------------------------

The pipeline can be run in steps or it can be configured using a configuration file and setup in a single step. Scripts
will be generated to run all steps following clustering. The configuration file can be built from scratch or a template
can be output as shown below

.. code-block:: bash

    >>> mavis config --write template.cfg

This will create a template config file called template.cfg which can then be edited by the user. However this will be
a simple config with no library information. To generate a configuration file with the library information as well as
estimates for the fragment size parameters more inputs are required (see :ref:`generating the config file <example-generating-the-conf>` for more information).


.. _config-environment:


Environment Variables
---------------------------

Most of the default settings can be changed by using environment variables. The value given by the
environment variables will be used as the new default. Config or command-line parameters will still
override these settings.

All environment variables are prefixed with MAVIS and an underscore. Otherwise the variable name is the same
as that used for the command line parameter or config setting (uppercased). For example to change the default minimum mapping
quality used during the validate stage

.. code-block:: bash

    >>> export MAVIS_MIN_MAPPING_QUALITY=10


.. _resource-requirements:

Resource Requirements
----------------------------------

MAVIS has been tested on both unix and linux systems. For the standard pipeline, the validation stage is
the most computationally expensive. This will vary depending on the size of your input bam file and
the number of events input to be validated. There are a number of settings that can be adjusted to reduce
memory and cpu requirements depending on what the user is trying to analyze.  

Uninformative Filter
......................

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
...........................................

MAVIS chooses the number of jobs to split validate/annotate stages into based on
two settings: :term:`max_files` and :term:`min_clusters_per_file`.

For example, in the following situation say you have: 1000 clusters, ``max_files=10``, and ``min_clusters_per_file=10``. Then
MAVIS will set up 10 validation jobs each with 100 events.

However, if ``min_clusters_per_file=500``, then MAVIS would only set up 2 jobs each with 500 events. This is because
:term:`min_clusters_per_file` takes precedence over :term:`max_files`. 


.. |TOOLNAME| replace:: **MAVIS**
