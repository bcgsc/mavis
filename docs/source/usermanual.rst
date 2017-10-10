
.. figure:: _static/acronym.svg

About
---------

|TOOLNAME| is a pipeline to merge and validate input from different structural variant callers into a single report.
The pipeline consists four of main steps

- :ref:`cluster <mavis-cluster>`
- :ref:`validate <mavis-validate>`
- :ref:`annotate <mavis-annotate>`
- :ref:`pairing <mavis-pairing>`

|

--------------

|

Getting started
--------------------

see :ref:`Installation for developers <development-install>`


|

------

|



Running the Pipeline
.....................

The pipeline can be run calling the main script (see below) followed the pipeline step. The usage menu can be viewed
by running the without any arguments, or by giving the -h/--help option

**Example:**

.. code-block:: bash

    >>> mavis


Help sub-menus can be found by giving the pipeline step followed by no arguments or the -h options

.. code-block:: bash

    >>> mavis cluster -h

The most common use case is auto-generating a configuration file and then running the pipeline setup step.
The pipeline setup step will run clustering and create scripts for running the other steps.

.. code-block:: bash

    >>> mavis config .... -w config.cfg
    >>> mavis pipeline config.cfg -o /path/to/top/output_dir

This will create submission scripts as follows

.. code-block:: text

    output_dir/
    |-- library1/
    |   |-- validation/qsub.sh
    |   `-- annotation/qsub.sh
    |-- library2/
    |   |-- validation/qsub.sh
    |   `-- annotation/qsub.sh
    |-- pairing/qsub.sh
    `-- summary/qsub.sh

The qsub scripts are bash scripts meant for submission to an `SGE <http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html>`_
cluster. The summary job is held on the pairing job, the pairing job is held on all the annotation jobs,
and the annotation jobs are held on their validation jobs. This means that the scripts should be submitted in
the following order: (validation => annotation => pairing => summary)

.. code-block:: bash

    >>> ssh cluster_head_node
    >>> qsub output_dir/library1/validation/qsub.sh
    >>> qsub output_dir/library1/annotation/qsub.sh
    >>> qsub output_dir/library2/validation/qsub.sh
    >>> qsub output_dir/library2/annotation/qsub.sh
    >>> qsub output_dir/pairing/qsub.sh
    >>> qsub output_dir/summary/qsub.sh

|

----

|


Configuration and Settings
.............................


Generating a config file automatically
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The pipeline can be run in steps or it can be configured using a configuration file and setup in a single step. Scripts
will be generated to run all steps following clustering. The configuration file can be built from scratch or a template
can be output as shown below

.. code-block:: bash

    >>> mavis config --write template.cfg

This will create a template config file called template.cfg which can then be edited by the user. However this will be
a simple config with no library information. To generate a configuration file with the library information as well as
estimates for the fragment size parameters more inputs are required.

A simple example with a single library would look like this (see below)

.. code-block:: bash

    >>> mavis config --write output.cfg \
        --library Library1 genome diseased /path/to/bam/file/library1.bam False

This creates a configuration file but is still missing some information before it can be run by the pipeline, the input
files containing the breakpoint pairs. So a more complete example is shown below

.. code-block:: bash

    >>> mavis config --write output.cfg \
        --library Library1 genome diseased /path/to/bam/file/library1.bam False \
        --library Library2 genome normal /path/to/bam/file/library2.bam False \
        --input /path/to/bpp/file Library1 Library2 \
        --input /path/to/other/bpp/file Library1 Library2

In the above example Library1 is the tumour genome and Library2 is the normal genome. The same input files are
used for both

Manually creating the configuration File
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

While not recommended, the configuration file can also be built manually. The minimum required inputs are the library
configuration sections. There must be at least one library section and the library section must at minimum have the
following attributes given (see below).

.. code-block:: python

    [Library1]
    protocol = genome
    bam_file = /path/to/bam/file/library1.bam
    read_length = 125
    median_fragment_size = 435
    stdev_fragment_size = 100
    stranded_bam = False
    inputs = /path/to/bpp/file
    disease_status = diseased


Environment Variables
,,,,,,,,,,,,,,,,,,,,,,,,,

Most of the default settings can be changed by using environment variables. The value given by the
environment variables will be used as the new default. Config or command-line parameters will still
override these settings.

All environment variables are prefixed with MAVIS and an underscore. Otherwise the variable name is the same
as that used for the command line parameter or config setting (uppercased). For example to change the default minimum mapping
quality used during the validate stage

.. code-block:: bash

    >>> export MAVIS_MIN_MAPPING_QUALITY=10





.. |TOOLNAME| replace:: **MAVIS**
