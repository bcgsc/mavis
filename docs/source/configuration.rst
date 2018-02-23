
Configuration and Settings
---------------------------

.. _pipeline-config:

Pipeline Configuration File
...............................


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
        --library Library1 genome diseased False /path/to/bam/file/library1.bam

This creates a configuration file but is still missing some information before it can be run by the pipeline, the input
files containing the breakpoint pairs. So a more complete example is shown below

.. code-block:: bash

    >>> mavis config --write output.cfg \
        --library Library1 genome diseased False /path/to/bam/file/library1.bam \
        --library Library2 genome normal False /path/to/bam/file/library2.bam \
        --assign Library1 /path/to/bpp/file /path/to/other/bpp/file \
        --assign Library2 /path/to/bpp/file /path/to/other/bpp/file

In the above example Library1 is the tumour genome and Library2 is the normal genome. The same input files are
used for both

Conversion During Configuration
++++++++++++++++++++++++++++++++++

Raw SV caller output files can be set up for conversion during the configuration stage. This will require using
the ``--convert`` option (for natively supported tools/formats). The alias given to the ``--convert`` option can
be used in the same way a path for an already converted file would be, and passed as input

For example if the input files in the above example were DELLY output files

.. code-block:: bash

    >>> mavis config --write output.cfg \
        --library Library1 genome diseased False /path/to/bam/file/library1.bam \
        --library Library2 genome normal False /path/to/bam/file/library2.bam \
        --assign Library1 alias1 other_alias \
        --assign Library2 alias1 other_alias \
        --convert alias1 /path/to/bpp/file delly False \
        --convert other_alias /path/to/other/bpp/file delly False



.. _config-environment:


Environment Variables
.......................

Most of the default settings can be changed by using environment variables. The value given by the
environment variables will be used as the new default. Config or command-line parameters will still
override these settings.

All environment variables are prefixed with MAVIS and an underscore. Otherwise the variable name is the same
as that used for the command line parameter or config setting (uppercased). For example to change the default minimum mapping
quality used during the validate stage

.. code-block:: bash

    >>> export MAVIS_MIN_MAPPING_QUALITY=10





.. |TOOLNAME| replace:: **MAVIS**
