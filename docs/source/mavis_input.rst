MAVIS standard input file format
==================================

.. _mavis-input-format:

The requirements are described in `JIRA <https://www.bcgsc.ca/jira/browse/APA-618>`_ and are listed below.
These pertain to the input files from the various tools you want to merge. The expected input columns are given
below. All columns must be given except: individual breakpoint strand columns do not need to be given if
the input is not stranded and opposing_strands has been specified

Required Columns
,,,,,,,,,,,,,,,,,

- :term:`break1_chromosome`
- :term:`break1_position_start`
- :term:`break1_position_end` (can be the same as break1_position_start)
- :term:`break2_chromosome`
- :term:`break2_position_start`
- :term:`break2_position_end` (can be the same as break2_position_start)


Optional Columns
,,,,,,,,,,,,,,,,,

Optional Columns that are not given as input will be added with default (or command line parameter options) during
the clustering stage of MAVIS as some are required for subsequent pipeline steps

- :term:`break1_strand` (defaults to not-specified during clustering)
- :term:`break1_orientation` (expanded to all possible values during clustering)
- :term:`break2_strand` (defaults to not-specified during clustering)
- :term:`break2_orientation` (expanded to all possible values during clustering)
- :term:`opposing_strands` (expanded to all possible values during clustering)
- :term:`stranded` (defaults to False during clustering)
- :term:`library` (defaults to command line library parameter during clustering)
- :term:`protocol` (defaults to command line protocol parameter during clustering)
- :term:`tools` (defaults to an empty string during clustering)


The different pipeline steps of MAVIS have different input column requirements. These are summarized below (for the
pipeline steps which can act as the pipeline start)

+-----------------------------------+-----------+-----------+-----------+
| column name                       | cluster   | annotate  | validate  |
+===================================+===========+===========+===========+
| :term:`break1_chromosome`         | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break1_position_start`     | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break1_position_end`       | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break2_chromosome`         | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break2_position_start`     | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break2_position_end`       | X         | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break1_strand`             |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break1_orientation`        |           | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break2_strand`             |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`break2_orientation`        |           | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`opposing_strands`          |           | X         | X         |
+-----------------------------------+-----------+-----------+-----------+
| :term:`stranded`                  |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`library`                   |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`protocol`                  |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`tools`                     |           |           |           |
+-----------------------------------+-----------+-----------+-----------+
| :term:`event_type`                |           | X         |           |
+-----------------------------------+-----------+-----------+-----------+

Some native tool outputs are :ref:`supported <supported-sv-callers>` and have built in methods to convert to the above format. Any unsupported
tools can be used as long as the user converts the tools native output to match the above format.
