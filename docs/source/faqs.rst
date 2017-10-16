Frequently Asked Questions
=============================

How can I use an unsupported tool?
-----------------------------------

MAVIS supports a lot of SVCallers natively meaning that it can read their output files directly using built-in
conversion options. However tools are evolving and being created constantly. To allow the user to stay up-to-date
with the latest tools MAVIS defines a :ref:`standard input file <mavis-input-format>`.
Using this standard input format means that users who wish to use an unsupported input filetype do not need to edit
the MAVIS codebase to do so.

**What to do when the target tool doesn't output all the necessary information?**

Don't worry, this is the case for a lot of tools. MAVIS accepts unknowns for this reason for some of the required
columns. These unknown or not-specified values are then expanded to all possible combinations during the clustering
step. For this reason it is sometimes helpful to use the :term:`tracking_id` column to track your calls through the MAVIS
pipeline.


How can I track my SV calls?
------------------------------

There are a lot of steps to the MAVIS pipeline and calls may be collapsed or expanded throughout the process. To ensure
you can trace your original calls through to the final output, MAVIS uses UUID identifiers assigned at the clustering stage
which are mapped to your original calls through the assignment mapping file which is output during clustering.


How to build the reference annotations input file?
-----------------------------------------------------

Instructions on generating annotations from ensembl can be found on the :ref:`reference annotations page <generate-reference-annotations>`.
It is also possible to use other sources/databases for the annotations but is left up to the user to convert them to the expected
format. See :ref:`an example here <reference-files-annotations>`.
