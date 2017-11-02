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

See :ref:`Writing A Custom Conversion Script <custom-conversion>` for more details 


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

Can I run mavis locally?
----------------------------

Although the default pipeline builds scripts to submit to a compute cluster, The same scripts can also be run locally. They should be executed
in the same order as the main submit script would submit them: validation, annotation, pairing, then summary. For example if your MAVIS
output looked something like this

.. code:: text

    .
    |-- lib/
    |   |-- validate/.../submit.sh
    |   `-- annotate/.../submit.sh
    |-- pairing/submit.sh
    `-- summary/submit.sh

you would run them in the following order

.. code:: bash

    bash ./lib/validate/.../submit.sh; \
        bash ./lib/annotate/.../submit.sh; \
        bash ./pairing/submit.sh; \
        bash ./summary/submit.sh

As MAVIS splits clusters into multiple jobs by default there may be many scripts to bash. If want to avoid this (an your machine has sufficient memory)
then set ``MAVIS_MAX_FILES=1`` to restrict the number of jobs/cluster files to 1.


Why Does My Event Have 0 Spanning Reads?
------------------------------------------

When looking at the evidence columns of the MAVIS output files it is important to look first at the :term:`call_method`
column. Some evidence types do not apply to certain :term:`call methods <call_method>` so they will always be 0 or None.
For example, if a small indel is called by contig, we would first look at the :term:`contig_remapped_reads` column. This
is the evidence that was used in determining whether or not to include the call in the output file.

If the :term:`break1_split_reads` column is 0 and the call method is not by :term:`split reads <split read>` it does not mean there is low evidence.


I See Split Reads in :term:`IGV`, Why Does MAVIS Call 0 Split Reads?
--------------------------------------------------------------------

If the event was called by :term:`contig` for example, the breakpoint positions will be based on the alignment of the contig.
Only split reads which exactly match this breakpoint will be given as evidence by split reads.

If the event is not an exact breakpoint call, only flanking evidence will be given


Why is the Breakpoint Called by MAVIS Different From What I See in :term:`IGV`?
---------------------------------------------------------------------------------

MAVIS normalizes the read alignments before calling events. This is especially important in repeat regions.
Aligners like :term:`bwa mem <bwa>` align deletions to the start of a repeat span, whereas MAVIS follows the
:term:`hgvs` standard and aligns deletions to the end of a repeat span.
