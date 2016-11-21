Glossary
===========

..  glossary::
    :sorted:

    reference template
        can be a chromosome or alternate assembly. For example chr1

    strand
        refers to the positive/forward or the negative/reverse strand in the template reference. Defined in the code by the enum-like object :py:attr:`~structural_variant.constants.STRAND`

    orientation
        always specified with respect to the positive/forward strand on the reference. Defined in the code by the enum-like object :py:attr:`~structural_variant.constants.ORIENT`


Read Terms
------------

.. image:: _static/svmerge_defn.svg

..  glossary::
    :sorted:

    insert size
        the distance between two paired reads

    read length
        the length of a given read

    fragment size
        the distance between (and including) a read-pair

Types of read evidence
------------------------

.. image:: _static/svmerge_defn_reads.svg

..  glossary::
    :sorted:

    split reads
        reads that span a breakpoint such that part of the read is aligned and part is soft-clipped

    spanning reads
        reads which span a set of breakpoints. Generally applicable to small indels. Spanning reads are generally a subset of split reads

    flanking reads
        refers the read pairs that flank the breakpoints

    half-mapped
        half mapped reads are reads where one partner is unmapped. Next to a breakpoint these may occur is there is a large amount of untemplated sequence inserted
