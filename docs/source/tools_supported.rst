.. _supported-sv-callers:


Supported SV Callers
================================

MAVIS supports output from a wide-variety of SV callers. Assumptions are made for each tool based on interpretation of
the output and the publications for each tool. The tools and versions currently supported are given below. Versions listed
indicate the version of the tool for which output files have been tested as input into MAVIS

.. list-table::
    :header-rows: 1

    *   - Name
        - Version(s)
        - MAVIS input
        - Publication
    *   - `BreakDancer <https://github.com/genome/breakdancer>`_
        - ``1.4.5``
        -
        - [Chen-2009]_
    *   - `Chimerascan <https://code.google.com/archive/p/chimerascan>`_
        - ``0.4.5``
        - ``*.bedpe``
        - [Iyer-2011]_
    *   - `DeFuse <https://bitbucket.org/dranew/defuse>`_
        - ``0.6.2``
        - ``results/results.classify.tsv``
        - [McPherson-2011]_
    *   - `DELLY <https://github.com/dellytools/delly>`_
        - ``0.6.1`` ``0.7.3``
        - ``combined.vcf`` (converted from bcf)
        - [Rausch-2012]_
    *   - `Manta <https://github.com/Illumina/manta>`_
        - ``1.0.0``
        - ``{diploidSV,somaticSV}.vcf``
        - [Chen-2016]_
    *   - `Pindel <https://github.com/genome/pindel>`_
        - ``0.2.5b9``
        -
        - [Ye-2009]_
    *   - `Trans-ABySS <https://github.com/bcgsc/transabyss>`_
        - ``1.4.8 (custom)``
        - ``fusions/*.tsv``
        - [Robertson-2010]_

.. note::

    The trans-abyss version used was an in-house dev version. However the output columns are compatible with 1.4.8 as that
    was the version branched from


:ref:`DELLY <Rausch-2012>` Post-processing
---------------------------------------------

Some post-processing on the delly output files is generally done prior to input. The output bcf files are converted to a vcf file

.. code:: bash

    bcftools concat -f /path/to/file/with/vcf/list --allow-overlaps --output-type v --output combined.vcf


.. _custom-conversion:

Writing A Custom Conversion Script
-----------------------------------

Logic Example - :ref:`Chimerascan <Iyer-2011>`
++++++++++++++++++++++++++++++++++++++++++++++++++++++


The following is a description of how the conversion script for :ref:`Chimerascan <Iyer-2011>` was generated. While this is a built-in
conversion command now, the logic could also have been put in an external script. As mentioned above, there are a number of
assumptions that had to be made about the tools output to convert it to the :ref:`standard mavis format <mavis-input-format>`.
Assumptions were then verified by reviewing at a series of called events in :term:`IGV`. In the current example,
:ref:`Chimerascan <Iyer-2011>` output has six columns of interest that were used in the conversion

- start3p
- end3p
- strand3p
- start5p
- end5p
- strand5p

The above columns describe two segments which are joined. MAVIS requires the position of the join. It was assumed that the
segments are always joined as a :term:`sense fusion`. Using this assumption there are four logical cases to determine the position of the breakpoints.

i.e. the first case would be: If both strands are positive, then the end of the five-prime segment (end5p)
is the first breakpoint and the start of the three-prime segment is the second breakpoint

The logic for all cases is shown in the code below


.. literalinclude:: ./../../mavis/tools.py
    :pyobject: _parse_chimerascan
    :emphasize-lines: 10-22

Calling A Custom Conversion Script
+++++++++++++++++++++++++++++++++++++

Custom conversion scripts can be specified during :ref:`automatic config generation <pipeline-config>` using the
``--external_conversion`` option.

.. note::

    Any external conversion scripts must take a ``-o`` option which requires a single
    outputfile argument. This outputfile must be the converted file output by the script.
    Additionally, the conversion script must be specified by its full path name and have executeable permissions.

In the following example the user has created a custom conversion script ``my_convert_script.py`` which they
are passing an input file named ``my_input1.txt``.

.. code:: bash

    mavis config --external_conversion my_converted_input1 "my_convert_script.py my_input1.txt ... "

This will then be called during the pipeline step as

.. code:: bash

    my_convert_script.py my_input1.txt ... -o /path/to/output/dir/converted_inputs/my_converted_input1.tab


You can also re-use the same conversion script if you have multiple inputs to convert simply by specifying a different alias

.. code:: bash

    mavis config \
        --external_conversion my_converted_input1 "my_convert_script.py my_input1.txt" \
        --external_conversion my_converted_input2 "my_convert_script.py my_input2.txt"
