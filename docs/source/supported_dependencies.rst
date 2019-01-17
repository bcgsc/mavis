Supported Dependencies
=======================

MAVIS integrates with :ref:`SV callers <supported-sv-callers>`, :ref:`job schedulers <dep-job-schedulers>`, and :ref:`aligners <dep-aligners>`. While some of these dependencies are optional, all currently supported options are detailed below. The versions column
in the tables below list all the versions which were tested for each tool. Each version listed is known to be compatible with MAVIS.

.. _dep-job-schedulers:

Job Schedulers
--------------------

MAVIS can be run locally without a job scheduler (``MAVIS_SCHEDULER=LOCAL``) however, due to the computational resources generally required, it is recommended that you
use one of the supported schedulers listed below.


.. list-table::
    :header-rows: 1

    *   - Name
        - Version(s)
        - Environment Setting
    *   - :term:`TORQUE`
        - ``6.1.2``
        - ``MAVIS_SCHEDULER=TORQUE``
    *   - :term:`SGE`
        - ``8.1.8``
        - ``MAVIS_SCHEDULER=SGE``
    *   - :term:`SLURM`
        - ``17.02.1-2``
        - ``MAVIS_SCHEDULER=SLURM``

Users requiring support for other schedulers may make a request by `submitting an issue to our github page <https://github.com/bcgsc/mavis/issues>`__. Additionally, developers looking to
extend the functionality may submit a pull request (Please see the :ref:`guidelines for contributors <guideline-for-contributors>`).

MAVIS running locally uses the python :py:mod:`concurrent.futures` library to manage jobs.


....


.. _dep-aligners:

Aligners
-----------

Two aligners are supported :term:`bwa` and :term:`blat` (default).

.. list-table::
    :header-rows: 1

    *   - Name
        - Version(s)
        - Environment Setting
    *   - :term:`blat`
        - ``36x2`` ``36``
        - ``MAVIS_ALIGNER=blat``
    *   - :term:`bwa mem <bwa>`
        - ``0.7.15-r1140`` ``0.7.12``
        - ``MAVIS_ALIGNER='bwa mem'``

.. note::

    When setting the aligner you will also need to set the :term:`aligner_reference` to match


....


.. _supported-sv-callers:

SV Callers
-----------------------

MAVIS supports output from a wide-variety of :term:`SV` callers. Assumptions are made for each tool based on interpretation of
the output and the publications for each tool. The tools and versions currently supported are given below. Versions listed
indicate the version of the tool for which output files have been tested as input into MAVIS

MAVIS also supports a :ref:`general VCF input <general-vcf-inputs>`. It should be noted however that the tool tracked will only be listed as 'vcf' then.

.. list-table::
    :header-rows: 1

    *   - Name
        - Version(s)
        - MAVIS input
        - Publication
    *   - :term:`BreakDancer`
        - ``1.4.5``
        - ``Tools main output file(s)``
        - [Chen-2009]_
    *   - :term:`BreakSeq`
        - ``2.2``
        - ``work/breakseq.vcf.gz``
        - [Abyzov-2015]_
    *   - :term:`Chimerascan`
        - ``0.4.5``
        - ``*.bedpe``
        - [Iyer-2011]_
    *   - :term:`CNVnator`
        - ``0.3.3``
        - ``Tools main output file(s)``
        - [Abyzov-2011]_
    *   - :term:`DeFuse`
        - ``0.6.2``
        - ``results/results.classify.tsv``
        - [McPherson-2011]_
    *   - :term:`DELLY`
        - ``0.6.1`` ``0.7.3``
        - ``combined.vcf`` (converted from bcf)
        - [Rausch-2012]_
    *   - :term:`Manta`
        - ``1.0.0``
        - ``{diploidSV,somaticSV}.vcf``
        - [Chen-2016]_
    *   - :term:`Pindel`
        - ``0.2.5b9``
        - ``Tools main output file(s)``
        - [Ye-2009]_
    *   - :term:`Trans-ABySS`
        - ``1.4.8 (custom)``
        - ``{indels/events_novel_exons,fusions/*}.tsv``
        - [Robertson-2010]_
    *   - :term:`Strelka`
        - ``1.0.6``
        - ``passed.somatic.indels.vcf``
        - [Saunders-2012]_
    *   - :term:`STAR-Fusion`
        - ``1.4.0``
        - ``star-fusion.fusion_predictions.abridged.tsv``
        - [Haas-2017]_
.. note::

    :term:`Trans-ABySS`: The trans-abyss version used was an in-house dev version. However the output columns are compatible with 1.4.8 as that
    was the version branched from. Additionally, although indels can be used from both genome and transcriptome outputs of Trans-ABySS, it is
    reccommended to only use the genome indel calls as the transcriptome indels calls (for versions tested) introduce a very high number of
    false positives. This will slow down validation. It is much faster to simply use the genome indels for both genome and transcriptome.


:term:`DELLY` Post-processing
..............................................................

Some post-processing on the delly output files is generally done prior to input. The output BCF files are converted to a VCF file

.. code:: bash

    bcftools concat -f /path/to/file/with/vcf/list --allow-overlaps --output-type v --output combined.vcf


.. _custom-conversion:

Writing A Custom Conversion Script
..............................................................

Logic Example - :term:`Chimerascan`
+++++++++++++++++++++++++++++++++++++


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
..............................................................

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



.. _general-vcf-inputs:


General VCF inputs
..............................................................

Assuming that the tool outputting the VCF file follows standard conventions, then it is possible to use a general VCF
conversion that is not tool-specific. Given the wide variety in content for VCF files, MAVIS makes a number of
assumptions and the VCF conversion may not work for all VCFs. In general MAVIS follows the `VCF 4.2 specification <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_.
If the input tool you are using differs, it would be better to use a :ref:`custom conversion script <custom-conversion>`.

**Assumptions on non-standard INFO fields**

- ``PRECISE`` if given, Confidence intervals are ignored if given in favour of exact breakpoint calls using pos and END as the breakpoint positions
- ``CT`` values if given are representative of the breakpoint orientations.
- ``CHR2`` is given for all interchromosomal events
