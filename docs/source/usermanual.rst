
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



Input Files
....................

The requirements are described in `JIRA <https://www.bcgsc.ca/jira/browse/APA-618>`_ and are listed below.
These pertain to the input files from the various tools you want to merge. The expected input columns are given
below. All columns must be given except: individual breakpoint strand columns do not need to be given if
the input is not stranded and opposing_strands has been specified

Required Columns
,,,,,,,,,,,,,,,,,

- :term:`break1_chromosome`
- :term:`break1_position_start`
- :term:`break1_position_end`
- :term:`break1_strand`
- :term:`break1_orientation`
- :term:`break2_chromosome`
- :term:`break2_position_start`
- :term:`break2_position_end`
- :term:`break2_strand`
- :term:`break2_orientation`
- :term:`opposing_strands`
- :term:`stranded`
- :term:`library`
- :term:`protocol`
- :term:`tools`


Conversion scripts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The output of the calls to be merged from the various tools must first be put into a standard/common format. For the tools we commonly use
the locations to conversion scripts are listed below

+----------------------------------------------------------------------------+------------------------------------------------------------------------+
| tool name                                                                  | formatting script                                                      |
+============================================================================+========================================================================+
| `trans-abyss <http://www.bcgsc.ca/platform/bioinfo/software/trans-abyss>`_ | https://svn.bcgsc.ca/svn/SVIA/svmerge/trunk/tools/convert_ta.py        |
+----------------------------------------------------------------------------+------------------------------------------------------------------------+
| `DELLY <https://github.com/dellytools/delly>`_                             | https://svn.bcgsc.ca/svn/SVIA/delly/trunk/delly_vcf_2_tsv.py           |
+----------------------------------------------------------------------------+------------------------------------------------------------------------+
| `Manta <https://github.com/Illumina/manta>`_                               | https://svn.bcgsc.ca/svn/SVIA/manta/trunk/manta_svmerge.py             |
+----------------------------------------------------------------------------+------------------------------------------------------------------------+
| `deFUSE <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3098195/>`_          | https://svn.bcgsc.ca/svn/SVIA/deFUSE_scripts/trunk/deFUSE.svmerge.py   |
+----------------------------------------------------------------------------+------------------------------------------------------------------------+
| `chimerascan <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3187648/>`_     | https://svn.bcgsc.ca/svn/SVIA/chimerascan/trunk/chimerascan_svmerge.py |
+----------------------------------------------------------------------------+------------------------------------------------------------------------+



|

------

|



Reference Files
..................

There are several reference files that are required for full functionality of the |TOOLNAME| pipeline. If the same
reference file will be reused often then the user may find it helpful to set reasonable defaults. Default values
for any of the reference file arguments can be configured through ``MAVIS_`` prefixed environment variables.

+--------------------------------------------------------------+-----------------------------+-----------------------------+
| file                                                         | file type/format            | environment variable        |
+==============================================================+=============================+=============================+
| :ref:`reference genome <reference-files-reference-genome>`   | :term:`fasta`               | ``MAVIS_REFERENCE_GENOME``  |
+--------------------------------------------------------------+-----------------------------+-----------------------------+
| :ref:`annotations <reference-files-annotations>`             | :term:`JSON` or text/tabbed | ``MAVIS_ANNOTATIONS``       |
+--------------------------------------------------------------+-----------------------------+-----------------------------+
| :ref:`masking <reference-files-masking>`                     | text/tabbed                 | ``MAVIS_MASKING``           |
+--------------------------------------------------------------+-----------------------------+-----------------------------+
| :ref:`template metadata <reference-files-template-metadata>` | text/tabbed                 | ``MAVIS_TEMPLATE_METADATA`` |
+--------------------------------------------------------------+-----------------------------+-----------------------------+


If the environment variables above are set they will be used as the default values when any step of the pipeline
script is called (including generating the template config file)


.. _reference-files-reference-genome:

Reference Genome
,,,,,,,,,,,,,,,,,,,,,,,

These are the sequence files in fasta format that are used in aligning and generating the fusion sequences.

**Examples:**

- `UCSC hg19 chromosome fasta sequences <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/>`_

.. _reference-files-annotations:

Annotations
,,,,,,,,,,,,,,,,,,,,,,,

This is a custom file format. Essentially just a tabbed or :term:`JSON` file which contains the gene, transcript, exon,
translation and protein domain positional information

.. warning::

    the :func:`~mavis.annotate.file_io.load_reference_genes` will
    only load valid translations. If the cds sequence in the annotation is not
    a multiple of :attr:`~mavis.constants.CODON_SIZE` or if a
    reference genome (sequences) is given and the cds start and end are not
    M and * amino acids as expected the translation is not loaded

Example of the :term:`JSON` file structure can be seen below

.. code-block:: javascript

    [
        {
            "name": string,
            "start": int,
            "end": int
            "aliases": [string, string, ...],
            "transcripts": [
                {
                    "name": string,
                    "start": int,
                    "end": int,
                    "exons": [
                        {"start": int, "end": int, "name": string},
                        ...
                    ],
                    "cdna_coding_start": int,
                    "cdna_coding_end": int,
                    "domains": [
                        {
                            "name": string,
                            "regions": [
                                {"start" aa_start, "end": aa_end}
                            ],
                            "desc": string
                        },
                        ...
                    ]
                },
                ...
            ]
        },
        ...
    }

This reference file can be generated from any database with the necessary information.
There is a `basic perl script <https://svn.bcgsc.ca/svn/SVIA/svmerge/tools/generate_ensembl_json.pl>`_
to generate the :term:`JSON` file using a connection to the `Ensembl <http://uswest.ensembl.org/index.html>`_ perl api.

.. _reference-files-template-metadata:

Template Metadata
,,,,,,,,,,,,,,,,,,,,,,,,

This is the file which contains the band information for the chromosomes. This is only used during visualization.

**Examples:**

- `UCSC hg19 cytoband file <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>`_.

.. code-block:: text

    chr1    0       2300000 p36.33  gneg
    chr1    2300000 5400000 p36.32  gpos25
    chr1    5400000 7200000 p36.31  gneg
    chr1    7200000 9200000 p36.23  gpos25
    chr1    9200000 12700000        p36.22  gneg

.. _reference-files-masking:

Masking File
,,,,,,,,,,,,,,,,,,,,,,,

File which contains regions that we should ignore calls in. This can be used to filter out
regions with known false positives, bad mapping, centromeres, telomeres etc. An example is
shown below

.. code-block:: text

    #chr    start   end     name
    chr1    0       2300000 centromere
    chr1    9200000 12700000        telomere

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


.. |TOOLNAME| replace:: **MAVIS**
