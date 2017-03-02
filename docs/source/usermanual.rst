About
---------

|TOOLNAME| is a pipeline to merge and validate input from different structural variant callers into a single report.

- Breakpoint pair calls are read in from various tools and clustered by library/sample.
- The clustered calls are filtered by proximity to annotations (this is done as validating
  the evidence for individual calls is the most expensive part of the pipeline)
- The calls are validated against a paired-end read bam file where split, flanking, and spanning reads
  are gathered and analyzed. From each cluster new breakpoint pairs are called from the evidence.
- The pairs that passed validation are annotated (identical calls are merged). Fusion transcripts
  are predicted and figures are drawn for visualization.
- Pairs are compared between libraries to determine if they are equivalent events. This is where somatic vs
  germline and expressed vs not expressed can be determined.
- Finally a summary report is created to make the information more accessible to the end user

|

--------------

|

Getting started
--------------------

Install
....................

Installing Dependencies
,,,,,,,,,,,,,,,,,,,,,,,,,

to install dependencies, use the requirements.txt file

.. code-block:: bash

    pip install -r requirements.txt


Running Unit Tests
,,,,,,,,,,,,,,,,,,,,

to run the tests

.. code-block:: bash

    nosetests --with-coverage --cover-html --cover-html-dir=coverage --cover-package=mavis --cover-erase


Building the documentation
,,,,,,,,,,,,,,,,,,,,,,,,,,,

can build the documentation directly using make

.. code-block:: bash

    cd docs
    make html

this will generate html documentation that can be viewed in a browser by opening the index.html file


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
| `DELLY <https://github.com/dellytools/delly>`_                             |                                                                        |
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

There are several reference files that are required for full functionality of the |TOOLNAME| pipeline

Fasta sequence file(s)
,,,,,,,,,,,,,,,,,,,,,,,

These are the sequence files in fasta format that are used in aligning and generating the fusion sequences. Found here:
`UCSC hg19 chromosome fasta sequences <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/>`_

Reference Annotations
,,,,,,,,,,,,,,,,,,,,,,,

This is a custom file format. Essentially just a tabbed or json file which contains the gene, transcript, exon, translation and protein domain positional information

.. warning::

    the :func:`~mavis.annotate.file_io.load_reference_genes` will
    only load valid translations. If the cds sequence in the annotation is not
    a multiple of :attr:`~mavis.constants.CODON_SIZE` or if a
    reference genome (sequences) is given and the cds start and end are not
    M and * amino acids as expected the translation is not loaded

Example of the json format can be seen below

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
to generate the json file using a connection to the `Ensembl <http://uswest.ensembl.org/index.html>`_ perl api.



Template metadata file
,,,,,,,,,,,,,,,,,,,,,,,,

This is the file which contains the band information for the chromosomes.
Found here: `UCSC hg19 cytoband file <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>`_.
This is only used during visualization.

.. code-block:: text

    chr1    0       2300000 p36.33  gneg
    chr1    2300000 5400000 p36.32  gpos25
    chr1    5400000 7200000 p36.31  gneg
    chr1    7200000 9200000 p36.23  gpos25
    chr1    9200000 12700000        p36.22  gneg

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

The pipeline consists of five main scripts. The usage menus for any of the scripts can be viewed by running the
script with the -h/--help option

**Example:**

.. code-block:: bash

    python sv_merge.py -h

.. toctree::
    :maxdepth: 1

    auto/sv_merge
    auto/sv_annotate
    auto/sv_validate
    auto/sv_pair
    auto/sv_summary

There are some parameters that need to be computed from the bam files. This can generally be done by running the
profile_bam.py script found in the tools directory

.. code-block:: bash

    >>> python tools/profile_bam.py /path/to/bam/file -c 16
    profiling chr 16

    FINAL
    average                    396.72
    average stdev              98.89
    median                     383
    median distrib[0.80] stdev 59.56
    median distrib[0.90] stdev 72.80
    median distrib[0.95] stdev 82.20
    median distrib[0.99] stdev 93.94
    median distrib[1.00] stdev 99.84

generally giving it a single chromosome will be enough reads but it can be given as many chromosomes/templates as
required. This script calculates the median insert size and then the standard deviation (wrt to the median not mean)
from all or a portion of the distribution of insert sizes

.. |TOOLNAME| replace:: **MAVIS**
