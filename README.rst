SVMerge User Manual
====================

pipeline to merge and validate input from different structural variant callers into a single report


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

    nosetests --with-coverage --cover-html --cover-html-dir=coverage --cover-package=structural_variant --cover-erase


Building the documentation
,,,,,,,,,,,,,,,,,,,,,,,,,,,

can build the documentation directly using make

.. code-block:: bash

    cd docs
    make html

this will generate html documentation that can be viewed in a browser by opening the index.html file


Input Files
....................

The requirements are described in `JIRA <https://www.bcgsc.ca/jira/browse/APA-618>`_ and are listed below.
These pertain to the input files from the various tools you want to merge

::

    | column name       | description                                                                       |
    |-------------------|-----------------------------------------------------------------------------------|
    | start_position    | range for the first breakpoint position                                           |
    | start_strand      | the reference strand aligned to                                                   |
    | start_orientation | the orientation (L or R) retained at the first breakpoint wrt the positive strand |
    | end_chromosome    |                                                                                   |
    | end_position      | range for the second breakpoint position                                          |
    | end_strand        | the reference strand aligned to                                                   |
    | end_orientation   | the orientation (L or R) retained at the second breakpoint wrt the positive strand|
    | protocol          | genome or transcriptome                                                           |
    | library           |                                                                                   |
    | tool_version      | tool name and tool version number joined with an underscore, no spaces            |
    | opposing_strand   | boolean value to describe if the breakpoints are on opposing strands              |
    | stranded          | boolean, if True then read1 in read pairs is assume to be the proper strand       |

Available Pre-formatting scripts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

- `convert_ta.py <https://svn.bcgsc.ca/svn/SVIA/sv_compile/tags/0.0.1/tools/convert_ta.py>`_

Reference Files
..................

There are several reference files that are required for full functionality of the svmerge pipeline

Fasta sequence file(s)
,,,,,,,,,,,,,,,,,,,,,,,

These are the sequence files in fasta format that are used in aligning and generating the fusion sequences. For hg19 these files can be found here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr<#>.fa.gz

Reference Annotations
,,,,,,,,,,,,,,,,,,,,,,,

This is a custom file format. essentially just a tabbed file which contains the gene, transcript, exon, translation and protein domain positional information

Template metadata file
,,,,,,,,,,,,,,,,,,,,,,,,

This is the file which contains the band information for the chromosomes. This is only used in drawing templates which is an optional step. Therefore this file is also optional. You can file the file for hg19 here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz


Running the Pipeline
.....................

The pipeline consists of three main scripts. The usage menus for any of the scripts can be viewed by running the
script with the -h/--help option

**Example:**

.. code-block:: bash

    python sv_merge.py -h

.. toctree::
    :maxdepth: 1

    sv_merge
    sv_annotate
    sv_validate


Theory and Models
--------------------

A breakpoint is defined by the reference template (i.e. chromosome), position (or range of positions) on the template,
:term:`orientation`, and :term:`strand`.

One of the most confusing parts about working with :term:`contig` and paired-end reads is relating them to the
breakpoint so that you can determine which types will support an event. For convenience We have shown the expected
:term:`strand` and :term:`orientation` of both :term:`contig` and read-pair supporting evidence side-by-side for the
major event types

.. figure:: _static/svmerge_read_pairs_vs_contigs_evidence.svg
    :width: 100%

Gathering evidence from the bam file
......................................

we make some base assumptions with regards to paired-end read data

.. note::

    the distribution of insert sizes approximately follows a normal distribution

.. note::

    the most common insert size is the unmutated 'normal' fragment

with the above assumptions we take the median insert size to be the expected normal

Given that we expect mutations and therefore abnormal insert sizes we use a modified method to calculate the
**median standard deviation** (*s* in the equations below). We calculate the squared distance for each fragment
away from the median and then take a fraction of this to be 'normal' variation. So the most abnormal portion is
ignored, assuming it is supposed to be abnormal. This results in a calculation as follows, where the original set
Y is the set of insert sizes from the bam file and f is the fraction of insert sizes assumed to be normal

.. math::

    l = |Y| \\
    y_m = \text{median value of Y}\\

    X = \left\{ x_i \mid x_i = (y_i - y_m)^2 \mid x_i \leq x_{i+1}\right\} \\

    s = \sqrt{\sum_{i=0}^{||l \cdot f||}{x_i}}

Using the above equation we can generate a modified version of the standard deviation (s above) as shown in the figure
below (stdev). This gives us an idea of when to judge an insert size as abnormal and where we expect our normal read
pairs insert sizes to fall.

.. figure:: _static/svmerge_insert_size_distrb_fractions.svg
    :width: 100%

    Distribution of insert sizes (absolute values) of proper read pairs, and different normal distribution fits using
    the above equation. The different coloured curves are computed with different parameters. black: the standard
    calculation using all data points and the mean as centre; dark green: median as centre and a fraction of f=0.80;
    light green: median as centre, f=0.90; light blue: median and f=0.95; dark blue: median and f=1.00.

As we can see from the distribution above the median approximates the distribution centre better than the mean,
likely because it is more resistant to outliers.

.. figure::  _static/svmerge_insert_size_distrb.svg
    :width: 100%

    Distribution of insert sizes (absolute values) of proper read pairs. In the above image the standard deviation
    (stdev) was calculated with respect to the median (383) using the fraction (f=0.99).



We use this in two ways

1. to find flanking evidence supporting deletions and insertions
2. to estimate the window size for where we will need to read from the bam when looking for evidence for a given event

The :py:func:`~structural_variant.validate.Evidence.generate_window` function uses the above concepts. The user will
define the :py:attr:`~structural_variant.validate.EvidenceSettings.median_insert_size` the
:py:attr:`~structural_variant.validate.EvidenceSettings.tdev_isize`, and the
:py:attr:`~structural_variant.validate.EvidenceSettings.stdev_count_abnormal` parameters defined in the
:class:`~structural_variant.validate.EvidenceSettings` class.

If the library has a transcriptome protocol this becomes a bit more complicated and we must take into account the
possible annotations when calculating the evidence window. see
:py:func:`~structural_variant.validate.Evidence.generate_transcriptome_window` for more

Classifying Events
.....................

the following decision tree is used in classifying events based on their breakpoints. Only valid combinations have
been shown

.. figure:: _static/svmerge_classification_tree.svg
    :width: 100%

    Classification Decision Tree. The above  diagram details the decsion logic for classifying events based on the
    orientation, strand and chromosomes or their respective breakpoints


Assembling Contigs
......................

During validation, for each breakpoint pair, we attempt to assemble a :term:`contig` to represent the sequence across
the breakpoints. This is assembled from the :term:`split reads` and mates of :term:`half-mapped` reads that have been
collected. The assembly uses a :term:`DeBruijn graph`.

Breakpoints can be called by multiple different :attr:`~structural_variant.constants.CALL_METHOD`.


Breakpoint sequence homology
..............................


Annotation
....................

We make the following assumptions when determining the annotations for each event

.. note::

    If both breakpoints are in the same gene, they must also be in the same transcript

.. note::

    If the breakpoint intervals overlap we do not annotate encompassed genes

.. note::

    Encompassed and 'nearest' genes are reported without respect to strand

There are specific question we want annotation to answer

- does the event result in a novel fusion transcript?
    - does this result in a protein?
    - does it retain the frame of the original proteins?
    - where does the breakpoint land in the original transcript?
    - if it is not in-frame where is the new truncation?
- what are the nearest genes outside the event (promoter swap?)
- what genes are encompassed within the event?


Splicing Model
.....................

After the events have been called and an annotation has been attached, we often want to predict information about the
putative fusion protein, which may be a product. In some cases, when a fusion transcript disrupts a splice-site, it is
not clear what the processed fusion transcript may be. SVMerge will calculate all possibilities according to the
following model.

.. figure:: _static/svmerge_splicing_model.svg
    :width: 100%

    Putative splicing scenarios. (A) a five-prime and the next three-prime splice sites are lost. (B) A five-prime
    splice site is lost. This brings about two splicing possibilities. Either the exon is skipped or the exon and
    proximal intron are retained. (C) A three-prime splice site is lost. (D) A three-prime splice site, and the next
    five-prime splice sites are lost.


Development
--------------------------

Guidelines for Contributors
.............................

- In general, follow pep8 style guides using a maximum line width of 120 characters
- docstrings should follow sphinx google code style as seen here http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
- any column name which may appear in any of the intermediate or final output files must be defined in :class:`~structural_variant.constants.COLUMNS`

TODO
...........................

.. todolist::


- Pairing events between libraries
