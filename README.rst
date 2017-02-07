|TOOLNAME| User Manual
========================

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

    nosetests --with-coverage --cover-html --cover-html-dir=coverage --cover-package=structural_variant --cover-erase


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


Available Pre-formatting scripts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

- `convert_ta.py <https://svn.bcgsc.ca/svn/SVIA/sv_compile/tags/0.0.1/tools/convert_ta.py>`_


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

    the :func:`~structural_variant.annotate.file_io.load_reference_genes` will
    only load valid translations. If the cds sequence in the annotation is not
    a multiple of :attr:`~structural_variant.constants.CODON_SIZE` or if a
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

This reference file can be generated from any database with the necessary information. There is a basic perl script to generate the json file using a connection to the `Ensembl <http://uswest.ensembl.org/index.html>`_ perl api.



Template metadata file
,,,,,,,,,,,,,,,,,,,,,,,,

This is the file which contains the band information for the chromosomes.
Found here: `UCSC hg19 cytoband file <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz>`_.
This is only used during visualization.


|

------

|



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

|
|

-----------------

|
|

Theory and Models
--------------------

Structural Variants in this tool are defined as a pair of breakpoints. A breakpoint is a genomic position
(interval) on some reference/template/chromosome which has a :term:`strand` and :term:`orientation`.
The orientation describes the portion of the reference that is retained.

.. figure:: _static/example_figure.svg
    :width: 100%

    Example output from the tool visulaizing an inverted translocation. (c1) the first template, chromosome 1. (cX) the second template
    chromosome X. (B1) the first breakpoint has a left orientation and retains the five prime portions of the gene, HUGO2. (B2) the
    second breakpoint also has a left orientation but retains the three prime portion of the gene, HUGO3. (T1) the original transcript
    for the first gene. (T2) the original transcript for the second gene. (F1) the fusion transcript
|

------

|

Paired-end Reads: Flanking evidence
......................................

One of the most confusing parts about working with :term:`contig` and paired-end reads is relating them to the
breakpoint so that you can determine which types will support an event. The flanking read types we outline here
are similarly described by `IGV <http://software.broadinstitute.org/software/igv/interpreting_insert_size>`_.
We have used similar coloring for the read pairs in the following diagrams to
facilitate ease of use for those already familiar with viewing bam files in IGV.

.. note::

    The major assumptions here are that the 'normal' read-pair is a read pair which has one read on the positive/forward
    strand and its partner on the negative/reverse strand. It is assumed that partners share a read name. As is the case for illumina reads


Deletion
,,,,,,,,,

For a deletion, we expect the flanking reads to be in the normal orientation but that the
insert size should be abnormal (for large deletions).

.. figure:: _static/read_pairs_deletion.svg
    :width: 100%

    Flanking read pair evidence for a deletion event. the read pairs will have a larger than expected insert size when mapped to the
    reference genome because in the mutant genome they are closer together, owing to the deletion event. (B1) the first breakpoint
    which has a left orientation (B2) the second breakpoint which has a right orientation. Both breakpoints would be on the positive
    strand (assuming that the input is stranded) which means that the first read in the pair would be on the positive strand and the
    second read in the pair would be on the negative/reverse strand.



Insertion
,,,,,,,,,,

.. figure:: _static/read_pairs_insertion.svg
    :width: 100%

    Flanking read pair evidence for an insertion event. the read pairs will have a smaller than expected insert size when mapped to the
    reference genome because in the mutant genome they are father apart, owing to the insertion event. (B1) the first breakpoint
    which has a left orientation (B2) the second breakpoint which has a right orientation. Both breakpoints would be on the positive
    strand (assuming that the input is stranded) which means that the first read in the pair would be on the positive strand and the
    second read in the pair would be on the negative/reverse strand.



Duplication
,,,,,,,,,,,,,

.. figure:: _static/read_pairs_duplication.svg
    :width: 100%

    Flanking read pair evidence for a tandem duplication event. the read pairs will have an abnormal orientation but still the
    same strands as the normal read pair. (B1) the first breakpoint will be on the positive strand and have a right orientation.
    (B2) the second breakpoint will be on the positive strand and have a left orientation.



Inversion
,,,,,,,,,,

.. figure:: _static/read_pairs_inversion_LL.svg
    :width: 100%

    Flanking read pair evidence for an inversion. Both breakpoints have left orientation.

.. figure:: _static/read_pairs_inversion_RR.svg
    :width: 100%

    Flanking read pair evidence for an inversion. Both breakpoints have right orientation.



Translocation
,,,,,,,,,,,,,,

.. figure:: _static/read_pairs_translocation_LR.svg
    :width: 100%

    Flanking read pair evidence for a translocation. (B1) the first breakpoint with a left orientation. (B2) the second breakpoint with a right orientation.

.. figure:: _static/read_pairs_translocation_RL.svg
    :width: 100%

    Flanking read pair evidence for a translocation. (B1) the first breakpoint with a right orientation. (B2) the second breakpoint with a left orientation.

Inverted Translocation
,,,,,,,,,,,,,,,,,,,,,,,,

.. figure:: _static/read_pairs_translocated_inversion_LL.svg
    :width: 100%

    Flanking read pair evidence for an inverted translocation. Both breakpoints have left orientation.

.. figure:: _static/read_pairs_translocated_inversion_RR.svg
    :width: 100%

    Flanking read pair evidence for an inverted translocation. Both breakpoints have right orientation.


|

------

|

Calculating the Evidence Window
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

.. code::

    from statistics import median
    import math

    inserts = [abs(read.template_length) for read in reads]  # the insert sizes of the reads
    f = 0.95 # fraction
    m = median(inserts) # get the median insert size value
    X = [math.pow(i - m, 2) for i in inserts]  # take the square error for each point
    end = int(round(len(X) * f))
    X = sorted(X)[0:end]
    stdev = math.sqrt(sum(X) / len(X))

Using the above calculations we can generate a modified version of the standard deviation (s above) as shown in the figure
below (stdev). This gives us an idea of when to judge an insert size as abnormal and where we expect our normal read
pairs insert sizes to fall.

.. figure:: _static/insert_size_distrb_fractions.svg
    :width: 100%

    Distribution of insert sizes (absolute values) of proper read pairs, and different normal distribution fits using
    the above equation. The different coloured curves are computed with different parameters. black: the standard
    calculation using all data points and the mean as centre; dark green: median as centre and a fraction of f=0.80;
    light green: median as centre, f=0.90; light blue: median and f=0.95; dark blue: median and f=1.00.

As we can see from the distribution above the median approximates the distribution centre better than the mean,
likely because it is more resistant to outliers.

.. figure::  _static/insert_size_distrb.svg
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


|

-----------------

|

Classifying Events
.....................

the following decision tree is used in classifying events based on their breakpoints. Only valid combinations have
been shown

.. figure:: _static/classification_tree.svg
    :width: 100%

    Classification Decision Tree. The above  diagram details the decsion logic for classifying events based on the
    orientation, strand and chromosomes or their respective breakpoints

|

-----------------

|

Assembling Contigs
......................

During validation, for each breakpoint pair, we attempt to assemble a :term:`contig` to represent the sequence across
the breakpoints. This is assembled from the :term:`split reads` and mates of :term:`half-mapped` reads that have been
collected. The assembly uses a :term:`DeBruijn graph`.

Breakpoints can be called by multiple different :attr:`~structural_variant.constants.CALL_METHOD`.

|

-----------------

|

Annotating Events
....................

We make the following assumptions when determining the annotations for each event

.. note::

    If both breakpoints are in the same gene, they must also be in the same transcript

.. note::

    If the breakpoint intervals overlap we do not annotate encompassed genes

.. note::

    Encompassed and 'nearest' genes are reported without respect to strand

There are specific questions we want annotation to answer. We collect gene level
annotations which describes things like what gene is near the breakpoint (useful
in the case of a potential promoter swap); what genes (besides the one selected)
also overlap the breakpoint; what genes are encompassed between the breakpoints
(for example in a deletion event the genes that would be deleted).

.. figure:: _static/annotations_summary.svg
    :width: 100%

    Gene level annotations at each breakpoint. Note: genes which fall
    between a breakpoint pair, encompassed genes, will not be present for
    interchromosomal events (translocations)

Next there are the fusion-product level annotations. If the event result ina fusion
transcript, the sequence of the fusion transcript is computed. This is translated
to a putative amino acid sequence from which protein metrics such as the possible
ORFs and domain sequences can be computed.

|

-----------------

|

Predicting Splicing Patterns
............................

After the events have been called and an annotation has been attached, we often want to predict information about the
putative fusion protein, which may be a product. In some cases, when a fusion transcript disrupts a splice-site, it is
not clear what the processed fusion transcript may be. |TOOLNAME| will calculate all possibilities according to the
following model.


.. figure:: _static/splicing_pattern_default.svg
    :width: 100%

    The default splicing pattern is a list of pairs of donor and acceptor splice sites

For a given list of non-abrogated splice sites (listed 5' to 3' on the strand of the transcript) donor splice sites
are paired with all following as seen below

.. figure:: _static/splicing_pattern_multiple_donors.svg
    :width: 100%

    Multiple abrogated acceptors sites. As one can see above this situation will result in 3 different splicing
    patterns depending on which donor is paired with the 2nd acceptor site

.. figure:: _static/splicing_pattern_multiple_acceptors.svg
    :width: 100%

    Multiple abrogated donor sites. As one can see above this situation will result in 3 different splicing
    patterns depending on which acceptor is paired with the 2nd donor site


More complex examples are drawn below. There are five classifications (:class:`~structural_variant.constants.SPLICE_TYPE`) for the different splicing patterns:

1. Retained intron (:class:`~structural_variant.constants.SPLICE_TYPE.RETAIN`)
2. Skipped exon (:attr:`~structural_variant.constants.SPLICE_TYPE.SKIP`)
3. Multiple retained introns (:attr:`~structural_variant.constants.SPLICE_TYPE.MULTI_RETAIN`)
4. Multiple skipped exons (:attr:`~structural_variant.constants.SPLICE_TYPE.MULTI_SKIP`)
5. Some combination of retained introns and skipped exons (:attr:`~structural_variant.constants.SPLICE_TYPE.COMPLEX`)

.. figure:: _static/splicing_model.svg
    :width: 100%

    Splicing scenarios

|

-----------------

|

Pairing Similar Events
.......................

After breakpoints have been called and annotated we often need to see if the same event was found in different samples. To do this we will need to compare events between genome and transcriptome libraries. For this, the following model is proposed. To compare events between different protocol (genome vs transcriptome) we use the annotation overlying the genome breakpoint and the splicing model we defined above to predict where we would expect to find the transcriptomic breakpoints. This gives rise to the following basic cases.

.. note::
    In all cases the predicted breakpoint is either the same as the genomic breakpoint, or it is the same as the nearest retained donor/acceptor to the breakpoint.

.. figure:: _static/breakpoint_prediction_exonic.svg
    :width: 100%

    (A-D) The breakpoint lands in an exon and the five prime portion of the transcript is retained. (A) The original
    splicing pattern showing the placement of the genomic breakpoint and the retained five prime portion. (B) The first
    splice site following the breakpoint is a donor and the second donor is used. (C) The first splice site following the
    breakpoint is a donor and the first donor is used. (D) The first slice site following the breakpoint is an acceptor. 
    (E-H) The breakpoint lands in an exon and the three prime portion of the transcript is retained. (E) The original
    splicing pattern showing the placement of the genomic breakpoint and the retained three prime portion. (F) The first
    splice site prior to the breakpoint is an acceptor and the first acceptor is used. (G) The first splice site prior to the
    breakpoint is an acceptor and the second acceptor is used. (H) The first slice site prior to the breakpoint is a donor



|
|

-----------------

|
|


Development
--------------------------

Guidelines for Contributors
.............................

- In general, follow `pep8 <https://www.python.org/dev/peps/pep-0008/>`_ style guides using a maximum line width of 120 characters
- docstrings should follow `sphinx google code style <http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_
- any column name which may appear in any of the intermediate or final output files must be defined in :class:`~structural_variant.constants.COLUMNS`

TODO
...........................

.. todolist::

.. |TOOLNAME| replace:: **SVMerge**
