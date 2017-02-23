
Theory and Models
--------------------

Structural Variants in this tool are defined as a pair of breakpoints. A breakpoint is a genomic position
(interval) on some reference/template/chromosome which has a :term:`strand` and :term:`orientation`.
The orientation describes the portion of the reference that is retained.

.. figure:: _static/example_figure.svg
    :width: 100%

    Example output from the tool visulaizing an inverted translocation. (c1) The first template, chromosome 1. (cX) The second template
    chromosome X. (B1) The first breakpoint has a left orientation and retains the five prime portions of the gene, HUGO2. (B2) The
    second breakpoint also has a left orientation but retains the three prime portion of the gene, HUGO3. (T1) The original transcript
    for the first gene. (T2) The original transcript for the second gene. (F1) The fusion transcript.
|

------

|

Paired-end Reads: Flanking evidence
......................................

One of the most confusing parts about working with contig and paired-end reads is relating them to the
breakpoint so that you can determine which types will support an event. The flanking read types we outline here
are similarly described by `IGV <http://software.broadinstitute.org/software/igv/interpreting_insert_size>`_.
We have used similar coloring for the read pairs in the following diagrams to
facilitate ease of use for those already familiar with viewing bam files in IGV.

.. note::

    The major assumptions here are that the 'normal' read-pair is a read pair which has one read on the positive/forward
    strand and its partner on the negative/reverse strand. It is assumed that partners share a read name, as is the case for Illumina reads.


Deletion
,,,,,,,,,

For a deletion, we expect the flanking reads to be in the normal orientation but that the
fragment size should be abnormal (for large deletions).

.. figure:: _static/read_pairs_deletion.svg
    :width: 100%

    Flanking read pair evidence for a deletion event. the read pairs will have a larger than expected fragment size when mapped to the
    reference genome because in the mutant genome they are closer together, owing to the deletion event. (B1) The first breakpoint
    which has a left orientation. (B2) The second breakpoint which has a right orientation. Both breakpoints would be on the positive
    strand (assuming that the input is stranded) which means that the first read in the pair would be on the positive strand and the
    second read in the pair would be on the negative/reverse strand.



Insertion
,,,,,,,,,,

.. figure:: _static/read_pairs_insertion.svg
    :width: 100%

    Flanking read pair evidence for an insertion event. The read pairs will have a smaller than expected fragment size when mapped to the
    reference genome because in the mutant genome they are father apart, owing to the insertion event. (B1) The first breakpoint
    which has a left orientation. (B2) The second breakpoint which has a right orientation. Both breakpoints would be on the positive
    strand (assuming that the input is stranded) which means that the first read in the pair would be on the positive strand and the
    second read in the pair would be on the negative/reverse strand.



Duplication
,,,,,,,,,,,,,

.. figure:: _static/read_pairs_duplication.svg
    :width: 100%

    Flanking read pair evidence for a tandem duplication event. The read pairs will have an abnormal orientation but still the
    same strands as the normal read pair. (B1) The first breakpoint will be on the positive strand and have a right orientation.
    (B2) The second breakpoint will be on the positive strand and have a left orientation.



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

.. figure:: _static/read_pair_definitions.svg
    :width: 100%

    Basic Terms used in describing read pairs are shown above: fragment size: the distance between the pair;
    read length: the length of the read; fragment size: the combined length of both reads and the fragment size

We make some base assumptions with regards to paired-end read data:

.. note::

    the distribution of fragment sizes approximately follows a normal distribution

.. note::

    the most common fragment size is the unmutated 'normal' fragment

With the above assumptions we take the median fragment size to be the expected normal.

Given that we expect mutations and therefore abnormal fragment sizes, we use a modified method to calculate the
**median standard deviation** (*s* in the equations below). We calculate the squared distance away from the 
median for each fragment and then take a fraction of this to be 'normal' variation. So the most abnormal portion is
ignored, assuming it is supposed to be abnormal. This results in a calculation as follows, where the original set
Y is the set of fragment sizes from the bam file and f is the fraction of fragment sizes assumed to be normal

.. code::

    from statistics import median
    import math

    fragments = [abs(read.template_length) for read in reads]  # the fragment sizes of the reads
    f = 0.95 # fraction
    m = median(fragments) # get the median fragment size value
    X = [math.pow(i - m, 2) for i in fragments]  # take the square error for each point
    end = int(round(len(X) * f))
    X = sorted(X)[0:end]
    stdev = math.sqrt(sum(X) / len(X))

Using the above calculations we can generate a modified version of the standard deviation (s above) as shown in the figure
below (stdev). This gives us an idea of when to judge an fragment size as abnormal and where we expect our normal read
pairs fragment sizes to fall.

.. figure:: _static/insert_size_distrb_fractions.svg
    :width: 100%

    Distribution of fragment sizes (absolute values) of proper read pairs, and different normal distribution fits using
    the above equation. The different coloured curves are computed with different parameters. black: the standard
    calculation using all data points and the mean as centre; dark green: median as centre and a fraction of f=0.80;
    light green: median as centre, f=0.90; light blue: median and f=0.95; dark blue: median and f=1.00.

As we can see from the distribution above the median approximates the distribution centre better than the mean,
likely because it is more resistant to outliers.

.. figure::  _static/insert_size_distrb.svg
    :width: 100%

    Distribution of fragment sizes (absolute values) of proper read pairs. In the above image the standard deviation
    (stdev) was calculated with respect to the median (383) using the fraction (f=0.99).



We use this in two ways

1. to find flanking evidence supporting deletions and insertions
2. to estimate the window size for where we will need to read from the bam when looking for evidence for a given event

The :py:func:`~structural_variant.validate.Evidence.generate_window` function uses the above concepts. The user will
define the :py:attr:`~structural_variant.validate.EvidenceSettings.median_fragment_size` the
:py:attr:`~structural_variant.validate.EvidenceSettings.stdev_fragment_size`, and the
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

The following decision tree is used in classifying events based on their breakpoints. Only valid combinations have
been shown.

.. figure:: _static/classification_tree.svg
    :width: 100%

    Classification Decision Tree. The above  diagram details the decsion logic for classifying events based on the
    orientation, strand and chromosomes or their respective breakpoints

|

-----------------

|

Assembling Contigs
......................

During validation, for each breakpoint pair, we attempt to assemble a contig to represent the sequence across
the breakpoints. This is assembled from the :term:`split reads` and mates of :term:`half-mapped` reads that have been
collected. The assembly uses a :term:`DeBruijn graph`.

Breakpoints can be called by multiple different :attr:`~structural_variant.constants.CALL_METHOD`.

|

------

|

Calling Breakpoints by Flanking Evidence
..........................................

Breakpoints are called by contig, split-read, or flanking pairs evidence. Contigs and split reads are used to call exact
breakpoints, where breakpoints called by flanking reads are generally assigned a probabalistic range.

The metrics used here are similar to those used in calculating the evidence window. We use the
:term:`max_expected_fragment_size` as the outer limit of how large the range can be. This is further refined taking
into account the range spanned by the :term:`flanking pairs` evidence and the position of the opposing breakpoint.

.. figure:: _static/call_breakpoint_by_flanking_reads.svg
    :width: 100%

    Calculation of the left-oriented breakpoint by flanking reads. Reads mapped to the breakpoint are shown in grey.
    The read on the right (black outline, no fill) demonstrates the read length used to narrow the right side bound of
    the estimated breakpoint interval.


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

Next there are the fusion-product level annotations. If the event result in a fusion
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
    breakpoint is a donor and the first donor is used. (D) The first splice site following the breakpoint is an acceptor.
    (E-H) The breakpoint lands in an exon and the three prime portion of the transcript is retained. (E) The original
    splicing pattern showing the placement of the genomic breakpoint and the retained three prime portion. (F) The first
    splice site prior to the breakpoint is an acceptor and the first acceptor is used. (G) The first splice site prior to the
    breakpoint is an acceptor and the second acceptor is used. (H) The first splice site prior to the breakpoint is a donor

|

------

|

Glossary
.............

..  glossary::
	:sorted:

    flanking read pair
        a pair of reads where one read maps to one side of a set of breakpoints and its mate maps to the other

    split read
        a read which aligns next to a breakpoint and is softclipped at one or more sides

    spanning read
        applies primarily to small structural variants. Reads which span both breakpoints

    half-mapped read
        a read whose mate is unaligned. Generally this refers to reads in the evidence stage that are mapped next to a breakpoint.

    breakpoint
         A breakpoint is a genomic position (interval) on some reference/template/chromosome which has a strand and orientation. The orientation describes the portion of the reference that is retained.

    event
        used interchangeably with :term:`structural variant`

    event type
        classification for a structural variant. see :term:`event_type`

    structural variant
        a genomic alteration that can be described by a pair of breakpoints and an :term:`event type`. The two breakpoints represent regions in the genome that are broken apart and reattached together.

    breakpoint pair
        :term:`structural variant` which has not been classified/given a type




.. |TOOLNAME| replace:: **MARVIN**
