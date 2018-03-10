Glossary
-----------


General Terms
+++++++++++++++++++++

.. glossary::
    :sorted:
    
    flanking read pair
        A pair of reads where one read maps to one side of a set of breakpoints and its mate maps to the other.

    split read
        A read which aligns next to a breakpoint and is softclipped at one or more sides.

    spanning read
        Applies primarily to small structural variants. Reads which span both breakpoints.

    half-mapped read
        A read whose mate is unaligned. Generally this refers to reads in the evidence stage that are mapped next to a breakpoint.

    breakpoint
        A breakpoint is a genomic position (interval) on some reference/template/chromosome which has a strand and orientation. The orientation describes the portion of the reference that is retained.

    event
        Used interchangeably with :term:`structural variant`.

    event type
        Classification for a structural variant. see :term:`event_type`.

    structural variant
        A genomic alteration that can be described by a pair of breakpoints and an :term:`event type`. The two breakpoints represent regions in the genome that are broken apart and reattached together.

    breakpoint pair
        Basic definition of a :term:`structural variant`. Does not automatically imply a classification/type.

    bed
        File format specification. See https://genome.ucsc.edu/FAQ/FAQformat#format1.

    IGV batch file
        This is a file format type defined by :term:`IGV` see
        `running IGV with a batch file <https://software.broadinstitute.org/software/igv/batch>`__.

    BAM
        File format specification. See https://genome.ucsc.edu/FAQ/FAQformat#format5.1.

    2bit
        File format specification. See https://genome.ucsc.edu/FAQ/FAQformat#format7.

    fasta
        File format specification. See https://genome.ucsc.edu/FAQ/FAQformat#format18.

    psl
        File format specification. See https://genome.ucsc.edu/FAQ/FAQformat#format2.

    pslx
        Extended format of a :term:`psl`.

    SVG
        SVG (Scalable vector graph) is an image format. see https://www.w3schools.com/graphics/svg_intro.asp.

    JSON
        JSON (JavaScript Object Notation) is a data file format. see https://www.w3schools.com/js/js_json_intro.asp.

    blat
        Alignment tool. see https://genome.ucsc.edu/FAQ/FAQblat.html.

    IGV
        Integrative Genomics Viewer is a visualization tool. see http://software.broadinstitute.org/software/igv.

    SGE
        Sun Grid Engine (SGE) is a job scheduling system for cluster management see http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html.

    SLURM
        SLURM is a job scheduling system for cluster management see https://slurm.schedmd.com/quickstart.html.

    BWA
        BWA is an alignement tool. See https://github.com/lh3/bwa

    :ref:`HGVS <den-Dunnen-2016>`
        Community based standard of reccommendations for variant notation. See http://varnomen.hgvs.org/

    BreakDancer
        BreakDancer is an SV caller. Soure for BreakDancer can be found `here <https://github.com/genome/breakdancer>`__ [Chen-2009]_ 

    Chimerascan
        Chimerascan is an SV caller. Source for Chimerascan can be found `here <https://code.google.com/archive/p/chimerascan>`__ [Iyer-2011]_
    
    DeFuse
        DeFuse is an SV caller. Source for DeFuse can be found `here <https://bitbucket.org/dranew/defuse>`__ [McPherson-2011]_

    DELLY
        DELLY is an SV caller. Source for DELLY can be found `here <https://github.com/dellytools/delly>`__ [Rausch-2012]_

    Manta 
        Manta is an SV caller. Source for Manta can be found `here <https://github.com/Illumina/manta>`__ [Chen-2016]_

    Pindel
        Pindel is an SV caller. Source for Pindel can be found `here <https://github.com/genome/pindel>`__ [Ye-2009]_

    Trans-ABySS
        Trans-ABySS is an SV caller. Source for Trans-ABySS can be found `here <https://github.com/bcgsc/transabyss>`__ [Robertson-2010]_

    SV
        Structural Variant


.. include:: config_settings_glossary.rst


Column Names
++++++++++++++

List of column names and their definitions. The types indicated here are the expected types in a row for a given column name.

.. glossary::
    :sorted:

    library
        Identifier for the library/source

    cluster_id
        Identifier for the merging/clustering step

    cluster_size
        :class:`int` - The number of breakpoint pair calls that were grouped in creating the cluster

    validation_id
        Identifier for the validation step

    annotation_id
        Identifier for the annotation step

    product_id
        Unique identifier of the final fusion including splicing and ORF decision from the annotation step

    event_type
        :class:`~mavis.constants.SVTYPE` - The classification of the event

    inferred_pairing
        A semi colon delimited of event identifiers i.e. <annotation_id>_<splicing pattern>_<cds start>_<cds end>
        which were paired to the current event based on predicted products

    pairing
        A semi colon delimited of event identifiers i.e. <annotation_id>_<splicing pattern>_<cds start>_<cds end>
        which were paired to the current event based on breakpoint positions

    gene1
        Gene for the current annotation at the first breakpoint

    gene1_direction
        :class:`~mavis.constants.PRIME` - The direction/prime of the gene

    gene2
        Gene for the current annotation at the second breakpoint

    gene2_direction
        :class:`~mavis.constants.PRIME` - The direction/prime of the gene. Has the following possible values

    gene1_aliases
        Other gene names associated with the current annotation at the first breakpoint

    gene2_aliases
        Other gene names associated with the current annotation at the second breakpoint

    gene_product_type
        :class:`~mavis.constants.GENE_PRODUCT_TYPE` - Describes if the putative fusion product will be sense or anti-sense

    fusion_cdna_coding_end
        Position wrt the 5' end of the fusion transcript where coding ends last base of the stop codon

    transcript1
        Transcript for the current annotation at the first breakpoint

    transcript2
        Transcript for the current annotation at the second breakpoint

    fusion_splicing_pattern
        :class:`~mavis.constants.SPLICE_TYPE` - Type of splicing pattern used to create the fusion cDNA.

    fusion_cdna_coding_start
        :class:`int` - Position wrt the 5' end of the fusion transcript where coding begins first base of the Met amino acid.

    fusion_cdna_coding_end
        :class:`int` - Position wrt the 5' end of the fusion transcript where coding ends last base of the stop codon

    fusion_mapped_domains
        ``JSON`` - List of domains in :term:`JSON` format where each domain start and end positions are given wrt to the fusion
        transcript and the mapping quality is the number of matching amino acid positions over the total
        number of amino acids. The sequence is the amino acid sequence of the domain on the reference/original
        transcript

    fusion_sequence_fasta_id
        The sequence identifier for the cdna sequence output fasta file

    fusion_sequence_fasta_file
        ``FILEPATH`` - Path to the corresponding fasta output file

    annotation_figure
        ``FILEPATH`` - File path to the svg drawing representing the annotation

    annotation_figure_legend
        ``JSON`` - :term:`JSON` data for the figure legend

    genes_encompassed
        Applies to intrachromosomal events only. List of genes which overlap any region that occurs between both
        breakpoints. For example in a deletion event these would be deleted genes.

    genes_overlapping_break1
        list of genes which overlap the first breakpoint

    genes_overlapping_break2
        list of genes which overlap the second breakpoint

    genes_proximal_to_break1
        list of genes near the breakpoint and the distance away from the breakpoint

    genes_proximal_to_break2
        list of genes near the breakpoint and the distance away from the breakpoint

    break1_chromosome
        :class:`str` - The name of the chromosome on which breakpoint 1 is situated

    break1_position_start
        :class:`int` - Start integer inclusive 1-based of the range representing breakpoint 1

    break1_position_end
        :class:`int` - End integer inclusive 1-based of the range representing breakpoint 1

    break1_orientation
        :class:`~mavis.constants.ORIENT` - The side of the breakpoint wrt the positive/forward strand that is retained.

    break1_strand
        :class:`~mavis.constants.STRAND` - The strand wrt to the reference positive/forward strand at this breakpoint.

    break1_seq
        :class:`str` - The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand

    break2_chromosome
        The name of the chromosome on which breakpoint 2 is situated

    break2_position_start
        :class:`int` - Start integer inclusive 1-based of the range representing breakpoint 2

    break2_position_end
        :class:`int` - End integer inclusive 1-based of the range representing breakpoint 2

    break2_orientation
        :class:`~mavis.constants.ORIENT` - The side of the breakpoint wrt the positive/forward strand that is retained.

    break2_strand
        :class:`~mavis.constants.STRAND` - The strand wrt to the reference positive/forward strand at this breakpoint.

    break2_seq
        :class:`str` - The sequence up to and including the breakpoint. Always given wrt to the positive/forward strand

    opposing_strands
        :class:`bool` - Specifies if breakpoints are on opposite strands wrt to the reference. Expects a boolean

    stranded
        :class:`bool` - Specifies if the sequencing protocol was strand specific or not. Expects a boolean

    protocol
        :class:`~mavis.constants.PROTOCOL` - Specifies the type of library

    tools
        The tools that called the event originally from the cluster step. Should be a semi-colon delimited list of
        <tool name>_<tool version>

    contigs_assembled
        :class:`int` - Number of contigs that were built from split read sequences

    contigs_aligned
        :class:`int` - Number of contigs that were able to align

    contig_alignment_query_name
        The query name for the contig alignment. Should match the 'read' name(s) in the .contigs.bam output file

    contig_seq
        :class:`str` - Sequence of the current contig wrt to the positive forward strand if not strand specific

    contig_remap_score
        :class:`float` - Score representing the number of sequences from the set of sequences given to the assembly
        algorithm that were aligned to the resulting contig with an acceptable scoring based on user-set thresholds.
        For any sequence its contribution to the score is divided by the number of mappings to give less weight to
        multimaps

    call_sequence_complexity
        :class:`float` - The minimum amount any two bases account for of the proportion of call sequence. An average for non-contig calls

    contig_remapped_reads
        :class:`int` - the number of reads from the input bam that map to the assembled contig

    contig_remapped_read_names
        read query names for the reads that were remapped. A -1 or -2 has been appended to the end of the name to
        indicate if this is the first or second read in the pair

    contig_alignment_score
        :class:`float` - A rank based on the alignment tool blat etc. of the alignment being used. An average if
        split alignments were used. Lower numbers indicate a better alignment. If it was the best alignment possible
        then this would be zero.

    contig_alignment_reference_start
        The reference start(s) <chr>:<position> of the contig alignment. Semi-colon delimited

    contig_alignment_cigar
        The cigar string(s) representing the contig alignment. Semi-colon delimited

    contig_remap_coverage
        :class:`float` - Fraction of the contig sequence which is covered by the remapped reads

    contig_build_score
        :class:`int` - Score representing the edge weights of all edges used in building the sequence

    contig_strand_specific
        :class:`bool` - A flag to indicate if it was possible to resolve the strand for this contig

    spanning_reads
        :class:`int` - the number of spanning reads which support the event

    spanning_read_names
        read query names of the spanning reads which support the current event

    call_method
        :class:`~mavis.constants.CALL_METHOD` - The method used to call the breakpoints

    flanking_pairs
        :class:`int` - Number of read-pairs where one read aligns to the first breakpoint window and the second read
        aligns to the other. The count here is based on the number of unique query names

    flanking_pairs_compatible
        :class:`int` - Number of flanking pairs of a compatible orientation type. This applies to insertions and
        duplications. Flanking pairs supporting an insertion will be compatible to a duplication and flanking pairs
        supporting a duplication will be compatible to an insertion (possibly indicating an internal translocation)

    flanking_median_fragment_size
        :class:`int` - The median fragment size of the flanking reads being used as evidence

    flanking_stdev_fragment_size
        :class:`float` - The standard deviation in fragment size of the flanking reads being used as evidence

    break1_split_reads
        :class:`int` - Number of split reads that call the exact breakpoint given

    break1_split_reads_forced
        :class:`int` - Number of split reads which were aligned to the opposite breakpoint window using a targeted
        alignment

    break2_split_reads
        :class:`int` - Number of split reads that call the exact breakpoint given

    break2_split_reads_forced
        :class:`int` - Number of split reads which were aligned to the opposite breakpoint window using a targeted
        alignment

    linking_split_reads
        :class:`int` - Number of split reads that align to both breakpoints

    untemplated_seq
        :class:`str` - The untemplated/novel sequence between the breakpoints

    break1_homologous_seq
        :class:`str` - Sequence in common at the first breakpoint and other side of the second breakpoint

    break2_homologous_seq
        :class:`str` - Sequence in common at the second breakpoint and other side of the first breakpoint

    break1_ewindow
        ``int-int`` - Window where evidence was gathered for the first breakpoint

    break1_ewindow_count
        :class:`int` - Number of reads processed/looked-at in the first evidence window

    break1_ewindow_practical_coverage
        :class:`float` - break2_ewindow_practical_coverage, break1_ewindow_count / len(break1_ewindow). Not the actual
        coverage as bins are sampled within and there is a read limit cutoff

    break2_ewindow
        ``int-int`` - Window where evidence was gathered for the second breakpoint

    break2_ewindow_count
        :class:`int` - Number of reads processed/looked-at in the second evidence window

    break2_ewindow_practical_coverage
        :class:`float` - break2_ewindow_practical_coverage, break2_ewindow_count / len(break2_ewindow). Not the actual
        coverage as bins are sampled within and there is a read limit cutoff

    raw_flanking_pairs
        :class:`int` - Number of flanking reads before calling the breakpoint. The count here is based on the number of
        unique query names

    raw_spanning_reads
        :class:`int` - Number of spanning reads collected during evidence collection before calling the breakpoint

    raw_break1_split_reads
        :class:`int` - Number of split reads before calling the breakpoint

    raw_break2_split_reads
        :class:`int` - Number of split reads before calling the breakpoint

    cdna_synon
        semi-colon delimited list of transcript ids which have an identical cdna sequence to the cdna sequence of the
        current fusion product

    protein_synon
        semi-colon delimited list of transcript ids which produce a translation with an identical amino-acid sequence
        to the current fusion product

    tracking_id
        column used to store input identifiers from the original SV calls. Used to track calls from the input files to
        the final outputs.

    fusion_protein_hgvs
        :class:`str` - Describes the fusion protein in HGVS notation. Will be None if the change is not an indel or is synonymous

    net_size
        ``int-int`` - The net size of an event. For translocations and inversion this will always be 0. For indels it will be
        negative for deletions and positive for insertions. It is a range to accommodate non-specific events.

    supplementary_call
        :class:`bool` - Flag to indicate if the current event was a supplementary call, meaning a call that was found as a result of
        validating another event.
