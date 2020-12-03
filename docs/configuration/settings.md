

# Configurable Settings
## aligner

**type**: `#!python mavis.align.SUPPORTED_ALIGNER`

**environment variable**: `MAVIS_ALIGNER`

**default**: `#!python 'blat'`

**accepted values**: `'bwa mem'`, `'blat'`


The aligner to use to map the contigs/reads back to the reference e.g blat or bwa
        

## aligner_reference

**type**: `#!python filepath`

**environment variable**: `MAVIS_ALIGNER_REFERENCE`

**default**: `#!python None`

Path to the aligner reference file used for aligning the contig sequences
        

## annotation_filters

**type**: `#!python str`

**environment variable**: `MAVIS_ANNOTATION_FILTERS`

**default**: `#!python 'choose_more_annotated,choose_transcripts_by_priority'`

A comma separated list of filters to apply to putative annotations
        

## annotation_memory

**type**: `#!python int`

**environment variable**: `MAVIS_ANNOTATION_MEMORY`

**default**: `#!python 12000`

Default memory limit (mb) for the annotation stage
        

## annotations

**type**: `#!python filepath`

**environment variable**: `MAVIS_ANNOTATIONS`

**default**: `#!python []`

Path to the reference annotations of genes, transcript, exons, domains, etc
        

## assembly_kmer_size

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_ASSEMBLY_KMER_SIZE`

**default**: `#!python 0.74`

The percent of the read length to make kmers for assembly
        

## assembly_max_paths

**type**: `#!python int`

**environment variable**: `MAVIS_ASSEMBLY_MAX_PATHS`

**default**: `#!python 8`

The maximum number of paths to resolve. this is used to limit when there is a messy assembly graph to resolve. the assembly will pre-calculate the number of paths (or putative assemblies) and stop if it is greater than the given setting
        

## assembly_min_edge_trim_weight

**type**: `#!python int`

**environment variable**: `MAVIS_ASSEMBLY_MIN_EDGE_TRIM_WEIGHT`

**default**: `#!python 3`

This is used to simplify the debruijn graph before path finding. edges with less than this frequency will be discarded if they are non-cutting, at a fork, or the end of a path
        

## assembly_min_exact_match_to_remap

**type**: `#!python int`

**environment variable**: `MAVIS_ASSEMBLY_MIN_EXACT_MATCH_TO_REMAP`

**default**: `#!python 15`

The minimum length of exact matches to initiate remapping a read to a contig
        

## assembly_min_remap_coverage

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_ASSEMBLY_MIN_REMAP_COVERAGE`

**default**: `#!python 0.9`

Minimum fraction of the contig sequence which the remapped sequences must align over
        

## assembly_min_remapped_seq

**type**: `#!python int`

**environment variable**: `MAVIS_ASSEMBLY_MIN_REMAPPED_SEQ`

**default**: `#!python 3`

The minimum input sequences that must remap for an assembled contig to be used
        

## assembly_min_uniq

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_ASSEMBLY_MIN_UNIQ`

**default**: `#!python 0.1`

Minimum percent uniq required to keep separate assembled contigs. if contigs are more similar then the lower scoring, then shorter, contig is dropped
        

## assembly_strand_concordance

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_ASSEMBLY_STRAND_CONCORDANCE`

**default**: `#!python 0.51`

When the number of remapped reads from each strand are compared, the ratio must be above this number to decide on the strand
        

## blat_limit_top_aln

**type**: `#!python int`

**environment variable**: `MAVIS_BLAT_LIMIT_TOP_ALN`

**default**: `#!python 10`

Number of results to return from blat (ranking based on score)
        

## blat_min_identity

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_BLAT_MIN_IDENTITY`

**default**: `#!python 0.9`

The minimum percent identity match required for blat results when aligning contigs
        

## breakpoint_color

**type**: `#!python str`

**environment variable**: `MAVIS_BREAKPOINT_COLOR`

**default**: `#!python '#000000'`

Breakpoint outline color
        

## call_error

**type**: `#!python int`

**environment variable**: `MAVIS_CALL_ERROR`

**default**: `#!python 10`

Buffer zone for the evidence window
        

## clean_aligner_files

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_CLEAN_ALIGNER_FILES`

**default**: `#!python False`

Remove the aligner output files after the validation stage is complete. not required for subsequent steps but can be useful in debugging and deep investigation of events
        

## cluster_initial_size_limit

**type**: `#!python int`

**environment variable**: `MAVIS_CLUSTER_INITIAL_SIZE_LIMIT`

**default**: `#!python 25`

The maximum cumulative size of both breakpoints for breakpoint pairs to be used in the initial clustering phase (combining based on overlap)
        

## cluster_radius

**type**: `#!python int`

**environment variable**: `MAVIS_CLUSTER_RADIUS`

**default**: `#!python 100`

Maximum distance allowed between paired breakpoint pairs
        

## concurrency_limit

**type**: `#!python int`

**environment variable**: `MAVIS_CONCURRENCY_LIMIT`

**default**: `#!python None`

The concurrency limit for tasks in any given job array or the number of concurrent processes allowed for a local run
        

## contig_aln_max_event_size

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_ALN_MAX_EVENT_SIZE`

**default**: `#!python 50`

Relates to determining breakpoints when pairing contig alignments. for any given read in a putative pair the soft clipping is extended to include any events of greater than this size. the softclipping is added to the side of the alignment as indicated by the breakpoint we are assigning pairs to
        

## contig_aln_merge_inner_anchor

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_ALN_MERGE_INNER_ANCHOR`

**default**: `#!python 20`

The minimum number of consecutive exact match base pairs to not merge events within a contig alignment
        

## contig_aln_merge_outer_anchor

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_ALN_MERGE_OUTER_ANCHOR`

**default**: `#!python 15`

Minimum consecutively aligned exact matches to anchor an end for merging internal events
        

## contig_aln_min_anchor_size

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_ALN_MIN_ANCHOR_SIZE`

**default**: `#!python 50`

The minimum number of aligned bases for a contig (m or =) in order to simplify. do not have to be consecutive
        

## contig_aln_min_extend_overlap

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_ALN_MIN_EXTEND_OVERLAP`

**default**: `#!python 10`

Minimum number of bases the query coverage interval must be extended by in order to pair alignments as a single split alignment
        

## contig_aln_min_query_consumption

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_CONTIG_ALN_MIN_QUERY_CONSUMPTION`

**default**: `#!python 0.9`

Minimum fraction of the original query sequence that must be used by the read(s) of the alignment
        

## contig_aln_min_score

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_CONTIG_ALN_MIN_SCORE`

**default**: `#!python 0.9`

Minimum score for a contig to be used as evidence in a call by contig
        

## contig_call_distance

**type**: `#!python int`

**environment variable**: `MAVIS_CONTIG_CALL_DISTANCE`

**default**: `#!python 10`

The maximum distance allowed between breakpoint pairs (called by contig) in order for them to pair
        

## dgv_annotation

**type**: `#!python filepath`

**environment variable**: `MAVIS_DGV_ANNOTATION`

**default**: `#!python []`

Path to the dgv reference processed to look like the cytoband file
        

## domain_color

**type**: `#!python str`

**environment variable**: `MAVIS_DOMAIN_COLOR`

**default**: `#!python '#ccccb3'`

Domain fill color
        

## domain_mismatch_color

**type**: `#!python str`

**environment variable**: `MAVIS_DOMAIN_MISMATCH_COLOR`

**default**: `#!python '#b2182b'`

Domain fill color on 0%% match
        

## domain_name_regex_filter

**type**: `#!python str`

**environment variable**: `MAVIS_DOMAIN_NAME_REGEX_FILTER`

**default**: `#!python '^PF\\d+$'`

The regular expression used to select domains to be displayed (filtered by name)
        

## domain_scaffold_color

**type**: `#!python str`

**environment variable**: `MAVIS_DOMAIN_SCAFFOLD_COLOR`

**default**: `#!python '#000000'`

The color of the domain scaffold
        

## draw_fusions_only

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_DRAW_FUSIONS_ONLY`

**default**: `#!python True`

Flag to indicate if events which do not produce a fusion transcript should produce illustrations
        

## draw_non_synonymous_cdna_only

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_DRAW_NON_SYNONYMOUS_CDNA_ONLY`

**default**: `#!python True`

Flag to indicate if events which are synonymous at the cdna level should produce illustrations
        

## drawing_width_iter_increase

**type**: `#!python int`

**environment variable**: `MAVIS_DRAWING_WIDTH_ITER_INCREASE`

**default**: `#!python 500`

The amount (in  pixels) by which to increase the drawing width upon failure to fit
        

## exon_min_focus_size

**type**: `#!python int`

**environment variable**: `MAVIS_EXON_MIN_FOCUS_SIZE`

**default**: `#!python 10`

Minimum size of an exon for it to be granted a label or min exon width
        

## fetch_min_bin_size

**type**: `#!python int`

**environment variable**: `MAVIS_FETCH_MIN_BIN_SIZE`

**default**: `#!python 50`

The minimum size of any bin for reading from a bam file. increasing this number will result in smaller bins being merged or less bins being created (depending on the fetch method)
        

## fetch_reads_bins

**type**: `#!python int`

**environment variable**: `MAVIS_FETCH_READS_BINS`

**default**: `#!python 5`

Number of bins to split an evidence window into to ensure more even sampling of high coverage regions
        

## fetch_reads_limit

**type**: `#!python int`

**environment variable**: `MAVIS_FETCH_READS_LIMIT`

**default**: `#!python 3000`

Maximum number of reads, cap, to loop over for any given evidence window
        

## filter_cdna_synon

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_FILTER_CDNA_SYNON`

**default**: `#!python True`

Filter all annotations synonymous at the cdna level
        

## filter_min_complexity

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_FILTER_MIN_COMPLEXITY`

**default**: `#!python 0.2`

Filter event calls based on call sequence complexity
        

## filter_min_flanking_reads

**type**: `#!python int`

**environment variable**: `MAVIS_FILTER_MIN_FLANKING_READS`

**default**: `#!python 10`

Minimum number of flanking pairs for a call by flanking pairs
        

## filter_min_linking_split_reads

**type**: `#!python int`

**environment variable**: `MAVIS_FILTER_MIN_LINKING_SPLIT_READS`

**default**: `#!python 1`

Minimum number of linking split reads for a call by split reads
        

## filter_min_remapped_reads

**type**: `#!python int`

**environment variable**: `MAVIS_FILTER_MIN_REMAPPED_READS`

**default**: `#!python 5`

Minimum number of remapped reads for a call by contig
        

## filter_min_spanning_reads

**type**: `#!python int`

**environment variable**: `MAVIS_FILTER_MIN_SPANNING_READS`

**default**: `#!python 5`

Minimum number of spanning reads for a call by spanning reads
        

## filter_min_split_reads

**type**: `#!python int`

**environment variable**: `MAVIS_FILTER_MIN_SPLIT_READS`

**default**: `#!python 5`

Minimum number of split reads for a call by split reads
        

## filter_protein_synon

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_FILTER_PROTEIN_SYNON`

**default**: `#!python False`

Filter all annotations synonymous at the protein level
        

## filter_secondary_alignments

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_FILTER_SECONDARY_ALIGNMENTS`

**default**: `#!python True`

Filter secondary alignments when gathering read evidence
        

## filter_trans_homopolymers

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_FILTER_TRANS_HOMOPOLYMERS`

**default**: `#!python True`

Filter all single bp ins/del/dup events that are in a homopolymer region of at least 3 bps and are not paired to a genomic event
        

## flanking_call_distance

**type**: `#!python int`

**environment variable**: `MAVIS_FLANKING_CALL_DISTANCE`

**default**: `#!python 50`

The maximum distance allowed between breakpoint pairs (called by flanking pairs) in order for them to pair
        

## fuzzy_mismatch_number

**type**: `#!python int`

**environment variable**: `MAVIS_FUZZY_MISMATCH_NUMBER`

**default**: `#!python 1`

The number of events/mismatches allowed to be considered a fuzzy match
        

## gene1_color

**type**: `#!python str`

**environment variable**: `MAVIS_GENE1_COLOR`

**default**: `#!python '#657e91'`

The color of genes near the first gene
        

## gene1_color_selected

**type**: `#!python str`

**environment variable**: `MAVIS_GENE1_COLOR_SELECTED`

**default**: `#!python '#518dc5'`

The color of the first gene
        

## gene2_color

**type**: `#!python str`

**environment variable**: `MAVIS_GENE2_COLOR`

**default**: `#!python '#325556'`

The color of genes near the second gene
        

## gene2_color_selected

**type**: `#!python str`

**environment variable**: `MAVIS_GENE2_COLOR_SELECTED`

**default**: `#!python '#4c9677'`

The color of the second gene
        

## import_env

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_IMPORT_ENV`

**default**: `#!python True`

Flag to import environment variables
        

## input_call_distance

**type**: `#!python int`

**environment variable**: `MAVIS_INPUT_CALL_DISTANCE`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by input tools, not validated) in order for them to pair
        

## label_color

**type**: `#!python str`

**environment variable**: `MAVIS_LABEL_COLOR`

**default**: `#!python '#000000'`

The label color
        

## limit_to_chr

**type**: `#!python str`

**environment variable**: `MAVIS_LIMIT_TO_CHR`

**default**: `#!python ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']`

A list of chromosome names to use. breakpointpairs on other chromosomes will be filteredout. for example '1 2 3 4' would filter out events/breakpoint pairs on any chromosomes but 1, 2, 3, and 4
        

## mail_type

**type**: `#!python mavis.schedule.constants.MAIL_TYPE`

**environment variable**: `MAVIS_MAIL_TYPE`

**default**: `#!python 'NONE'`

**accepted values**: `'BEGIN'`, `'END'`, `'FAIL'`, `'ALL'`, `'NONE'`


When to notify the mail_user (if given)
        

## mail_user

**type**: `#!python str`

**environment variable**: `MAVIS_MAIL_USER`

**default**: `#!python ''`

User(s) to send notifications to
        

## mask_fill

**type**: `#!python str`

**environment variable**: `MAVIS_MASK_FILL`

**default**: `#!python '#ffffff'`

Color of mask (for deleted region etc.)
        

## mask_opacity

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_MASK_OPACITY`

**default**: `#!python 0.7`

Opacity of the mask layer
        

## masking

**type**: `#!python filepath`

**environment variable**: `MAVIS_MASKING`

**default**: `#!python []`

File containing regions for which input events overlapping them are dropped prior to validation
        

## max_drawing_retries

**type**: `#!python int`

**environment variable**: `MAVIS_MAX_DRAWING_RETRIES`

**default**: `#!python 5`

The maximum number of retries for attempting a drawing. each iteration the width is extended. if it is still insufficient after this number a gene-level only drawing will be output
        

## max_files

**type**: `#!python int`

**environment variable**: `MAVIS_MAX_FILES`

**default**: `#!python 200`

The maximum number of files to output from clustering/splitting
        

## max_orf_cap

**type**: `#!python int`

**environment variable**: `MAVIS_MAX_ORF_CAP`

**default**: `#!python 3`

The maximum number of orfs to return (best putative orfs will be retained)
        

## max_proximity

**type**: `#!python int`

**environment variable**: `MAVIS_MAX_PROXIMITY`

**default**: `#!python 5000`

The maximum distance away from an annotation before the region in considered to be uninformative
        

## max_sc_preceeding_anchor

**type**: `#!python int`

**environment variable**: `MAVIS_MAX_SC_PRECEEDING_ANCHOR`

**default**: `#!python 6`

When remapping a softclipped read this determines the amount of softclipping allowed on the side opposite of where we expect it. for example for a softclipped read on a breakpoint with a left orientation this limits the amount of softclipping that is allowed on the right. if this is set to none then there is no limit on softclipping
        

## memory_limit

**type**: `#!python int`

**environment variable**: `MAVIS_MEMORY_LIMIT`

**default**: `#!python 16000`

The maximum number of megabytes (mb) any given job is allowed
        

## min_anchor_exact

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_ANCHOR_EXACT`

**default**: `#!python 6`

Applies to re-aligning softclipped reads to the opposing breakpoint. the minimum number of consecutive exact matches to anchor a read to initiate targeted realignment
        

## min_anchor_fuzzy

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_ANCHOR_FUZZY`

**default**: `#!python 10`

Applies to re-aligning softclipped reads to the opposing breakpoint. the minimum length of a fuzzy match to anchor a read to initiate targeted realignment
        

## min_anchor_match

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_MIN_ANCHOR_MATCH`

**default**: `#!python 0.9`

Minimum percent match for a read to be kept as evidence
        

## min_call_complexity

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_MIN_CALL_COMPLEXITY`

**default**: `#!python 0.1`

The minimum complexity score for a call sequence. is an average for non-contig calls. filters low complexity contigs before alignment. see [contig_complexity](#contig_complexity)
        

## min_clusters_per_file

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_CLUSTERS_PER_FILE`

**default**: `#!python 50`

The minimum number of breakpoint pairs to output to a file
        

## min_domain_mapping_match

**type**: `#!python float_fraction`

**environment variable**: `MAVIS_MIN_DOMAIN_MAPPING_MATCH`

**default**: `#!python 0.9`

A number between 0 and 1 representing the minimum percent match a domain must map to the fusion transcript to be displayed
        

## min_double_aligned_to_estimate_insertion_size

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_DOUBLE_ALIGNED_TO_ESTIMATE_INSERTION_SIZE`

**default**: `#!python 2`

The minimum number of reads which map soft-clipped to both breakpoints to assume the size of the untemplated sequence between the breakpoints is at most the read length - 2 * min_softclipping
        

## min_flanking_pairs_resolution

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_FLANKING_PAIRS_RESOLUTION`

**default**: `#!python 10`

The minimum number of flanking reads required to call a breakpoint by flanking evidence
        

## min_linking_split_reads

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_LINKING_SPLIT_READS`

**default**: `#!python 2`

The minimum number of split reads which aligned to both breakpoints
        

## min_mapping_quality

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_MAPPING_QUALITY`

**default**: `#!python 5`

The minimum mapping quality of reads to be used as evidence
        

## min_non_target_aligned_split_reads

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_NON_TARGET_ALIGNED_SPLIT_READS`

**default**: `#!python 1`

The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local alignment to the target region to call a breakpoint by split read evidence
        

## min_orf_size

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_ORF_SIZE`

**default**: `#!python 300`

The minimum length (in base pairs) to retain a putative open reading frame (orf)
        

## min_sample_size_to_apply_percentage

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_SAMPLE_SIZE_TO_APPLY_PERCENTAGE`

**default**: `#!python 10`

Minimum number of aligned bases to compute a match percent. if there are less than this number of aligned bases (match or mismatch) the percent comparator is not used
        

## min_softclipping

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_SOFTCLIPPING`

**default**: `#!python 6`

Minimum number of soft-clipped bases required for a read to be used as soft-clipped evidence
        

## min_spanning_reads_resolution

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_SPANNING_READS_RESOLUTION`

**default**: `#!python 5`

Minimum number of spanning reads required to call an event by spanning evidence
        

## min_splits_reads_resolution

**type**: `#!python int`

**environment variable**: `MAVIS_MIN_SPLITS_READS_RESOLUTION`

**default**: `#!python 3`

Minimum number of split reads required to call a breakpoint by split reads
        

## novel_exon_color

**type**: `#!python str`

**environment variable**: `MAVIS_NOVEL_EXON_COLOR`

**default**: `#!python '#5D3F6A'`

Novel exon fill color
        

## outer_window_min_event_size

**type**: `#!python int`

**environment variable**: `MAVIS_OUTER_WINDOW_MIN_EVENT_SIZE`

**default**: `#!python 125`

The minimum size of an event in order for flanking read evidence to be collected
        

## queue

**type**: `#!python str`

**environment variable**: `MAVIS_QUEUE`

**default**: `#!python ''`

The queue jobs are to be submitted to
        

## reference_genome

**type**: `#!python filepath`

**environment variable**: `MAVIS_REFERENCE_GENOME`

**default**: `#!python []`

Path to the human reference genome fasta file
        

## remote_head_ssh

**type**: `#!python str`

**environment variable**: `MAVIS_REMOTE_HEAD_SSH`

**default**: `#!python ''`

Ssh target for remote scheduler commands
        

## scaffold_color

**type**: `#!python str`

**environment variable**: `MAVIS_SCAFFOLD_COLOR`

**default**: `#!python '#000000'`

The color used for the gene/transcripts scaffolds
        

## scheduler

**type**: `#!python mavis.schedule.constants.SCHEDULER`

**environment variable**: `MAVIS_SCHEDULER`

**default**: `#!python 'SLURM'`

**accepted values**: `'SGE'`, `'SLURM'`, `'TORQUE'`, `'LOCAL'`


The scheduler being used
        

## spanning_call_distance

**type**: `#!python int`

**environment variable**: `MAVIS_SPANNING_CALL_DISTANCE`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by spanning reads) in order for them to pair
        

## splice_color

**type**: `#!python str`

**environment variable**: `MAVIS_SPLICE_COLOR`

**default**: `#!python '#000000'`

Splicing lines color
        

## split_call_distance

**type**: `#!python int`

**environment variable**: `MAVIS_SPLIT_CALL_DISTANCE`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by split reads) in order for them to pair
        

## stdev_count_abnormal

**type**: `#!python float`

**environment variable**: `MAVIS_STDEV_COUNT_ABNORMAL`

**default**: `#!python 3.0`

The number of standard deviations away from the normal considered expected and therefore not qualifying as flanking reads
        

## strand_determining_read

**type**: `#!python int`

**environment variable**: `MAVIS_STRAND_DETERMINING_READ`

**default**: `#!python 2`

1 or 2. the read in the pair which determines if (assuming a stranded protocol) the first or second read in the pair matches the strand sequenced
        

## template_metadata

**type**: `#!python filepath`

**environment variable**: `MAVIS_TEMPLATE_METADATA`

**default**: `#!python []`

File containing the cytoband template information. used for illustrations only
        

## time_limit

**type**: `#!python int`

**environment variable**: `MAVIS_TIME_LIMIT`

**default**: `#!python 57600`

The time in seconds any given jobs is allowed
        

## trans_fetch_reads_limit

**type**: `#!python int`

**environment variable**: `MAVIS_TRANS_FETCH_READS_LIMIT`

**default**: `#!python 12000`

Related to [fetch_reads_limit](#fetch_reads_limit). overrides fetch_reads_limit for transcriptome libraries when set. if this has a value of none then fetch_reads_limit will be used for transcriptome libraries instead
        

## trans_min_mapping_quality

**type**: `#!python int`

**environment variable**: `MAVIS_TRANS_MIN_MAPPING_QUALITY`

**default**: `#!python 0`

Related to [min_mapping_quality](#min_mapping_quality). overrides the min_mapping_quality if the library is a transcriptome and this is set to any number not none. if this value is none, min_mapping_quality is used for transcriptomes aswell as genomes
        

## trans_validation_memory

**type**: `#!python int`

**environment variable**: `MAVIS_TRANS_VALIDATION_MEMORY`

**default**: `#!python 18000`

Default memory limit (mb) for the validation stage (for transcriptomes)
        

## uninformative_filter

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_UNINFORMATIVE_FILTER`

**default**: `#!python False`

Flag that determines if breakpoint pairs which are not within max_proximity to any annotations are filtered out prior to clustering
        

## validation_memory

**type**: `#!python int`

**environment variable**: `MAVIS_VALIDATION_MEMORY`

**default**: `#!python 16000`

Default memory limit (mb) for the validation stage
        

## width

**type**: `#!python int`

**environment variable**: `MAVIS_WIDTH`

**default**: `#!python 1000`

The drawing width in pixels
        

## write_evidence_files

**type**: `#!python cast_boolean`

**environment variable**: `MAVIS_WRITE_EVIDENCE_FILES`

**default**: `#!python True`

Write the intermediate bam and bed files containing the raw evidence collected and contigs aligned. not required for subsequent steps but can be useful in debugging and deep investigation of events
        

