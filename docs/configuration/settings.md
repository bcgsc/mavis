

# Configurable Settings
## annotate.annotation_filters

**type**: `#!python List[str]`

**default**: `#!python ['choose_more_annotated', 'choose_transcripts_by_priority']`

A comma separated list of filters to apply to putative annotations

**schema definition**:
```json
{
    "items": {
        "enum": [
            "choose_more_annotated",
            "choose_transcripts_by_priority"
        ],
        "type": "string"
    },
    "type": "array"
}
```


## annotate.draw_fusions_only

**type**: `#!python boolean`

**default**: `#!python True`

Flag to indicate if events which do not produce a fusion transcript should produce illustrations

**schema definition**:
```json
{
    "type": "boolean"
}
```


## annotate.draw_non_synonymous_cdna_only

**type**: `#!python boolean`

**default**: `#!python True`

Flag to indicate if events which are synonymous at the cdna level should produce illustrations

**schema definition**:
```json
{
    "type": "boolean"
}
```


## annotate.max_orf_cap

**type**: `#!python int`

**default**: `#!python 3`

The maximum number of orfs to return (best putative orfs will be retained)

**schema definition**:
```json
{
    "type": "integer"
}
```


## annotate.min_domain_mapping_match

**type**: `#!python number`

**default**: `#!python 0.9`

A number between 0 and 1 representing the minimum percent match a domain must map to the fusion transcript to be displayed

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## annotate.min_orf_size

**type**: `#!python int`

**default**: `#!python 300`

The minimum length (in base pairs) to retain a putative open reading frame (orf)

**schema definition**:
```json
{
    "type": "integer"
}
```


## bam_stats.distribution_fraction

**type**: `#!python number`

**default**: `#!python 0.97`

the proportion of the distribution to use in computing stdev

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0.01,
    "type": "number"
}
```


## bam_stats.sample_bin_size

**type**: `#!python int`

**default**: `#!python 1000`

how large to make the sample bin (in bp)

**schema definition**:
```json
{
    "type": "integer"
}
```


## bam_stats.sample_cap

**type**: `#!python int`

**default**: `#!python 1000`

maximum number of reads to collect for any given sample region

**schema definition**:
```json
{
    "type": "integer"
}
```


## bam_stats.sample_size

**type**: `#!python int`

**default**: `#!python 500`

the number of genes/bins to compute stats over

**schema definition**:
```json
{
    "type": "integer"
}
```


## cluster.cluster_initial_size_limit

**type**: `#!python int`

**default**: `#!python 25`

The maximum cumulative size of both breakpoints for breakpoint pairs to be used in the initial clustering phase (combining based on overlap)

**schema definition**:
```json
{
    "type": "integer"
}
```


## cluster.cluster_radius

**type**: `#!python int`

**default**: `#!python 100`

Maximum distance allowed between paired breakpoint pairs

**schema definition**:
```json
{
    "type": "integer"
}
```


## cluster.limit_to_chr

**type**: `#!python Union[List, null]`

**default**: `#!python ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']`

A list of chromosome names to use. breakpointpairs on other chromosomes will be filteredout. for example '1 2 3 4' would filter out events/breakpoint pairs on any chromosomes but 1, 2, 3, and 4

**schema definition**:
```json
{
    "items": {
        "type": "string"
    },
    "type": [
        "array",
        "null"
    ]
}
```


## cluster.max_files

**type**: `#!python int`

**default**: `#!python 200`

The maximum number of files to output from clustering/splitting

**schema definition**:
```json
{
    "minimum": 1,
    "type": "integer"
}
```


## cluster.max_proximity

**type**: `#!python int`

**default**: `#!python 5000`

The maximum distance away from an annotation before the region in considered to be uninformative

**schema definition**:
```json
{
    "type": "integer"
}
```


## cluster.min_clusters_per_file

**type**: `#!python int`

**default**: `#!python 50`

The minimum number of breakpoint pairs to output to a file

**schema definition**:
```json
{
    "minimum": 1,
    "type": "integer"
}
```


## cluster.split_only

**type**: `#!python boolean`

**default**: `#!python False`

just split the input files, do not merge input breakpoints into clusters

**schema definition**:
```json
{
    "type": "boolean"
}
```


## cluster.uninformative_filter

**type**: `#!python boolean`

**default**: `#!python False`

Flag that determines if breakpoint pairs which are not within max_proximity to any annotations are filtered out prior to clustering

**schema definition**:
```json
{
    "type": "boolean"
}
```


## illustrate.breakpoint_color

**type**: `#!python str`

**default**: `#!python '#000000'`

Breakpoint outline color

**schema definition**:
```json
{
    "type": "string"
}
```


## illustrate.domain_color

**type**: `#!python str`

**default**: `#!python '#ccccb3'`

Domain fill color

**schema definition**:
```json
{
    "type": "string"
}
```


## illustrate.domain_mismatch_color

**type**: `#!python str`

**default**: `#!python '#b2182b'`

Domain fill color on 0%% match

**schema definition**:
```json
{
    "type": "string"
}
```


## illustrate.domain_name_regex_filter

**type**: `#!python str`

**default**: `#!python '^PF\\d+$'`

The regular expression used to select domains to be displayed (filtered by name)

**schema definition**:
```json
{
    "type": "string"
}
```


## illustrate.domain_scaffold_color

**type**: `#!python str`

**default**: `#!python '#000000'`

The color of the domain scaffold

**schema definition**:
```json
{
    "type": "string"
}
```


## illustrate.drawing_width_iter_increase

**type**: `#!python int`

**default**: `#!python 500`

The amount (in  pixels) by which to increase the drawing width upon failure to fit

**schema definition**:
```json
{
    "type": "integer"
}
```


## illustrate.exon_min_focus_size

**type**: `#!python int`

**default**: `#!python 10`

Minimum size of an exon for it to be granted a label or min exon width

**schema definition**:
```json
{
    "type": "integer"
}
```


## illustrate.gene1_color

**type**: `#!python str`

**default**: `#!python '#657e91'`

The color of genes near the first gene

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.gene1_color_selected

**type**: `#!python str`

**default**: `#!python '#518dc5'`

The color of the first gene

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.gene2_color

**type**: `#!python str`

**default**: `#!python '#325556'`

The color of genes near the second gene

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.gene2_color_selected

**type**: `#!python str`

**default**: `#!python '#4c9677'`

The color of the second gene

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.label_color

**type**: `#!python str`

**default**: `#!python '#000000'`

The label color

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.mask_fill

**type**: `#!python str`

**default**: `#!python '#ffffff'`

Color of mask (for deleted region etc.)

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.mask_opacity

**type**: `#!python number`

**default**: `#!python 0.7`

Opacity of the mask layer

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## illustrate.max_drawing_retries

**type**: `#!python int`

**default**: `#!python 5`

The maximum number of retries for attempting a drawing. each iteration the width is extended. if it is still insufficient after this number a gene-level only drawing will be output

**schema definition**:
```json
{
    "type": "integer"
}
```


## illustrate.novel_exon_color

**type**: `#!python str`

**default**: `#!python '#5D3F6A'`

Novel exon fill color

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.scaffold_color

**type**: `#!python str`

**default**: `#!python '#000000'`

The color used for the gene/transcripts scaffolds

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.splice_color

**type**: `#!python str`

**default**: `#!python '#000000'`

Splicing lines color

**schema definition**:
```json
{
    "pattern": "^#[a-zA-Z0-9]{6}",
    "type": "string"
}
```


## illustrate.width

**type**: `#!python int`

**default**: `#!python 1000`

The drawing width in pixels

**schema definition**:
```json
{
    "type": "integer"
}
```


## log

**type**: `#!python str`

**default**: `#!python None`



**schema definition**:
```json
{
    "type": "string"
}
```


## log_level

**type**: `#!python str`

**default**: `#!python 'INFO'`



**schema definition**:
```json
{
    "enum": [
        "INFO",
        "DEBUG"
    ],
    "type": "string"
}
```


## output_dir

**type**: `#!python str`

**default**: `#!python None`



**schema definition**:
```json
{
    "type": "string"
}
```


## pairing.contig_call_distance

**type**: `#!python int`

**default**: `#!python 10`

The maximum distance allowed between breakpoint pairs (called by contig) in order for them to pair

**schema definition**:
```json
{
    "type": "integer"
}
```


## pairing.flanking_call_distance

**type**: `#!python int`

**default**: `#!python 50`

The maximum distance allowed between breakpoint pairs (called by flanking pairs) in order for them to pair

**schema definition**:
```json
{
    "type": "integer"
}
```


## pairing.input_call_distance

**type**: `#!python int`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by input tools, not validated) in order for them to pair

**schema definition**:
```json
{
    "type": "integer"
}
```


## pairing.spanning_call_distance

**type**: `#!python int`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by spanning reads) in order for them to pair

**schema definition**:
```json
{
    "type": "integer"
}
```


## pairing.split_call_distance

**type**: `#!python int`

**default**: `#!python 20`

The maximum distance allowed between breakpoint pairs (called by split reads) in order for them to pair

**schema definition**:
```json
{
    "type": "integer"
}
```


## reference.aligner_reference

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        "tests/data/mock_reference_genome.2bit"
    ],
    "items": {
        "type": "string"
    },
    "maxItems": 1,
    "minItems": 1,
    "type": "array"
}
```


## reference.annotations

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        "tests/data/mock_annotations.json"
    ],
    "items": {
        "type": "string"
    },
    "minItems": 1,
    "type": "array"
}
```


## reference.dgv_annotation

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        [
            "tests/data/mock_dgv_annotation.txt"
        ]
    ],
    "items": {
        "type": "string"
    },
    "minItems": 1,
    "type": "array"
}
```


## reference.masking

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        [
            "tests/data/mock_masking.tab"
        ]
    ],
    "items": {
        "type": "string"
    },
    "minItems": 1,
    "type": "array"
}
```


## reference.reference_genome

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        [
            "tests/data/mock_reference_genome.fa"
        ]
    ],
    "items": {
        "type": "string"
    },
    "minItems": 1,
    "type": "array"
}
```


## reference.template_metadata

**type**: `#!python List[str]`

**default**: `#!python None`



**schema definition**:
```json
{
    "examples": [
        [
            "tests/data/cytoBand.txt"
        ]
    ],
    "items": {
        "type": "string"
    },
    "minItems": 1,
    "type": "array"
}
```


## skip_stage.validate

**type**: `#!python boolean`

**default**: `#!python False`

skip the validation stage of the MAVIS pipeline

**schema definition**:
```json
{
    "type": "boolean"
}
```


## summary.filter_cdna_synon

**type**: `#!python boolean`

**default**: `#!python True`

Filter all annotations synonymous at the cdna level

**schema definition**:
```json
{
    "type": "boolean"
}
```


## summary.filter_min_complexity

**type**: `#!python number`

**default**: `#!python 0.2`

Filter event calls based on call sequence complexity

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## summary.filter_min_flanking_reads

**type**: `#!python int`

**default**: `#!python 10`

Minimum number of flanking pairs for a call by flanking pairs

**schema definition**:
```json
{
    "type": "integer"
}
```


## summary.filter_min_linking_split_reads

**type**: `#!python int`

**default**: `#!python 1`

Minimum number of linking split reads for a call by split reads

**schema definition**:
```json
{
    "type": "integer"
}
```


## summary.filter_min_remapped_reads

**type**: `#!python int`

**default**: `#!python 5`

Minimum number of remapped reads for a call by contig

**schema definition**:
```json
{
    "type": "integer"
}
```


## summary.filter_min_spanning_reads

**type**: `#!python int`

**default**: `#!python 5`

Minimum number of spanning reads for a call by spanning reads

**schema definition**:
```json
{
    "type": "integer"
}
```


## summary.filter_min_split_reads

**type**: `#!python int`

**default**: `#!python 5`

Minimum number of split reads for a call by split reads

**schema definition**:
```json
{
    "type": "integer"
}
```


## summary.filter_protein_synon

**type**: `#!python boolean`

**default**: `#!python False`

Filter all annotations synonymous at the protein level

**schema definition**:
```json
{
    "type": "boolean"
}
```


## summary.filter_trans_homopolymers

**type**: `#!python boolean`

**default**: `#!python True`

Filter all single bp ins/del/dup events that are in a homopolymer region of at least 3 bps and are not paired to a genomic event

**schema definition**:
```json
{
    "type": "boolean"
}
```


## validate.aligner

**type**: `#!python str`

**default**: `#!python 'blat'`

The aligner to use to map the contigs/reads back to the reference e.g blat or bwa

**schema definition**:
```json
{
    "enum": [
        "bwa mem",
        "blat"
    ],
    "type": "string"
}
```


## validate.assembly_kmer_size

**type**: `#!python number`

**default**: `#!python 0.74`

The percent of the read length to make kmers for assembly

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.assembly_max_paths

**type**: `#!python int`

**default**: `#!python 8`

The maximum number of paths to resolve. this is used to limit when there is a messy assembly graph to resolve. the assembly will pre-calculate the number of paths (or putative assemblies) and stop if it is greater than the given setting

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.assembly_min_edge_trim_weight

**type**: `#!python int`

**default**: `#!python 3`

This is used to simplify the debruijn graph before path finding. edges with less than this frequency will be discarded if they are non-cutting, at a fork, or the end of a path

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.assembly_min_exact_match_to_remap

**type**: `#!python int`

**default**: `#!python 15`

The minimum length of exact matches to initiate remapping a read to a contig

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.assembly_min_remap_coverage

**type**: `#!python number`

**default**: `#!python 0.9`

Minimum fraction of the contig sequence which the remapped sequences must align over

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.assembly_min_remapped_seq

**type**: `#!python int`

**default**: `#!python 3`

The minimum input sequences that must remap for an assembled contig to be used

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.assembly_min_uniq

**type**: `#!python number`

**default**: `#!python 0.1`

Minimum percent uniq required to keep separate assembled contigs. if contigs are more similar then the lower scoring, then shorter, contig is dropped

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.assembly_strand_concordance

**type**: `#!python number`

**default**: `#!python 0.51`

When the number of remapped reads from each strand are compared, the ratio must be above this number to decide on the strand

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.blat_limit_top_aln

**type**: `#!python int`

**default**: `#!python 10`

Number of results to return from blat (ranking based on score)

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.blat_min_identity

**type**: `#!python number`

**default**: `#!python 0.9`

The minimum percent identity match required for blat results when aligning contigs

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.call_error

**type**: `#!python int`

**default**: `#!python 10`

Buffer zone for the evidence window

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.clean_aligner_files

**type**: `#!python boolean`

**default**: `#!python False`

Remove the aligner output files after the validation stage is complete. not required for subsequent steps but can be useful in debugging and deep investigation of events

**schema definition**:
```json
{
    "type": "boolean"
}
```


## validate.contig_aln_max_event_size

**type**: `#!python int`

**default**: `#!python 50`

Relates to determining breakpoints when pairing contig alignments. for any given read in a putative pair the soft clipping is extended to include any events of greater than this size. the softclipping is added to the side of the alignment as indicated by the breakpoint we are assigning pairs to

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.contig_aln_merge_inner_anchor

**type**: `#!python int`

**default**: `#!python 20`

The minimum number of consecutive exact match base pairs to not merge events within a contig alignment

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.contig_aln_merge_outer_anchor

**type**: `#!python int`

**default**: `#!python 15`

Minimum consecutively aligned exact matches to anchor an end for merging internal events

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.contig_aln_min_anchor_size

**type**: `#!python int`

**default**: `#!python 50`

The minimum number of aligned bases for a contig (m or =) in order to simplify. do not have to be consecutive

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.contig_aln_min_extend_overlap

**type**: `#!python int`

**default**: `#!python 10`

Minimum number of bases the query coverage interval must be extended by in order to pair alignments as a single split alignment

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.contig_aln_min_query_consumption

**type**: `#!python number`

**default**: `#!python 0.9`

Minimum fraction of the original query sequence that must be used by the read(s) of the alignment

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.contig_aln_min_score

**type**: `#!python number`

**default**: `#!python 0.9`

Minimum score for a contig to be used as evidence in a call by contig

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.fetch_min_bin_size

**type**: `#!python int`

**default**: `#!python 50`

The minimum size of any bin for reading from a bam file. increasing this number will result in smaller bins being merged or less bins being created (depending on the fetch method)

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.fetch_reads_bins

**type**: `#!python int`

**default**: `#!python 5`

Number of bins to split an evidence window into to ensure more even sampling of high coverage regions

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.fetch_reads_limit

**type**: `#!python int`

**default**: `#!python 3000`

Maximum number of reads, cap, to loop over for any given evidence window

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.filter_secondary_alignments

**type**: `#!python boolean`

**default**: `#!python True`

Filter secondary alignments when gathering read evidence

**schema definition**:
```json
{
    "type": "boolean"
}
```


## validate.fuzzy_mismatch_number

**type**: `#!python int`

**default**: `#!python 1`

The number of events/mismatches allowed to be considered a fuzzy match

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.max_sc_preceeding_anchor

**type**: `#!python int`

**default**: `#!python 6`

When remapping a softclipped read this determines the amount of softclipping allowed on the side opposite of where we expect it. for example for a softclipped read on a breakpoint with a left orientation this limits the amount of softclipping that is allowed on the right. if this is set to none then there is no limit on softclipping

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_anchor_exact

**type**: `#!python int`

**default**: `#!python 6`

Applies to re-aligning softclipped reads to the opposing breakpoint. the minimum number of consecutive exact matches to anchor a read to initiate targeted realignment

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_anchor_fuzzy

**type**: `#!python int`

**default**: `#!python 10`

Applies to re-aligning softclipped reads to the opposing breakpoint. the minimum length of a fuzzy match to anchor a read to initiate targeted realignment

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_anchor_match

**type**: `#!python number`

**default**: `#!python 0.9`

Minimum percent match for a read to be kept as evidence

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.min_call_complexity

**type**: `#!python number`

**default**: `#!python 0.1`

The minimum complexity score for a call sequence. is an average for non-contig calls. filters low complexity contigs before alignment. see [contig_complexity](#contig_complexity)

**schema definition**:
```json
{
    "maximum": 1,
    "minimum": 0,
    "type": "number"
}
```


## validate.min_double_aligned_to_estimate_insertion_size

**type**: `#!python int`

**default**: `#!python 2`

The minimum number of reads which map soft-clipped to both breakpoints to assume the size of the untemplated sequence between the breakpoints is at most the read length - 2 * min_softclipping

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_flanking_pairs_resolution

**type**: `#!python int`

**default**: `#!python 10`

The minimum number of flanking reads required to call a breakpoint by flanking evidence

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_linking_split_reads

**type**: `#!python int`

**default**: `#!python 2`

The minimum number of split reads which aligned to both breakpoints

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_mapping_quality

**type**: `#!python int`

**default**: `#!python 5`

The minimum mapping quality of reads to be used as evidence

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_non_target_aligned_split_reads

**type**: `#!python int`

**default**: `#!python 1`

The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local alignment to the target region to call a breakpoint by split read evidence

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_sample_size_to_apply_percentage

**type**: `#!python int`

**default**: `#!python 10`

Minimum number of aligned bases to compute a match percent. if there are less than this number of aligned bases (match or mismatch) the percent comparator is not used

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_softclipping

**type**: `#!python int`

**default**: `#!python 6`

Minimum number of soft-clipped bases required for a read to be used as soft-clipped evidence

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_spanning_reads_resolution

**type**: `#!python int`

**default**: `#!python 5`

Minimum number of spanning reads required to call an event by spanning evidence

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.min_splits_reads_resolution

**type**: `#!python int`

**default**: `#!python 3`

Minimum number of split reads required to call a breakpoint by split reads

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.outer_window_min_event_size

**type**: `#!python int`

**default**: `#!python 125`

The minimum size of an event in order for flanking read evidence to be collected

**schema definition**:
```json
{
    "type": "integer"
}
```


## validate.stdev_count_abnormal

**type**: `#!python number`

**default**: `#!python 3`

The number of standard deviations away from the normal considered expected and therefore not qualifying as flanking reads

**schema definition**:
```json
{
    "type": "number"
}
```


## validate.trans_fetch_reads_limit

**type**: `#!python Union[int, null]`

**default**: `#!python 12000`

Related to [fetch_reads_limit](#fetch_reads_limit). overrides fetch_reads_limit for transcriptome libraries when set. if this has a value of none then fetch_reads_limit will be used for transcriptome libraries instead

**schema definition**:
```json
{
    "type": [
        "integer",
        "null"
    ]
}
```


## validate.trans_min_mapping_quality

**type**: `#!python Union[int, null]`

**default**: `#!python 0`

Related to [min_mapping_quality](#min_mapping_quality). overrides the min_mapping_quality if the library is a transcriptome and this is set to any number not none. if this value is none, min_mapping_quality is used for transcriptomes aswell as genomes

**schema definition**:
```json
{
    "type": [
        "integer",
        "null"
    ]
}
```


## validate.write_evidence_files

**type**: `#!python boolean`

**default**: `#!python True`

Write the intermediate bam and bed files containing the raw evidence collected and contigs aligned. not required for subsequent steps but can be useful in debugging and deep investigation of events

**schema definition**:
```json
{
    "type": "boolean"
}
```


