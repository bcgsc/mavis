from ..constants import float_fraction
from ..align import SUPPORTED_ALIGNER
from ..util import WeakMavisNamespace

PASS_FILENAME = 'validation-passed.tab'

DEFAULTS = WeakMavisNamespace()
"""
- :term:`aligner`
- :term:`assembly_kmer_size`
- :term:`assembly_max_paths`
- :term:`assembly_min_edge_trim_weight`
- :term:`assembly_min_exact_match_to_remap`
- :term:`assembly_min_remap_coverage`
- :term:`assembly_min_remapped_seq`
- :term:`assembly_min_uniq`
- :term:`assembly_strand_concordance`
- :term:`blat_limit_top_aln`
- :term:`blat_min_identity`
- :term:`call_error`
- :term:`contig_aln_max_event_size`
- :term:`contig_aln_merge_inner_anchor`
- :term:`contig_aln_merge_outer_anchor`
- :term:`contig_aln_min_anchor_size`
- :term:`contig_aln_min_extend_overlap`
- :term:`contig_aln_min_query_consumption`
- :term:`contig_aln_min_score`
- :term:`fetch_min_bin_size`
- :term:`fetch_reads_bins`
- :term:`fetch_reads_limit`
- :term:`filter_secondary_alignments`
- :term:`fuzzy_mismatch_number`
- :term:`max_sc_preceeding_anchor`
- :term:`min_anchor_exact`
- :term:`min_anchor_fuzzy`
- :term:`min_anchor_match`
- :term:`min_double_aligned_to_estimate_insertion_size`
- :term:`min_flanking_pairs_resolution`
- :term:`min_linking_split_reads`
- :term:`min_mapping_quality`
- :term:`min_non_target_aligned_split_reads`
- :term:`min_sample_size_to_apply_percentage`
- :term:`min_softclipping`
- :term:`min_spanning_reads_resolution`
- :term:`min_splits_reads_resolution`
- :term:`outer_window_min_event_size`
- :term:`stdev_count_abnormal`
- :term:`strand_determining_read`

"""
DEFAULTS.add(
    'min_call_complexity', 0.10, cast_type=float_fraction,
    defn='The minimum complexity score for a call sequence. Is an average for non-contig calls. Filters '
         'low complexity contigs before alignment. see :term:`contig_complexity`')
DEFAULTS.add(
    'aligner', SUPPORTED_ALIGNER.BLAT, cast_type=SUPPORTED_ALIGNER,
    defn='the aligner to use to map the contigs/reads back to the reference e.g blat or bwa')
DEFAULTS.add(
    'assembly_kmer_size', 0.74, cast_type=float_fraction,
    defn='The percent of the read length to make kmers for assembly')
DEFAULTS.add(
    'assembly_max_paths', 8,
    defn='the maximum number of paths to resolve. This is used to limit when there is a messy assembly graph to '
    'resolve. The assembly will pre-calculate the number of paths (or putative assemblies) and stop if it is greater '
    'than the given setting.')
DEFAULTS.add(
    'assembly_min_uniq', 0.10, cast_type=float_fraction,
    defn='Minimum percent uniq required to keep separate assembled contigs. If contigs are more similar then the lower scoring, then shorter, contig is dropped')
DEFAULTS.add(
    'assembly_min_exact_match_to_remap', 15,
    defn='The minimum length of exact matches to initiate remapping a read to a contig')
DEFAULTS.add(
    'assembly_min_edge_trim_weight', 3,
    defn='this is used to simplify the DeBruijn graph before path finding. Edges with less than this frequency will '
    'be discarded if they are non-cutting, at a fork, or the end of a path')
DEFAULTS.add(
    'assembly_min_remap_coverage', 0.9, cast_type=float_fraction,
    defn='Minimum fraction of the contig sequence which the remapped sequences must align over')
DEFAULTS.add(
    'assembly_min_remapped_seq', 3,
    defn='The minimum input sequences that must remap for an assembled contig to be used')
DEFAULTS.add(
    'assembly_strand_concordance', 0.51, cast_type=float_fraction,
    defn='When the number of remapped reads from each strand are compared, the ratio must be above this number to '
    'decide on the strand')
DEFAULTS.add(
    'blat_min_identity', 0.9, cast_type=float_fraction,
    defn='The minimum percent identity match required for blat results when aligning contigs')
DEFAULTS.add(
    'blat_limit_top_aln', 10,
    defn='Number of results to return from blat (ranking based on score)')
DEFAULTS.add(
    'call_error', 10,
    defn='buffer zone for the evidence window')
DEFAULTS.add(
    'contig_aln_max_event_size', 50,
    defn='relates to determining breakpoints when pairing contig alignments. For any given read in a putative pair '
    'the soft clipping is extended to include any events of greater than this size. The softclipping is added to the '
    'side of the alignment as indicated by the breakpoint we are assigning pairs to')
DEFAULTS.add(
    'contig_aln_merge_inner_anchor', 20,
    defn='the minimum number of consecutive exact match base pairs to not merge events within a contig alignment')
DEFAULTS.add(
    'contig_aln_merge_outer_anchor', 15,
    defn='minimum consecutively aligned exact matches to anchor an end for merging internal events')
DEFAULTS.add(
    'contig_aln_min_anchor_size', 50,
    defn='the minimum number of aligned bases for a contig (M or =) in order to simplify. Do not have to be consecutive.')
DEFAULTS.add(
    'contig_aln_min_query_consumption', 0.9, cast_type=float_fraction,
    defn='minimum fraction of the original query sequence that must be used by the read(s) of the alignment')
DEFAULTS.add(
    'contig_aln_min_extend_overlap', 10,
    defn='minimum number of bases the query coverage interval must be extended by in order to pair alignments as a single split alignment')
DEFAULTS.add(
    'contig_aln_min_score', 0.9, cast_type=float_fraction, defn='minimum score for a contig to be used as evidence in a call by contig')
DEFAULTS.add(
    'fetch_min_bin_size', 50,
    defn='the minimum size of any bin for reading from a bam file. Increasing this number will result in smaller bins '
    'being merged or less bins being created (depending on the fetch method)')
DEFAULTS.add(
    'fetch_reads_bins', 5,
    defn='number of bins to split an evidence window into to ensure more even sampling of high coverage regions')
DEFAULTS.add(
    'fetch_reads_limit', 3000,
    defn='maximum number of reads, cap, to loop over for any given evidence window')
DEFAULTS.add(
    'trans_fetch_reads_limit', 12000, cast_type=int, nullable=True,
    defn='Related to :term:`fetch_reads_limit`. Overrides fetch_reads_limit for transcriptome libraries when set. '
    'If this has a value of None then fetch_reads_limit will be used for transcriptome libraries instead')
DEFAULTS.add(
    'filter_secondary_alignments', True,
    defn='filter secondary alignments when gathering read evidence')
DEFAULTS.add(
    'fuzzy_mismatch_number', 1,
    defn='The number of events/mismatches allowed to be considered a fuzzy match')
DEFAULTS.add(
    'max_sc_preceeding_anchor', 6,
    defn='when remapping a softclipped read this determines the amount of softclipping allowed on the side opposite of '
    'where we expect it. For example for a softclipped read on a breakpoint with a left orientation this limits the '
    'amount of softclipping that is allowed on the right. If this is set to None then there is no limit on softclipping')
DEFAULTS.add(
    'min_anchor_exact', 6, defn='Applies to re-aligning softclipped reads to the opposing breakpoint. The minimum '
    'number of consecutive exact matches to anchor a read to initiate targeted realignment')
DEFAULTS.add(
    'min_anchor_fuzzy', 10, defn='Applies to re-aligning softclipped reads to the opposing breakpoint. The minimum '
    'length of a fuzzy match to anchor a read to initiate targeted realignment')
DEFAULTS.add(
    'min_anchor_match', 0.9, cast_type=float_fraction,
    defn='Minimum percent match for a read to be kept as evidence')
DEFAULTS.add(
    'min_double_aligned_to_estimate_insertion_size', 2,
    defn='The minimum number of reads which map soft-clipped to both breakpoints to assume the size of the '
    'untemplated sequence between the breakpoints is at most the read length - 2 * min_softclipping')
DEFAULTS.add(
    'min_flanking_pairs_resolution', 10,
    defn='the minimum number of flanking reads required to call a breakpoint by flanking evidence')
DEFAULTS.add(
    'min_linking_split_reads', 2,
    defn='The minimum number of split reads which aligned to both breakpoints')
DEFAULTS.add(
    'min_mapping_quality', 5,
    defn='the minimum mapping quality of reads to be used as evidence')
DEFAULTS.add(
    'trans_min_mapping_quality', 0, cast_type=int, nullable=True,
    defn='Related to :term:`min_mapping_quality`. Overrides the min_mapping_quality if the library is a transcriptome '
    'and this is set to any number not None. If this value is None, min_mapping_quality is used for transcriptomes as'
    'well as genomes')
DEFAULTS.add(
    'min_non_target_aligned_split_reads', 1,
    defn='The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local '
    'alignment to the target region to call a breakpoint by split read evidence')
DEFAULTS.add(
    'min_sample_size_to_apply_percentage', 10, defn='Minimum number of aligned bases to compute a match percent. '
    'If there are less than this number of aligned bases (match or mismatch) the percent comparator is not used')
DEFAULTS.add(
    'min_softclipping', 6,
    defn='minimum number of soft-clipped bases required for a read to be used as soft-clipped evidence')
DEFAULTS.add(
    'min_spanning_reads_resolution', 5, defn='Minimum number of spanning reads required to call an event by spanning evidence')
DEFAULTS.add(
    'min_splits_reads_resolution', 3,
    defn='minimum number of split reads required to call a breakpoint by split reads')
DEFAULTS.add(
    'stdev_count_abnormal', 3.0,
    defn='the number of standard deviations away from the normal considered expected and therefore not qualifying as '
    'flanking reads')
DEFAULTS.add(
    'strand_determining_read', 2,
    defn='1 or 2. The read in the pair which determines if (assuming a stranded protocol) the first or second read in '
    'the pair matches the strand sequenced')
DEFAULTS.add(
    'outer_window_min_event_size', 125,
    defn='the minimum size of an event in order for flanking read evidence to be collected')
DEFAULTS.add(
    'write_evidence_files', True, defn='write the intermediate bam and bed files containing the raw evidence collected and '
    'contigs aligned. Not required for subsequent steps but can be useful in debugging and deep investigation of events')
DEFAULTS.add(
    'clean_aligner_files', False, defn='Remove the aligner output files after the validation stage is complete. Not'
    ' required for subsequent steps but can be useful in debugging and deep investigation of events')
