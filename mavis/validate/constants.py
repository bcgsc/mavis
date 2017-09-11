from ..util import MavisNamespace

DEFAULTS = MavisNamespace(
    aligner='blat',
    assembly_include_flanking_pairs=True,
    assembly_include_half_mapped_reads=True,
    assembly_max_kmer_size=None,
    assembly_max_kmer_strict=True,
    assembly_max_paths=4,
    assembly_min_edge_weight=2,
    assembly_min_exact_match_to_remap=4,
    assembly_min_nc_edge_weight=4,
    assembly_min_remap_coverage=0.9,
    assembly_min_remapped_seq=3,
    assembly_min_tgt_to_exclude_half_map=7,
    assembly_strand_concordance=0.51,
    blat_min_identity=0.9,
    blat_limit_top_aln=10,
    call_error=10,
    consensus_req=3,
    contig_aln_max_event_size=50,
    contig_aln_merge_inner_anchor=20,
    contig_aln_merge_outer_anchor=15,
    contig_aln_min_anchor_size=50,
    contig_aln_min_query_consumption=0.9,
    fetch_reads_bins=5,
    fetch_reads_limit=10000,
    fetch_method_individual=True,
    fetch_min_bin_size=50,
    filter_secondary_alignments=True,
    fuzzy_mismatch_number=1,
    max_sc_preceeding_anchor=6,
    min_anchor_exact=6,
    min_anchor_fuzzy=10,
    min_anchor_match=0.9,
    min_double_aligned_to_estimate_insertion_size=2,
    min_flanking_pairs_resolution=10,
    min_linking_split_reads=2,
    min_mapping_quality=5,
    min_non_target_aligned_split_reads=1,
    min_sample_size_to_apply_percentage=10,
    min_softclipping=6,
    min_spanning_reads_resolution=5,
    min_splits_reads_resolution=3,
    sc_extension_stop=5,
    stdev_count_abnormal=3.0,
    strand_determining_read=2,
    outer_window_min_event_size=125
)
""":class:`MavisNamespace`: holds the settings for computations with the Evidence objects

.. glossary::
    :sorted:

    aligner
        the aligner to use to map the contigs/reads back to the reference e.g blat or bwa

    read_length
        length of reads in the bam file

    stdev_fragment_size
        the standard deviation in insert sizes of paired end reads

    median_fragment_size
        the median insert size of paired end reads

    stdev_count_abnormal
        the number of standard deviations away from the normal considered expected and therefore not qualifying as
        flanking reads

    call_error
        buffer zone for the evidence window

    min_splits_reads_resolution
        minimum number of split reads required to call a breakpoint by split reads

    min_softclipping
        minimum number of soft-clipped bases required for a read to be used as soft-clipped evidence

    min_mapping_quality
        the minimum mapping quality of reads to be used as evidence

    fetch_reads_limit
        maximum number of reads, cap, to loop over for any given evidence window

    fetch_reads_bins
        number of bins to split an evidence window into to ensure more even sampling of high coverage regions

    fetch_min_bin_size
        the minimum size of any bin for reading from a bam file. Increasing this number will result in smaller bins
        being merged or less bins being created (depending on the fetch method)

    filter_secondary_alignments
        filter secondary alignments when gathering read evidence

    min_non_target_aligned_split_reads
        The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local alignment
        to the target region to call a breakpoint by split read evidence

    min_double_aligned_to_estimate_insertion_size
        The minimum number of reads which map soft-clipped to both breakpoints to assume the size of the untemplated
        sequence between the breakpoints is at most the read length - 2 * min_softclipping

    min_linking_split_reads
        The minimum number of split reads which aligned to both breakpoints

    min_flanking_pairs_resolution
        the minimum number of flanking reads required to call a breakpoint by flanking evidence

    assembly_max_paths
        the maximum number of paths to resolve. This is used to limit when there is a messy assembly graph to resolve.
        The assembly will pre-calculate the number of paths (or putative assemblies) and stop if it is greater than
        the given setting.

    assembly_max_kmer_size
        the minimum between this and the smallest length input sequence is used as the kmer size for assembling
        the DeBruijn Graph. If this is not set the default is the 75% of the minimum length input sequence

    assembly_max_kmer_strict
        if set to True then any sequences input to the assembly algorithm that cannot create a kmer of this size will
        be discard. However, if this is set to False, then the kmer size will be reduced accordingly and all input
        sequences will be used in the assembly algorithm

    assembly_min_nc_edge_weight
        Discards all non-cutting edges with a weight/frequency less than this from the DeBruijn graph

    assembly_min_edge_weight
        Discards all edges with a weight/frequency less than this from the DeBruijn graph

    assembly_min_remapped_seq
        The minimum input sequences that must remap for an assembled contig to be used

    assembly_min_remap_coverage
        minimum fraction of the contig sequence which the remapped sequences must align over

    assembly_strand_concordance
        when the number of remapped reads from each strand are compared, the ratio must be above this number to decide
        on the strand

    assembly_include_flanking_pairs
        if true then when the split reads are assembled, any flanking read pairs will also be added

    strand_determining_read
        1 or 2. The read in the pair which determines if (assuming a stranded protocol) the first or second read in the
        pair matches the strand sequenced

    max_sc_preceeding_anchor
        when remapping a softclipped read this determines the amount of softclipping allowed on the side opposite
        of where we expect it. For example for a softclipped read on a breakpoint with a left orientation this limits
        the amount of softclipping that is allowed on the right. If this is set to None then there is no limit on
        softclipping

    blat_min_identity
        the minimum percent identity match required for blat results when aligning contigs

    blat_min_percent_of_max_score
        filter based on blat score values. The best match, highest score, is used as the top value and other results
        must be at least this fraction of the maximum score or they are filtered out

    contig_aln_min_query_consumption
        minimum fraction of the original query sequence that must be used by the read(s) of the alignment

    contig_aln_max_event_size
        relates to determining breakpoints when pairing contig alignments. For any given read in a putative pair the
        soft clipping is extended to include any events of greater than this size. The softclipping is added to the
        side of the alignment as indicated by the breakpoint we are assigning pairs to

    contig_aln_min_anchor_size
        the minimum number of aligned bases for a contig (M or =) in order to simplify. Do not have to be consecutive.

    contig_aln_min_exact_block_event_merge
        the minimum number of exact matches in a read cigar string to stop merging on "events". Events meaning
        insertions, deletions, and mismatches

    outer_window_min_event_size
        the minimum size of an event in order for flanking read evidence to be collected

    contig_aln_merge_inner_anchor
        the minimum number of consecutive exact match base pairs to not merge events within a contig alignment
"""
