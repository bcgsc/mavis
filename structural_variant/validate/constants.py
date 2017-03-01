from argparse import Namespace

VALIDATION_DEFAULTS = Namespace(
    assembly_include_flanking_pairs=True,
    assembly_include_half_mapped_reads=True,
    assembly_max_paths=20,
    assembly_min_edge_weight=3,
    assembly_min_remap=3,
    assembly_min_tgt_to_exclude_half_map=7,
    assembly_strand_concordance=0.51,
    assembly_max_kmer_size=None,
    assembly_max_kmer_strict=True,
    call_error=10,
    consensus_req=3,
    fetch_reads_bins=3,
    fetch_reads_limit=3000,
    filter_secondary_alignments=True,
    fuzzy_mismatch_number=1,
    max_sc_preceeding_anchor=6,
    min_anchor_exact=6,
    min_anchor_fuzzy=10,
    min_anchor_match=0.9,
    min_flanking_pairs_resolution=3,
    min_linking_split_reads=2,
    min_mapping_quality=20,
    min_non_target_aligned_split_reads=1,
    min_sample_size_to_apply_percentage=10,
    min_softclipping=6,
    min_splits_reads_resolution=3,
    min_double_aligned_to_estimate_insertion_size=2,
    sc_extension_stop=5,
    stdev_count_abnormal=3,
    strand_determining_read=2
)
""":class:`~argparse.Namespace`: holds the settings for computations with the Evidence objects

.. glossary::
    :sorted:

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

    filter_secondary_alignments
        filter secondary alignments when gathering read evidence

    assembly_min_edge_weight
        when building the initial deBruijn graph edge weights are determined by the frequency of the kmer they represent
        in all the input sequences. The parameter here discards edges to the kmer they represent in all the input
        sequences. The parameter here discards edges to simply the graph if they have a weight less than specified

    assembly_min_remap
        The minimum input sequences that must remap for an assembly to be used

    min_non_target_aligned_split_reads
        The minimum number of split reads aligned to a breakpoint by the input bam and no forced by local alignment
        to the target region to call a breakpoint by split read evidence

    min_double_aligned_to_estimate_insertion_size
        The minimum number of reads which map soft-clipped to both breakpoints to assume the size of the untemplated 
        sequence between the breakpoints is at most the read length - 2 * min_softclipping

    min_linking_split_reads
        The minimum number of split reads which aligned to both breakpoints

    min_flanking_pairs_resolutionraw_break1_half_mapped_reads
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
"""

