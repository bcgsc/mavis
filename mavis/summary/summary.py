from .constants import PAIRING_STATE
from ..breakpoint import Breakpoint, BreakpointPair
from ..constants import CALL_METHOD, COLUMNS, DISEASE_STATUS, PROTOCOL, SVTYPE
from ..interval import Interval
from ..pairing.pairing import pair_by_distance, product_key
from ..util import get_connected_components


def filter_by_annotations(bpp_list, best_transcripts):
    """
    Args:
        bpp_list (list of BreakpointPair): list of pairs to filter
        best_transcripts (:class `dict` of :any:`Transcript` by :class:`str`): the best transcripts of the annotations
          based on their names

    """
    strings = []
    for bpp in bpp_list:
        for attr in ['gene1', 'gene2', 'transcript1', 'transcript2']:
            if bpp.data[attr] is not None:
                strings.append(bpp.data[attr])
    string_ranks = {s: i for i, s in enumerate(sorted(strings))}
    string_ranks[None] = len(strings)

    def sort_key(bpp):
        if bpp.fusion_cdna_coding_start is None:
            result = [1, 0]
        else:
            result = [0, -1 * (int(bpp.fusion_cdna_coding_end) - int(bpp.fusion_cdna_coding_start))]

        result.extend([
            0 if bpp.transcript1 in best_transcripts else 1,
            0 if bpp.transcript2 in best_transcripts else 1,
            sum([bpp.transcript1 is None, bpp.transcript2 is None]),
            string_ranks[bpp.gene1], string_ranks[bpp.gene2],
            string_ranks[bpp.transcript1], string_ranks[bpp.transcript2]
        ])
        return tuple(result)
    bpp_list = sorted(bpp_list, key=sort_key)
    result = []
    removed = []
    for bpp in bpp_list:
        if sort_key(bpp) == sort_key(bpp_list[0]):
            result.append(bpp)
        else:
            removed.append(bpp)
    return result, removed


def filter_by_call_method(bpp_list):
    """
    Filters a set of breakpoint pairs to returns the call with the most evidence.
    Prefers contig evidence over spanning over split over flanking, etc.
    """
    # ranking scores of the methods (more is better)
    def sort_key(bpp):
        key = [bpp.data.get(col, 0) if bpp.data.get(col, 0) is not None else 0 for col in [
            'contig_remapped_reads',
            'contig_alignment_score',
            'spanning_reads',
            'break1_split_reads',
            'break2_split_reads',
            'linking_split_reads',
            'flanking_pairs'
        ]]
        return tuple(key)
    if not bpp_list:
        return bpp_list
    bpp_list = sorted(bpp_list, key=sort_key, reverse=True)

    # filter to the top ranked method
    result = []
    removed = []
    for bpp in bpp_list:
        if sort_key(bpp) == sort_key(bpp_list[0]):
            result.append(bpp)
        else:
            removed.append(bpp)
    return result, removed


def group_events(events):
    """
    group events together and join data attributes
    """
    # take the outer regions of the breakpoints
    first = events[0]
    new_bpp = BreakpointPair(
        Breakpoint(
            first.break1.chr,
            min([b.break1.start for b in events]),
            max([b.break1.end for b in events]),
            orient=first.break1.orient,
            strand=first.break1.strand),
        Breakpoint(
            first.break2.chr,
            min([b.break2.start for b in events]),
            max([b.break2.end for b in events]),
            orient=first.break2.orient,
            strand=first.break2.strand),
        opposing_strands=first.opposing_strands,
        stranded=first.stranded
    )
    data_columns = set()
    for bpp in events:
        data_columns.update(bpp.data.keys())
        if any([
            bpp.break1.chr != new_bpp.break1.chr,
            bpp.break2.chr != new_bpp.break2.chr,
            bpp.break1.orient != new_bpp.break1.orient,
            bpp.break2.orient != new_bpp.break2.orient,
            bpp.opposing_strands != new_bpp.opposing_strands,
            bpp.break1.strand != new_bpp.break1.strand,
            bpp.break2.strand != new_bpp.break2.strand
        ]):
            raise AssertionError('cannot group events differing on key elements', bpp, new_bpp)

    # Note: There are some attributes that shouldn't be lost if different, currently appending the information
    # The evidence could be better off as a max instead of a join
    list_columns = {
        COLUMNS.contig_seq, COLUMNS.call_method,
        COLUMNS.break1_split_reads, COLUMNS.break2_split_reads, COLUMNS.contig_alignment_score,
        COLUMNS.spanning_reads, COLUMNS.flanking_pairs, COLUMNS.tools,
        COLUMNS.product_id, COLUMNS.event_type, COLUMNS.annotation_id,
        COLUMNS.pairing, COLUMNS.annotation_figure,
        COLUMNS.contig_remapped_reads, COLUMNS.tools,
        COLUMNS.tracking_id
    }
    for col in data_columns:
        new_data = sorted(list({bpp.data[col] for bpp in events}), key=lambda x: str(x))
        if len(new_data) == 1:
            new_bpp.data[col] = new_data[0]
        elif col in list_columns:
            new_bpp.data[col] = ';'.join([str(v) for v in new_data])

    untemplated_seq = {bpp.untemplated_seq for bpp in events}
    if len(untemplated_seq) == 1:
        new_bpp.untemplated_seq = list(untemplated_seq)[0]

    return new_bpp


def group_by_distance(calls, distances):
    """
    groups a set of calls based on their proximity. Returns a new list of calls where close calls have been merged
    """
    mapping = {}
    for call in calls:
        mapping.setdefault(product_key(call), []).append(call)
    pairing = pair_by_distance(calls, distances, against_self=True)
    # merge all the 'close-enough' pairs
    grouped_calls = []
    removed_calls = []
    for component in get_connected_components(pairing):
        if len(component) == 1:
            grouped_calls.extend(mapping[component.pop()])
        else:
            pairs = []
            for key in component:
                pairs.extend(mapping[key])
            grouped_calls.append(group_events(pairs))
            removed_calls.extend(pairs)
    return grouped_calls, removed_calls


def annotate_dgv(bpps, dgv_regions_by_reference_name, distance=0):
    """
    given a list of bpps and a dgv reference, annotate the events that are within the set distance of both breakpoints

    Args:
        bpps (list) : the list of BreakpointPair objects
        dgv_regions_by_reference_name (dict) : the dgv reference regions file loaded by load_masking_regions
        distance (int) : the minimum distance required to match a dgv event with a breakpoint
    """
    for chrom in dgv_regions_by_reference_name:
        dgv_regions_by_reference_name[chrom] = sorted(dgv_regions_by_reference_name[chrom], key=lambda x: x.start)

    lowest_resolution = max([len(b.break1) for b in bpps])  # only need start res

    # only look at the bpps that dgv events could pair to, Intrachromosomal
    for bpp in [b for b in bpps if not b.interchromosomal and b.break1.chr in dgv_regions_by_reference_name]:
        for dgv_region in dgv_regions_by_reference_name[bpp.break1.chr]:
            dist = abs(Interval.dist(Interval(dgv_region.start), bpp.break1))
            if dist > lowest_resolution + distance:
                break
            elif dist > distance or abs(Interval.dist(Interval(dgv_region.end), bpp.break2)) > distance:
                continue
            refname = dgv_region.reference_object
            try:
                refname = dgv_region.reference_object.name
            except AttributeError:
                pass
            bpp.data['dgv'] = '{}({}:{}-{})'.format(dgv_region.name, refname, dgv_region.start, dgv_region.end)


def get_pairing_state(current_protocol, current_disease_state, other_protocol, other_disease_state, is_matched=False, inferred_is_matched=False):
    """
    given two libraries, returns the appropriate descriptor for their matched state

    Args:
        current_protocol (PROTOCOL): the protocol of the current library
        current_disease_state (DISEASE_STATUS): the disease status of the current library
        other_protocol (PROTOCOL): protocol of the library being comparing to
        other_disease_state (DISEASE_STATUS): disease status of the library being compared to
        is_matched (bool): True if the libraries are paired

    Returns:
        (PAIRING_STATE): descriptor of the pairing of the two libraries
    """
    PROTOCOL.enforce(current_protocol)
    PROTOCOL.enforce(other_protocol)
    DISEASE_STATUS.enforce(current_disease_state)
    DISEASE_STATUS.enforce(other_disease_state)

    curr = (current_protocol, current_disease_state)
    other = (other_protocol, other_disease_state)

    dg = (PROTOCOL.GENOME, DISEASE_STATUS.DISEASED)
    dt = (PROTOCOL.TRANS, DISEASE_STATUS.DISEASED)
    ng = (PROTOCOL.GENOME, DISEASE_STATUS.NORMAL)

    if current_protocol != other_protocol:
        is_matched = is_matched or inferred_is_matched

    if curr in {dg, dt} and other == ng:
        return PAIRING_STATE.GERMLINE if is_matched else PAIRING_STATE.SOMATIC
    elif curr in {dg, ng} and other == dt:
        return PAIRING_STATE.EXP if is_matched else PAIRING_STATE.NO_EXP
    elif curr == dt and other == dg:
        return PAIRING_STATE.GENOMIC if is_matched else PAIRING_STATE.NO_GENOMIC
    else:
        return PAIRING_STATE.MATCH if is_matched else PAIRING_STATE.NO_MATCH


def filter_by_evidence(
    bpps,
    filter_min_remapped_reads=5,
    filter_min_spanning_reads=5,
    filter_min_flanking_reads=10,
    filter_min_split_reads=5,
    filter_min_linking_split_reads=1
):
    filtered = []
    removed = []
    for bpp in bpps:
        if bpp.call_method == CALL_METHOD.CONTIG:
            # inherently the breakpoints have been linked
            if int(bpp.contig_remapped_reads) < filter_min_remapped_reads:
                removed.append(bpp)
                continue
        elif bpp.call_method == CALL_METHOD.SPAN:
            if bpp.spanning_reads < filter_min_spanning_reads:
                removed.append(bpp)
                continue
        elif bpp.call_method == CALL_METHOD.SPLIT:
            linking_split_reads = bpp.linking_split_reads
            if bpp.event_type == SVTYPE.INS:
                linking_split_reads += bpp.flanking_pairs
            if any([
                bpp.break1_split_reads + bpp.break1_split_reads_forced < filter_min_split_reads,
                bpp.break2_split_reads + bpp.break2_split_reads_forced < filter_min_split_reads,
                linking_split_reads < filter_min_linking_split_reads,
                bpp.break1_split_reads < 1,
                bpp.break2_split_reads < 1
            ]):
                removed.append(bpp)
                continue
        elif bpp.call_method == CALL_METHOD.FLANK:
            if bpp.flanking_pairs < filter_min_flanking_reads:
                removed.append(bpp)
                continue
        elif bpp.call_method != CALL_METHOD.INPUT:
            raise AssertionError('unexpected value for call_method: {}'.format(
                bpp.call_method))
        filtered.append(bpp)
    return filtered, removed
