import os
import sys
import itertools
from Bio import SeqIO


# local modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from ..constants import PROTOCOL, COLUMNS, CALL_METHOD, SVTYPE, SPLICE_TYPE
from ..pipeline.util import read_inputs, output_tabbed_file, log
from . import equivalent_events


def main(inputs, output, flanking_call_distance, split_call_distance, contig_call_distance, **kwargs):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        flanking_call_distance (int): pairing distance for pairing with an event called by :term:`flanking read pair`
        split_call_distance (int): pairing distance for pairing with an event called by :term:`split read`
        contig_call_distance (int): pairing distance for pairing with an event called by contig or :term:`spanning read`
    """
    # load the file
    DISTANCES = {
        CALL_METHOD.FLANK: flanking_call_distance,
        CALL_METHOD.SPLIT: split_call_distance,
        CALL_METHOD.CONTIG: contig_call_distance
    }

    bpps = []
    bpps.extend(read_inputs(
        inputs,
        require=[
            COLUMNS.cluster_id,
            COLUMNS.validation_id,
            COLUMNS.annotation_id,
            COLUMNS.library,
            COLUMNS.fusion_cdna_coding_start,
            COLUMNS.fusion_cdna_coding_end,
            COLUMNS.fusion_sequence_fasta_id,
            COLUMNS.fusion_sequence_fasta_file
        ],
        in_={
            COLUMNS.protocol: PROTOCOL,
            COLUMNS.event_type: SVTYPE,
            COLUMNS.break1_call_method: CALL_METHOD,
            COLUMNS.break2_call_method: CALL_METHOD,
            COLUMNS.fusion_splicing_pattern: SPLICE_TYPE.values() + [None, 'None']
        },
        add={
            COLUMNS.fusion_cdna_coding_start: None,
            COLUMNS.fusion_cdna_coding_end: None,
            COLUMNS.fusion_sequence_fasta_id: None,
            COLUMNS.fusion_sequence_fasta_file: None,
            COLUMNS.fusion_splicing_pattern: None
        }
    ))
    log('read {} breakpoint pairs'.format(len(bpps)))
    libraries = set()

    SEQUENCES = dict()
    sequence_files = set()
    for bpp in bpps:
        if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
            sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
        libraries.add(bpp.data[COLUMNS.library])
    log('pairing between', len(libraries), 'libraries')
    for f in sorted(list(sequence_files)):
        log('loading:', f)
        with open(f, 'rU') as fh:
            temp = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
            for k in temp:
                if k in SEQUENCES:
                    raise AssertionError('sequence identifiers are not unique', k)
            SEQUENCES.update(temp)

    TRANSCRIPTS = dict()

    for chr, genes in args.annotations[1].items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in TRANSCRIPTS:
                    raise AssertionError('transcript is not unique', gene, t)
                TRANSCRIPTS[t.name] = t

    calls_by_lib = dict()
    pairkey_bpp_mapping = dict()
    pairing = dict()

    for bpp in bpps:
        lib_key = (bpp.data[COLUMNS.library], bpp.data[COLUMNS.protocol])
        pair_key = [
            bpp.data[COLUMNS.library],
            bpp.data[COLUMNS.protocol],
            bpp.data[COLUMNS.annotation_id],
            bpp.data[COLUMNS.fusion_splicing_pattern],
            bpp.data[COLUMNS.fusion_cdna_coding_start],
            bpp.data[COLUMNS.fusion_cdna_coding_end]
        ]
        pair_key = '_'.join([str(k) for k in pair_key if k is not None])
        bpp.data[COLUMNS.product_id] = pair_key
        calls_by_lib.setdefault(lib_key, set())
        calls_by_lib[lib_key].add(pair_key)

        if pair_key in calls_by_lib:
            raise KeyError('duplicate bpp is not unique within lib', pair_key, bpp, bpp.data)
        pairkey_bpp_mapping[pair_key] = bpp
        pairing[pair_key] = set()

    # pairwise comparison of breakpoints between all libraries
    for lib1, lib2 in itertools.combinations(calls_by_lib.keys(), 2):
        # for each two libraries pair all calls
        log(len(calls_by_lib[lib1]) * len(calls_by_lib[lib2]), 'comparison(s) between', lib1, 'and', lib2)
        for pkey1, pkey2 in itertools.product(calls_by_lib[lib1], calls_by_lib[lib2]):
            if equivalent_events(
                pairkey_bpp_mapping[pkey1],
                pairkey_bpp_mapping[pkey2],
                DISTANCES=DISTANCES,
                TRANSCRIPTS=TRANSCRIPTS,
                SEQUENCES=SEQUENCES
            ):
                pairing[pkey1].add(pkey2)
                pairing[pkey2].add(pkey1)
    
    for pkey, pairs in pairing.items():
        bpp = pairkey_bpp_mapping[pkey]
        # filter any matches where genes match but transcripts do not
        filtered = []
        for pkey2 in pairs:
            pair_bpp = pairkey_bpp_mapping[pkey2]
            if bpp.data[COLUMNS.gene1] and bpp.data[COLUMNS.gene1] == pair_bpp.data[COLUMNS.gene1]:
                if bpp.data[COLUMNS.transcript1] != pair_bpp.data[COLUMNS.transcript1]:
                    continue
            if bpp.data[COLUMNS.gene2] and bpp.data[COLUMNS.gene2] == pair_bpp.data[COLUMNS.gene2]:
                if bpp.data[COLUMNS.transcript2] != pair_bpp.data[COLUMNS.transcript2]:
                    continue
            filtered.append(pkey2)
        bpp.data[COLUMNS.pairing] = ';'.join(sorted(filtered))

    fname = os.path.join(
        output,
        'mavis_paired_{}.tab'.format('_'.join(sorted([l for l, p in calls_by_lib])))
    )
    output_tabbed_file(bpps, fname)
