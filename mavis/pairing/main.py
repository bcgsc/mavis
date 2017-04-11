import os
import itertools
from Bio import SeqIO
from ..constants import PROTOCOL, COLUMNS, CALL_METHOD, SVTYPE, SPLICE_TYPE
from .constants import DEFAULTS
from ..util import read_inputs, output_tabbed_file, log
from . import equivalent_events


def main(
    inputs, output, annotations, product_fasta_sequence_files=None,
    flanking_call_distance=DEFAULTS.flanking_call_distance,
    split_call_distance=DEFAULTS.split_call_distance,
    contig_call_distance=DEFAULTS.contig_call_distance,
    **kwargs
):
    """
    Args:
        inputs (:class:`List` of :class:`str`): list of input files to read
        output (str): path to the output directory
        flanking_call_distance (int): pairing distance for pairing with an event called by :term:`flanking read pair`
        split_call_distance (int): pairing distance for pairing with an event called by :term:`split read`
        contig_call_distance (int): pairing distance for pairing with an event called by contig or :term:`spanning read`
    """
    if product_fasta_sequence_files is None:
        product_fasta_sequence_files = set()
    else:
        product_fasta_sequence_files = set(product_fasta_sequence_files)
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

    product_sequences = dict()
    # get all sequence ids and all files to read from
    for bpp in bpps:
        if bpp.data[COLUMNS.fusion_sequence_fasta_file]:
            product_fasta_sequence_files.add(bpp.data[COLUMNS.fusion_sequence_fasta_file])
        if bpp.data[COLUMNS.fusion_sequence_fasta_id]:
            product_sequences[bpp.data[COLUMNS.fusion_sequence_fasta_id]] = None
        libraries.add(bpp.data[COLUMNS.library])
    
    # load sequences from all files detected
    log('detected', len(libraries), 'libraries')
    for fname in sorted(list(product_fasta_sequence_files)):
        log('loading:', fname)
        try:
            with open(fname, 'rU') as fh:
                temp = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
                for fid, fseq in temp.items():
                    if fid in product_sequences and product_sequences[fid] is not None and \
                            product_sequences[fid] != fseq:
                        raise AssertionError('sequence identifiers are not unique', fid, fseq, product_sequences[fid])
                    product_sequences[fid] = fseq
        except IOError as err:
            log('failed for open input file', err)
    
    # ensure that all sequences have been found
    for seqid, seq in product_sequences.items():
        if seq is None:
            raise KeyError('failed to find sequence for the product', seqid, seq)
        product_sequences[seqid] = str(seq.seq)
    
    # load all transcripts
    reference_transcripts = dict()
    for chr, genes in annotations.items():
        for gene in genes:
            for t in gene.transcripts:
                if t.name in reference_transcripts:
                    raise KeyError('transcript name is not unique', gene, t)
                reference_transcripts[t.name] = t

    # map the calls by library and ensure there are no name/key conflicts
    calls_by_lib = dict()
    bpp_by_product_key = dict()
    pairings = dict()
    
    for bpp in bpps:
        lib = bpp.data[COLUMNS.library]
        product_key = '_'.join([str(v) for v in [
            bpp.data[COLUMNS.library],
            bpp.data[COLUMNS.protocol],
            bpp.data[COLUMNS.annotation_id],
            bpp.data[COLUMNS.fusion_splicing_pattern],
            bpp.data[COLUMNS.fusion_cdna_coding_start],
            bpp.data[COLUMNS.fusion_cdna_coding_end]
        ]])
        bpp.data[COLUMNS.product_id] = product_key
        calls_by_lib.setdefault(lib, set())
        
        pairings[product_key] = set()
        
        if product_key in calls_by_lib[lib]:
            raise KeyError('duplicate bpp is not unique within lib', lib, product_key, bpp, bpp.data)
        else:
            calls_by_lib[lib].add(product_key)
        bpp_by_product_key[product_key] = bpp
    
    # pairwise comparison of breakpoints between all libraries
    for lib1, lib2 in itertools.combinations(calls_by_lib.keys(), 2):
        # for each two libraries pair all calls
        calls1 = calls_by_lib[lib1]
        calls2 = calls_by_lib[lib2]

        log(len(calls1) * len(calls2), 'comparison(s) between', lib1, 'and', lib2)
        
        for product_key1, product_key2 in itertools.product(calls1, calls2):
            print('testing pairing', product_key1, product_key2)
            if equivalent_events(
                bpp_by_product_key[product_key1],
                bpp_by_product_key[product_key2],
                DISTANCES=DISTANCES,
                reference_transcripts=reference_transcripts,
                product_sequences=product_sequences
            ):
                print('successful pairing')
                pairings[product_key1].add(product_key2)
                pairings[product_key2].add(product_key1)
            else:
                print('did not pair')
    
    for product_key, paired_product_keys in pairings.items():
        bpp = bpp_by_product_key[product_key]
        
        # filter any matches where genes match but transcripts do not
        filtered = []
        for paired_product_key in paired_product_keys:
            paired_bpp = bpp_by_product_key[paired_product_key]
            
            if bpp.data[COLUMNS.gene1] and bpp.data[COLUMNS.gene1] == paired_bpp.data[COLUMNS.gene1]:
                if bpp.data[COLUMNS.transcript1] != paired_bpp.data[COLUMNS.transcript1]:
                    continue
            if bpp.data[COLUMNS.gene2] and bpp.data[COLUMNS.gene2] == paired_bpp.data[COLUMNS.gene2]:
                if bpp.data[COLUMNS.transcript2] != paired_bpp.data[COLUMNS.transcript2]:
                    continue
            filtered.append(paired_product_key)
        bpp.data[COLUMNS.pairing] = ';'.join(sorted(filtered))

    fname = os.path.join(
        output,
        'mavis_paired_{}.tab'.format('_'.join(sorted(list(libraries))))
    )
    output_tabbed_file(bpps, fname)
