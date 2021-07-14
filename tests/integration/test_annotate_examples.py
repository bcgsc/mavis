from mavis.annotate.fusion import FusionTranscript
from mavis.annotate.variant import (
    Annotation,
    IndelCall,
    annotate_events,
    call_protein_indel,
    flatten_fusion_transcript,
)
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, PROTOCOL, SPLICE_TYPE, STRAND, SVTYPE

from ..util import long_running_test
from . import MockLongString, MockObject, get_example_genes


def get_best(gene):
    for transcript in gene.unspliced_transcripts:
        if transcript.is_best_transcript:
            return transcript
    raise KeyError('no best transcript for gene', gene)


class TestNDUFA12:
    def test_annotate_events_synonymous(self):
        gene = get_example_genes()['NDUFA12']
        reference_annotations = {gene.chr: [gene]}
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }

        for gene_list in reference_annotations.values():
            for gene in gene_list:
                for t in gene.transcripts:
                    print(t)
        b1 = Breakpoint(gene.chr, 95344068, orient=ORIENT.LEFT, strand=STRAND.NS)
        b2 = Breakpoint(gene.chr, 95344379, orient=ORIENT.RIGHT, strand=STRAND.NS)
        bpp = BreakpointPair(
            b1,
            b2,
            stranded=False,
            opposing_strands=False,
            event_type=SVTYPE.DEL,
            protocol=PROTOCOL.GENOME,
            untemplated_seq='',
        )
        annotations = annotate_events(
            [bpp], reference_genome=reference_genome, annotations=reference_annotations
        )
        ann = annotations[0]
        for a in annotations:
            print(a, a.fusion, a.fusion.transcripts)
            print(a.transcript1, a.transcript1.transcripts)
        fseq = ann.fusion.transcripts[0].get_seq()
        refseq = ann.transcript1.transcripts[0].get_seq(reference_genome)
        assert fseq == refseq
        assert len(annotations) == 1


class TestARID1B:
    def test_small_duplication(self):
        gene = get_example_genes()['ARID1B']
        reference_annotations = {gene.chr: [gene]}
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }
        best = get_best(gene)

        bpp = BreakpointPair(
            Breakpoint('6', 157100005, strand='+', orient='R'),
            Breakpoint('6', 157100007, strand='+', orient='L'),
            event_type=SVTYPE.DUP,
            untemplated_seq='',
            protocol=PROTOCOL.GENOME,
        )
        # annotate the breakpoint with the gene
        annotations = annotate_events(
            [bpp], reference_genome=reference_genome, annotations=reference_annotations
        )
        assert len(annotations) == 1

        ann = Annotation(bpp, transcript1=best, transcript2=best)
        ft = FusionTranscript.build(
            ann,
            reference_genome,
            min_orf_size=300,
            max_orf_cap=10,
            min_domain_mapping_match=0.9,
        )
        ref_tx = best.translations[0]
        fusion_tx = ft.translations[0]

        # compare the fusion translation to the refernece translation to create the protein notation
        ref_aa_seq = ref_tx.get_aa_seq(reference_genome)
        call = IndelCall(ref_aa_seq, fusion_tx.get_aa_seq())
        assert call.is_dup

        notation = call_protein_indel(ref_tx, fusion_tx, reference_genome)
        print(notation)
        assert notation == 'ENST00000346085:p.G319dupG'


class TestSVEP1:
    @long_running_test
    def test_annotate_small_intronic_inversion(self):
        gene = get_example_genes()['SVEP1']
        reference_annotations = {gene.chr: [gene]}
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }
        best = get_best(gene)

        bpp = BreakpointPair(
            Breakpoint('9', 113152627, 113152627, orient='L'),
            Breakpoint('9', 113152635, 113152635, orient='L'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq='',
        )
        annotations = annotate_events(
            [bpp], reference_genome=reference_genome, annotations=reference_annotations
        )
        for a in annotations:
            print(a, a.transcript1, a.transcript2)
        assert len(annotations) == 1
        ann = annotations[0]
        assert ann.transcript1 == best
        assert ann.transcript2 == best
        refseq = best.transcripts[0].get_seq(reference_genome)
        assert len(ann.fusion.transcripts) == 1
        assert ann.fusion.transcripts[0].get_seq() == refseq

    @long_running_test
    def test_build_single_transcript_inversion(self):
        gene = get_example_genes()['SVEP1']
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }
        best = get_best(gene)

        bpp = BreakpointPair(
            Breakpoint('9', 113152627, 113152627, orient='L'),
            Breakpoint('9', 113152635, 113152635, orient='L'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq='',
        )
        ann = Annotation(bpp, transcript1=best, transcript2=best)
        ft = FusionTranscript.build(
            ann,
            reference_genome,
            min_orf_size=300,
            max_orf_cap=10,
            min_domain_mapping_match=0.9,
        )
        refseq = best.transcripts[0].get_seq(reference_genome)
        assert len(ft.transcripts) == 1
        assert ft.transcripts[0].get_seq() == refseq


class TestPRKCB:
    def test_retained_intron(self):
        gene = get_example_genes()['PRKCB']
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }
        best = get_best(gene)

        bpp = BreakpointPair(
            Breakpoint('16', 23957049, orient='L'),
            Breakpoint('16', 23957050, orient='R'),
            opposing_strands=False,
            stranded=False,
            event_type=SVTYPE.INS,
            protocol=PROTOCOL.TRANS,
            untemplated_seq='A',
        )
        ann = Annotation(bpp, transcript1=best, transcript2=best)
        ft = FusionTranscript.build(
            ann,
            reference_genome,
            min_orf_size=300,
            max_orf_cap=10,
            min_domain_mapping_match=0.9,
        )
        assert len(ft.transcripts) == 1
        print(ft.transcripts[0].splicing_pattern)
        print(best.transcripts[0].splicing_pattern)
        assert ft.transcripts[0].splicing_pattern.splice_type == SPLICE_TYPE.RETAIN


class TestDSTYK:
    def test_build_single_transcript_inversion_reverse_strand(self):
        print(get_example_genes().keys())
        gene = get_example_genes()['DSTYK']
        reference_genome = {
            gene.chr: MockObject(seq=MockLongString(gene.seq, offset=gene.start - 1))
        }
        best = get_best(gene)

        # 1:205178631R 1:205178835R inversion
        bpp = BreakpointPair(
            Breakpoint('1', 205178631, orient='R'),
            Breakpoint('1', 205178835, orient='R'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq='',
        )
        ann = Annotation(bpp, transcript1=best, transcript2=best)
        ft = FusionTranscript.build(
            ann,
            reference_genome,
            min_orf_size=300,
            max_orf_cap=10,
            min_domain_mapping_match=0.9,
        )
        print(ft.exons)
        print(ft.break1, ft.break2)
        for ex in ft.exons:
            print(
                ex,
                len(ex),
                '==>',
                ft.exon_mapping.get(ex.position, None),
                len(ft.exon_mapping.get(ex.position, None)),
                ft.exon_number(ex),
            )
        # refseq = best.transcripts[0].get_seq(reference_genome)
        assert len(ft.transcripts) == 1
        assert ft.break1 == 1860
        assert ft.break2 == 2065
        flatten_fusion_transcript(ft.transcripts[0])  # test no error
