import os
import unittest

from mavis.annotate.variant import annotate_events, Annotation, flatten_fusion_transcript, call_protein_indel, IndelCall
from mavis.annotate.fusion import FusionTranscript
from mavis.annotate.constants import SPLICE_TYPE
from mavis.breakpoint import Breakpoint, BreakpointPair
from mavis.constants import ORIENT, PROTOCOL, STRAND, SVTYPE

from . import get_example_genes, MockObject, MockLongString


def get_best(gene):
    for transcript in gene.unspliced_transcripts:
        if transcript.is_best_transcript:
            return transcript
    raise KeyError('no best transcript for gene', gene)


class TestNDUFA12(unittest.TestCase):

    def setUp(self):
        print(get_example_genes().keys())
        self.gene = get_example_genes()['NDUFA12']
        self.reference_annotations = {self.gene.chr: [self.gene]}
        self.reference_genome = {self.gene.chr: MockObject(
            seq=MockLongString(self.gene.seq, offset=self.gene.start - 1)
        )}
        self.best = get_best(self.gene)

    def test_annotate_events_synonymous(self):
        for gene_list in self.reference_annotations.values():
            for gene in gene_list:
                for t in gene.transcripts:
                    print(t)
        b1 = Breakpoint(self.gene.chr, 95344068, orient=ORIENT.LEFT, strand=STRAND.NS)
        b2 = Breakpoint(self.gene.chr, 95344379, orient=ORIENT.RIGHT, strand=STRAND.NS)
        bpp = BreakpointPair(
            b1, b2, stranded=False, opposing_strands=False, event_type=SVTYPE.DEL, protocol=PROTOCOL.GENOME,
            untemplated_seq='')
        annotations = annotate_events([bpp], reference_genome=self.reference_genome, annotations=self.reference_annotations)
        ann = annotations[0]
        for a in annotations:
            print(a, a.fusion, a.fusion.transcripts)
            print(a.transcript1, a.transcript1.transcripts)
        fseq = ann.fusion.transcripts[0].get_seq()
        refseq = ann.transcript1.transcripts[0].get_seq(self.reference_genome)
        self.assertEqual(refseq, fseq)
        self.assertEqual(1, len(annotations))


class TestARID1B(unittest.TestCase):
    def setUp(self):
        self.gene = get_example_genes()['ARID1B']
        self.reference_annotations = {self.gene.chr: [self.gene]}
        self.reference_genome = {self.gene.chr: MockObject(
            seq=MockLongString(self.gene.seq, offset=self.gene.start - 1)
        )}
        self.best = get_best(self.gene)

    def test_small_duplication(self):
        bpp = BreakpointPair(
            Breakpoint('6', 157100005, strand='+', orient='R'),
            Breakpoint('6', 157100007, strand='+', orient='L'),
            event_type=SVTYPE.DUP,
            untemplated_seq='',
            protocol=PROTOCOL.GENOME
        )
        # annotate the breakpoint with the gene
        annotations = annotate_events([bpp], reference_genome=self.reference_genome, annotations=self.reference_annotations)
        self.assertEqual(1, len(annotations))

        ann = Annotation(bpp, transcript1=self.best, transcript2=self.best)
        ft = FusionTranscript.build(
            ann, self.reference_genome,
            min_orf_size=300, max_orf_cap=10, min_domain_mapping_match=0.9
        )
        ref_tx = self.best.translations[0]
        fusion_tx = ft.translations[0]

        # compare the fusion translation to the refernece translation to create the protein notation
        ref_aa_seq = ref_tx.get_aa_seq(self.reference_genome)
        call = IndelCall(ref_aa_seq, fusion_tx.get_aa_seq())
        self.assertTrue(call.is_dup)

        notation = call_protein_indel(ref_tx, fusion_tx, self.reference_genome)
        print(notation)
        self.assertEqual('ENST00000346085:p.G319dupG', notation)


class TestSVEP1(unittest.TestCase):

    def setUp(self):
        self.gene = get_example_genes()['SVEP1']
        self.reference_annotations = {self.gene.chr: [self.gene]}
        self.reference_genome = {self.gene.chr: MockObject(
            seq=MockLongString(self.gene.seq, offset=self.gene.start - 1)
        )}
        self.best = get_best(self.gene)

    def test_annotate_small_intronic_inversion(self):
        bpp = BreakpointPair(
            Breakpoint('9', 113152627, 113152627, orient='L'),
            Breakpoint('9', 113152635, 113152635, orient='L'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq=''
        )
        annotations = annotate_events([bpp], reference_genome=self.reference_genome, annotations=self.reference_annotations)
        for a in annotations:
            print(a, a.transcript1, a.transcript2)
        self.assertEqual(1, len(annotations))
        ann = annotations[0]
        self.assertEqual(self.best, ann.transcript1)
        self.assertEqual(self.best, ann.transcript2)
        refseq = self.best.transcripts[0].get_seq(self.reference_genome)
        self.assertEqual(1, len(ann.fusion.transcripts))
        self.assertEqual(refseq, ann.fusion.transcripts[0].get_seq())

    def test_build_single_transcript_inversion(self):
        bpp = BreakpointPair(
            Breakpoint('9', 113152627, 113152627, orient='L'),
            Breakpoint('9', 113152635, 113152635, orient='L'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq=''
        )
        ann = Annotation(bpp, transcript1=self.best, transcript2=self.best)
        ft = FusionTranscript.build(
            ann, self.reference_genome,
            min_orf_size=300, max_orf_cap=10, min_domain_mapping_match=0.9
        )
        refseq = self.best.transcripts[0].get_seq(self.reference_genome)
        self.assertEqual(1, len(ft.transcripts))
        self.assertEqual(refseq, ft.transcripts[0].get_seq())


class TestPRKCB(unittest.TestCase):
    def setUp(self):
        self.gene = get_example_genes()['PRKCB']
        self.reference_annotations = {self.gene.chr: [self.gene]}
        self.reference_genome = {self.gene.chr: MockObject(
            seq=MockLongString(self.gene.seq, offset=self.gene.start - 1)
        )}
        self.best = get_best(self.gene)

    def test_retained_intron(self):
        bpp = BreakpointPair(
            Breakpoint('16', 23957049, orient='L'),
            Breakpoint('16', 23957050, orient='R'),
            opposing_strands=False,
            stranded=False,
            event_type=SVTYPE.INS,
            protocol=PROTOCOL.TRANS,
            untemplated_seq='A'
        )
        ann = Annotation(bpp, transcript1=self.best, transcript2=self.best)
        ft = FusionTranscript.build(
            ann, self.reference_genome,
            min_orf_size=300, max_orf_cap=10, min_domain_mapping_match=0.9
        )
        self.assertEqual(1, len(ft.transcripts))
        print(ft.transcripts[0].splicing_pattern)
        print(self.best.transcripts[0].splicing_pattern)
        self.assertEqual(SPLICE_TYPE.RETAIN, ft.transcripts[0].splicing_pattern.splice_type)


class TestDSTYK(unittest.TestCase):
    def setUp(self):
        print(get_example_genes().keys())
        self.gene = get_example_genes()['DSTYK']
        self.reference_annotations = {self.gene.chr: [self.gene]}
        self.reference_genome = {self.gene.chr: MockObject(
            seq=MockLongString(self.gene.seq, offset=self.gene.start - 1)
        )}
        self.best = get_best(self.gene)

    def test_build_single_transcript_inversion_reverse_strand(self):
        # 1:205178631R 1:205178835R inversion
        bpp = BreakpointPair(
            Breakpoint('1', 205178631, orient='R'),
            Breakpoint('1', 205178835, orient='R'),
            opposing_strands=True,
            stranded=False,
            event_type=SVTYPE.INV,
            protocol=PROTOCOL.GENOME,
            untemplated_seq=''
        )
        ann = Annotation(bpp, transcript1=self.best, transcript2=self.best)
        ft = FusionTranscript.build(
            ann, self.reference_genome,
            min_orf_size=300, max_orf_cap=10, min_domain_mapping_match=0.9
        )
        print(ft.exons)
        print(ft.break1, ft.break2)
        for ex in ft.exons:
            print(ex, len(ex), '==>', ft.exon_mapping.get(ex.position, None), len(ft.exon_mapping.get(ex.position, None)), ft.exon_number(ex))
        # refseq = self.best.transcripts[0].get_seq(self.reference_genome)
        self.assertEqual(1, len(ft.transcripts))
        self.assertEqual(1860, ft.break1)
        self.assertEqual(2065, ft.break2)
        flatten_fusion_transcript(ft.transcripts[0])  # test no error
