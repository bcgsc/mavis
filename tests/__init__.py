from structural_variant.constants import CIGAR
from structural_variant.blat import BlatAlignedSegment
from structural_variant.annotate.genomic import usTranscript, Transcript
from structural_variant.annotate.protein import Translation
import os


filedir = os.path.join(os.path.dirname(__file__), 'files')
REFERENCE_GENOME_FILE = os.path.join(filedir, 'mock_reference_genome.fa')
REFERENCE_GENOME_FILE_2BIT = os.path.join(filedir, 'mock_reference_genome.2bit')
REFERENCE_ANNOTATIONS_FILE = os.path.join(filedir, 'mock_reference_annotations.tsv')
REFERENCE_ANNOTATIONS_FILE_JSON = os.path.join(filedir, 'mock_reference_annotations.json')
TEMPLATE_METADATA_FILE = os.path.join(filedir, 'cytoBand.txt')
BAM_INPUT = os.path.join(filedir, 'mini_mock_reads_for_events.sorted.bam')
FULL_BAM_INPUT = os.path.join(filedir, 'mock_reads_for_events.sorted.bam')
FULL_BASE_EVENTS = os.path.join(filedir, 'mock_sv_events.svmerge.tsv')
BASE_EVENTS = os.path.join(filedir, 'mini_mock_sv_events.svmerge.tsv')
BLAT_INPUT = os.path.join(filedir, 'blat_input.fa')
BLAT_OUTPUT = os.path.join(filedir, 'blat_output.pslx')


class MockRead:
    def __init__(
        self,
        query_name=None,
        reference_id=None,
        reference_start=None,
        reference_end=None,
        cigar=None,
        is_reverse=False,
        mate_is_reverse=True,
        next_reference_start=None,
        next_reference_id=None,
        reference_name=None,
        query_sequence=None,
        template_length=None,
        query_alignment_sequence=None,
        query_alignment_start=None,
        query_alignment_end=None,
        flag=None,
        tags=[],
        is_read1=True,
        is_paired=True,
        is_unmapped=False,
        mate_is_unmapped=False,
        **kwargs
    ):
        args = {
            is_reverse: False, mate_is_reverse: True,
            is_unmapped: False, mate_is_unmapped: False,
            next_reference_end: None, reference_end: None,
            reference_start: None, next_reference_start: None,
            template_length: None, query_sequence: None
        }
        kwargs.setdefault('is_read2', not kwargs.get('is_read1', True))
        kwargs.setdefault('is_read1', not kwargs['is_read2'])
        args.update(kwargs)

        for attr, val in args.items():
            setattr(self, attr, val)

        if self.reference_end is None and self.cigar and self.reference_start is not None:
            self.reference_end = reference_start + sum([f for v, f in cigar if v not in [CIGAR.S, CIGAR.I]])

        if self.query_alignment_sequence is None and self.cigar and self.query_sequence:
            s = 0 if self.cigar[0][0] != CIGAR.S else self.cigar[0][1]
            t = len(self.query_sequence)
            if self.cigar[-1][0] == CIGAR.S:
                t -= self.cigar[-1][1]
            self.query_alignment_sequence = self.query_sequence[s:t]
        if self.cigar and self.query_sequence:
            assert(len(self.query_sequence) == sum([f for v, f in self.cigar if v not in [CIGAR.H, CIGAR.N, CIGAR.D]]))
        if self.template_length is None and self.reference_end and self.next_reference_start:
            self.template_length = next_reference_start - reference_end
        else:
            self.template_length = template_length

        if flag:
            self.is_unmapped = bool(self.flag & int(0x4))
            self.mate_is_unmapped = bool(self.flag & int(0x8))
            self.is_reverse = bool(self.flag & int(0x10))
            self.mate_is_reverse = bool(self.flag & int(0x20))
            self.is_read1 = bool(self.flag & int(0x40))
            self.is_read2 = bool(self.flag & int(0x80))
            self.is_secondary = bool(self.flag & int(0x100))
            self.is_qcfail = bool(self.flag & int(0x200))
            self.is_supplementary = bool(self.flag & int(0x400))

    def query_coverage_interval(self):
        return BlatAlignedSegment.query_coverage_interval(self)

    def set_tag(self, tag, value, value_type=None, replace=True):
        new_tag=(tag,value)
        if not replace and new_tag in self.tags:
            self.tags.append(new_tag)
        else:
            self.tags.append(new_tag)

    def has_tag(self, tag):
        return tag in dict(self.tags).keys()

    def get_tag(self, tag):
        return dict(self.tags)[tag] if tag in dict(self.tags).keys() else False

class MockBamFileHandle:
    def __init__(self, chrom_to_tid={}):
        self.chrom_to_tid = chrom_to_tid

    def fetch(self, *pos):
        return []

    def get_tid(self, chrom):
        if chrom in self.chrom_to_tid:
            return self.chrom_to_tid[chrom]
        else:
            return -1

    def get_reference_name(self, input_tid):
        for chrom, tid in self.chrom_to_tid.items():
            if input_tid == tid:
                return chrom
        raise KeyError('invalid id')


class MockSeq:
    def __init__(self, seq=None):
        self.seq = seq


class MockString:
    def __init__(self, char=' '):
        self.char = char

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.char * (index.stop - index.start)
        else:
            return self.char


def build_transcript(gene, exons, cds_start, cds_end, domains, strand=None, is_best_transcript=False):
    ust = usTranscript(exons, gene=gene, strand=strand if strand is not None else gene.get_strand(), is_best_transcript=is_best_transcript)
    if gene is not None:
        gene.unspliced_transcripts.append(ust)

    for spl in ust.generate_splicing_patterns():
        t = Transcript(ust, spl)
        ust.spliced_transcripts.append(t)

        tx = Translation(cds_start, cds_end, t, domains=domains)
        t.translations.append(tx)

    return ust
