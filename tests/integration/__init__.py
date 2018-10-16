import os

from mavis.align import query_coverage_interval
from mavis.annotate.genomic import Transcript, PreTranscript
from mavis.annotate.file_io import load_annotations, load_reference_genome
from mavis.annotate.protein import Translation
from mavis.constants import CIGAR, NA_MAPPING_QUALITY

from ..util import DATA_DIR


ARGUMENT_ERROR = 2

RUN_FULL = int(os.environ.get('RUN_FULL', 1))
OUTPUT_SVG = int(os.environ.get('OUTPUT_SVG', 0))
_EXAMPLE_GENES = None


def get_example_genes():
    global _EXAMPLE_GENES
    if _EXAMPLE_GENES is None:
        _EXAMPLE_GENES = set_example_genes()
    return _EXAMPLE_GENES


def set_example_genes():
    result = {}
    genes = load_annotations(os.path.join(DATA_DIR, 'example_genes.json'))
    seqs = load_reference_genome(os.path.join(DATA_DIR, 'example_genes.fa'))
    for chr_genes in genes.values():
        for gene in chr_genes:
            if gene.name in seqs:
                gene.seq = str(seqs[gene.name].seq)
            result[gene.name] = gene
            if gene.aliases:
                for alias in gene.aliases:
                    result[alias] = gene
    print(result.keys())
    return result


class MockObject:
    def __init__(self, **kwargs):
        for arg, val in kwargs.items():
            setattr(self, arg, val)


class MockRead:
    def __init__(
        self,
        query_name=None,
        reference_id=0,
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
        query_alignment_length=None,
        flag=None,
        tags=[],
        is_read1=True,
        is_paired=True,
        is_unmapped=False,
        mate_is_unmapped=False,
        mapping_quality=NA_MAPPING_QUALITY,
        **kwargs
    ):
        for attr, val in kwargs.items():
            setattr(self, attr, val)
        self.mapping_quality = mapping_quality
        self.query_name = query_name
        self.reference_id = reference_id
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.cigar = cigar
        self.query_alignment_length = query_alignment_length
        self.query_sequence = query_sequence
        if self.reference_end is None and cigar and reference_start is not None:
            self.reference_end = reference_start + sum([f for v, f in cigar if v not in [CIGAR.S, CIGAR.I]])
        if not self.query_alignment_length:
            if cigar:
                self.query_alignment_length = sum([v for s, v in cigar if s not in [CIGAR.S, CIGAR.H]])
            elif self.query_sequence:
                self.query_alignment_length = len(self.query_sequence)
            else:
                self.query_alignment_length = 0
        self.is_reverse = is_reverse
        self.mate_is_reverse = mate_is_reverse
        self.next_reference_start = next_reference_start
        self.next_reference_id = next_reference_id
        self.reference_name = reference_name

        self.query_alignment_sequence = query_alignment_sequence
        self.query_alignment_start = query_alignment_start
        self.query_alignment_end = query_alignment_end
        self.flag = flag
        self.tags = tags
        if query_alignment_sequence is None and cigar and query_sequence:
            s = 0 if cigar[0][0] != CIGAR.S else cigar[0][1]
            t = len(query_sequence)
            if cigar[-1][0] == CIGAR.S:
                t -= cigar[-1][1]
            self.query_alignment_sequence = query_sequence[s:t]
        if cigar and query_sequence:
            if len(query_sequence) != sum([f for v, f in cigar if v not in [CIGAR.H, CIGAR.N, CIGAR.D]]):
                raise AssertionError('length of sequence does not match cigar', len(query_sequence), sum([f for v, f in cigar if v not in [CIGAR.H, CIGAR.N, CIGAR.D]]))
        if template_length is None and reference_end and next_reference_start:
            self.template_length = next_reference_start - reference_end
        else:
            self.template_length = template_length
        self.is_read1 = is_read1
        self.is_read2 = (not is_read1)
        self.is_paired = is_paired
        self.is_unmapped = is_unmapped
        self.mate_is_unmapped = mate_is_unmapped
        if self.reference_start and self.reference_end:
            if not cigar:
                self.cigar = [(CIGAR.M, self.reference_end - self.reference_start)]
            if not self.query_sequence:
                self.query_sequence = 'N' * (self.reference_end - self.reference_start)
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
        return query_coverage_interval(self)

    def set_tag(self, tag, value, value_type=None, replace=True):
        new_tag = (tag, value)
        if not replace and new_tag in self.tags:
            self.tags.append(new_tag)
        else:
            self.tags.append(new_tag)

    def has_tag(self, tag):
        return tag in dict(self.tags).keys()

    def get_tag(self, tag):
        return dict(self.tags)[tag] if tag in dict(self.tags).keys() else False

    def __str__(self):
        return '{}(ref_id={}, start={}, end={}, seq={})'.format(
            self.__class__.__name__, self.reference_id, self.reference_start, self.reference_end, self.query_sequence)


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


class MockString:
    def __init__(self, char=' '):
        self.char = char

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.char * (index.stop - index.start)
        else:
            return self.char


class MockLongString(str):
    def __new__(cls, *args, offset=0, **kw):
        s = str.__new__(cls, *args, **kw)
        setattr(s, 'offset', offset)
        return s

    def __getitem__(self, index):
        if isinstance(index, slice):
            index = slice(index.start - self.offset, index.stop - self.offset, index.step)
        else:
            index -= self.offset
        return str.__getitem__(self, index)

    def __len__(self):
        return self.offset + str.__len__(self)


def mock_read_pair(mock1, mock2):
    if mock1.reference_id != mock2.reference_id:
        mock1.template_length = 0
        mock2.template_length = 0
    mock1.next_reference_id = mock2.reference_id
    mock1.next_reference_start = mock2.reference_start
    mock1.mate_is_reverse = mock2.is_reverse
    mock1.is_paired = True
    mock1.is_read2 = not mock1.is_read1
    if mock1.template_length is None:
        mock1.template_length = mock2.reference_end - mock1.reference_start

    mock2.next_reference_id = mock1.reference_id
    mock2.next_reference_start = mock1.reference_start
    mock2.mate_is_reverse = mock1.is_reverse
    mock2.is_paired = True
    mock2.is_read1 = not mock1.is_read1
    mock2.is_read2 = not mock1.is_read2
    if mock2.query_name is None:
        mock2.query_name = mock1.query_name
    mock2.template_length = -1 * mock1.template_length
    return mock1, mock2


def build_transcript(gene, exons, cds_start, cds_end, domains, strand=None, is_best_transcript=False, name=None):
    pre_transcript = PreTranscript(exons, gene=gene, strand=strand if strand is not None else gene.get_strand(), is_best_transcript=is_best_transcript, name=name)
    if gene is not None:
        gene.unspliced_transcripts.append(pre_transcript)

    for spl in pre_transcript.generate_splicing_patterns():
        t = Transcript(pre_transcript, spl)
        pre_transcript.spliced_transcripts.append(t)

        tx = Translation(cds_start, cds_end, t, domains=domains)
        t.translations.append(tx)

    return pre_transcript
