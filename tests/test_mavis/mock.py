import os
import types
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from mavis.annotate.file_io import load_annotations, load_reference_genome
from mavis.annotate.genomic import PreTranscript, Transcript
from mavis.annotate.protein import Translation
from mavis.bam.cigar import (
    QUERY_ALIGNED_STATES,
    REFERENCE_ALIGNED_STATES,
    convert_cigar_to_string,
    convert_string_to_cigar,
)
from mavis.bam.read import SamRead
from mavis.constants import CIGAR, NA_MAPPING_QUALITY, PYSAM_READ_FLAGS

from ..util import get_data

ARGUMENT_ERROR = 2

RUN_FULL = int(os.environ.get('RUN_FULL', 1))
OUTPUT_SVG = int(os.environ.get('OUTPUT_SVG', 0))
_EXAMPLE_GENES = None


class Mock:
    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def bind_method(self, **kwargs):
        for attr, val in kwargs.items():
            val = types.MethodType(val, self)  # bind the method to self
            setattr(self, attr, val)

    def add_attr(self, attr, val):
        setattr(self, attr, val)

    def __contains__(self, item):
        if hasattr(self, item):
            return True
        return False


class MockFunction:
    def __init__(self, return_value):
        self.return_value = return_value

    def __call__(self, *pos, **kwargs):
        return self.return_value


class MockLongString:
    def __init__(self, string, offset):
        self.string = string
        self.offset = offset

    def __len__(self):
        return len(self.string) + self.offset

    def __getitem__(self, index):
        if not isinstance(index, slice):
            index = slice(index, index + 1)
        index = slice(index.start - self.offset, index.stop - self.offset, index.step)
        if index.start < 0:
            raise NotImplementedError('string portion not given')
        return self.string[index]


def get_example_genes():
    global _EXAMPLE_GENES
    if _EXAMPLE_GENES is None:
        _EXAMPLE_GENES = set_example_genes()
    return _EXAMPLE_GENES


def set_example_genes():
    result = {}
    genes = load_annotations(get_data('example_genes.json'))
    seqs = load_reference_genome(get_data('example_genes.fa'))
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


def flags_from_number(flag: int) -> Dict:
    return dict(
        is_unmapped=bool(flag & int(0x4)),
        mate_is_unmapped=bool(flag & int(0x8)),
        is_reverse=bool(flag & int(0x10)),
        mate_is_reverse=bool(flag & int(0x20)),
        is_read1=bool(flag & int(0x40)),
        is_secondary=bool(flag & int(0x100)),
        is_supplementary=bool(flag & int(0x400)),
    )


class MockString(str):
    def __getitem__(self, key):
        if isinstance(key, slice):
            size = key.stop - key.start
        else:
            size = 1
        return 'N' * size


@dataclass
class MockRead:
    cigarstring: str
    reference_start: int
    reference_name: str = ''
    query_sequence: str = MockString()
    is_reverse: bool = False
    query_name: str = ''
    is_read1: bool = False
    is_supplementary: bool = False
    is_secondary: bool = False
    is_paired: bool = True
    template_length: Optional[int] = None
    is_unmapped: bool = False
    mapping_quality: int = NA_MAPPING_QUALITY

    # mate flags
    next_reference_name: str = ''
    mate_is_reverse: bool = False
    next_reference_start: Optional[int] = None
    mate_is_unmapped: bool = False

    # custom flags for assembly reads
    alignment_rank: Optional[int] = None

    def __hash__(self):
        return hash(SamRead.__hash__(self))

    @property
    def query_length(self):
        return sum(
            [
                size
                for (cigar_state, size) in self.cigar
                if cigar_state in QUERY_ALIGNED_STATES - {CIGAR.H}
            ]
        )

    @property
    def is_proper_pair(self):
        return self.is_paired

    @property
    def reference_id(self):
        try:
            return int(self.reference_name) - 1
        except ValueError:
            return hash(str(self.reference_name))

    @property
    def next_reference_id(self):
        try:
            return int(self.next_reference_name) - 1
        except ValueError:
            return hash(str(self.next_reference_name))

    @property
    def seq(self):
        return self.query_sequence

    @property
    def cigar(self) -> List[Tuple[int, int]]:
        return convert_string_to_cigar(self.cigarstring)

    @cigar.setter
    def cigar(self, cigartuples: List[Tuple[int, int]]):
        self.cigarstring = convert_cigar_to_string(cigartuples)

    @property
    def reference_end(self):
        reference_pos = self.reference_start
        for state, size in self.cigar:
            if state in REFERENCE_ALIGNED_STATES:
                reference_pos += size
        return reference_pos

    @property
    def query_alignment_start(self):
        query_pos = 0
        for state, size in self.cigar:
            if state == CIGAR.H:
                continue
            elif state in QUERY_ALIGNED_STATES - REFERENCE_ALIGNED_STATES:
                query_pos += size
            elif state in QUERY_ALIGNED_STATES:
                return query_pos
        return None

    @property
    def query_alignment_end(self):
        query_pos = 0
        for state, size in self.cigar[::-1]:
            if state == CIGAR.H:
                continue
            elif state in QUERY_ALIGNED_STATES - REFERENCE_ALIGNED_STATES:
                query_pos += size
            elif state in QUERY_ALIGNED_STATES:
                return self.query_length - query_pos
        return None

    @property
    def query_alignment_length(self):
        return self.query_alignment_end - self.query_alignment_start

    @property
    def query_alignment_sequence(self):
        return self.query_sequence[self.query_alignment_start : self.query_alignment_end]

    def has_tag(self, tag_name: str) -> bool:
        return hasattr(self, tag_name)

    def get_tag(self, tag_name: str) -> bool:
        if self.has_tag(tag_name):
            return getattr(self, tag_name)
        return False

    def set_tag(self, tag_name: str, value, value_type='') -> bool:
        setattr(self, tag_name, value)

    # SamRead methods
    def key(self):
        return SamRead.key(self)

    @property
    def flag(self):
        flag = 0
        for flag_value, attr_name in [
            ('REVERSE', 'is_reverse'),
            ('MATE_REVERSE', 'mate_is_reverse'),
            ('MATE_UNMAPPED', 'mate_is_unmapped'),
            ('UNMAPPED', 'is_unmapped'),
            ('SECONDARY', 'is_secondary'),
            ('SUPPLEMENTARY', 'is_supplementary'),
        ]:
            if getattr(self, attr_name):
                flag |= PYSAM_READ_FLAGS[flag_value]

        if self.is_paired:
            if self.is_read1:
                flag |= PYSAM_READ_FLAGS.FIRST_IN_PAIR
            else:
                flag |= PYSAM_READ_FLAGS.LAST_IN_PAIR
        return flag

    @flag.setter
    def flag(self, value):
        for flag_value, attr_name in [
            ('REVERSE', 'is_reverse'),
            ('MATE_REVERSE', 'mate_is_reverse'),
            ('MATE_UNMAPPED', 'mate_is_unmapped'),
            ('UNMAPPED', 'is_unmapped'),
            ('FIRST_IN_PAIR', 'is_read1'),
            ('SECONDARY', 'is_secondary'),
            ('SUPPLEMENTARY', 'is_supplementary'),
        ]:
            if value & PYSAM_READ_FLAGS[flag_value]:
                setattr(self, attr_name, True)
        if value & PYSAM_READ_FLAGS.LAST_IN_PAIR:
            self.is_read1 = False


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


def mock_read_pair(mock1: MockRead, mock2: MockRead, proper_pair=True, reverse_order=False):
    # make sure pair flags are set
    mock1.is_paired = True
    mock2.is_paired = True

    mock1.is_read1 = not reverse_order
    mock2.is_read1 = reverse_order

    if not mock1.is_reverse and proper_pair:
        mock2.is_reverse = True

    if mock1.reference_name != mock2.reference_name:
        mock1.template_length = 0
        mock2.template_length = 0
    mock1.next_reference_name = mock2.reference_name
    mock1.next_reference_start = mock2.reference_start
    mock1.mate_is_reverse = mock2.is_reverse

    if mock1.template_length is None:
        # mock1.template_length = abs(mock1.reference_start - mock1.next_reference_start) + 1
        # if reverse_order:
        #     mock1.template_length *= -1
        mock1.template_length = mock1.next_reference_start - mock1.reference_start + 1

    mock2.next_reference_name = mock1.reference_name
    mock2.next_reference_start = mock1.reference_start
    mock2.mate_is_reverse = mock1.is_reverse
    if mock2.query_name is None:
        mock2.query_name = mock1.query_name
    mock2.template_length = -1 * mock1.template_length
    return mock1, mock2


def build_transcript(
    gene, exons, cds_start, cds_end, domains, strand=None, is_best_transcript=False, name=None
):
    pre_transcript = PreTranscript(
        exons,
        gene=gene,
        strand=strand if strand is not None else gene.get_strand(),
        is_best_transcript=is_best_transcript,
        name=name,
    )
    if gene is not None:
        gene.unspliced_transcripts.append(pre_transcript)

    for spl in pre_transcript.generate_splicing_patterns():
        t = Transcript(pre_transcript, spl)
        pre_transcript.spliced_transcripts.append(t)

        tx = Translation(cds_start, cds_end, t, domains=domains)
        t.translations.append(tx)

    return pre_transcript
