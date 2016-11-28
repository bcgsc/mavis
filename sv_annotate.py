
"""

"""
import itertools
import TSV
from structural_variant.constants import STRAND, SVTYPE
from structural_variant.interval import Interval
from structural_variant.breakpoint import BreakpointPair


class Annotation:
    """
    a fusion of two transcripts created by the associated breakpoint_pair
    will also hold the other annotations for overlapping and encompassed and nearest genes
    """
    def __init__(self, bpp, transcript1=None, transcript2=None, data={}):
        self.breakpoint_pair = bpp.copy()
        if transcript1 is not None:
            temp = bpp.break1 & transcript1
            self.breakpoint_pair.break1.start = temp[0]
            self.breakpoint_pair.break1.end = temp[1]
        if transcript2 is not None:
            temp = bpp.break2 & transcript2
            self.breakpoint_pair.break2.start = temp[0]
            self.breakpoint_pair.break2.end = temp[1]
        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.data = {}
        self.data.update(data)

        self.encompassed_genes = set()
        self.nearest_gene_break1 = set()
        self.nearest_gene_break2 = set()
        self.genes_at_break1 = set()
        self.genes_at_break2 = set()


def read_validation_file(self):
    pass

def gather_annotations(ref, bp):  # TODO
    """
    Args:
        ref (Dict[str,List[Gene]]): the list of reference genes hashed by chromosomes
        breakpoint_pairs (List[BreakpointPair]): breakpoint pairs we wish to annotate as events
    """
    annotations = []

    putative_break1_transcripts = set()
    putative_break2_transcripts = set()

    for gene in ref[bp.break1.chr]:
        if bp.stranded and gene.strand != bp.break1.strand:
            continue
        for t in gene.transcripts:
            if bp.break1.end < t.genomic_start or bp.break1.start > t.genomic_end:
                continue  # no overlap
            putative_break1_transcripts.add(t)
    for gene in ref[bp.break2.chr]:
        if bp.stranded and gene.strand != bp.break2.strand:
            continue
        for t in gene.transcripts:
            if bp.break2.end < t.genomic_start or bp.break2.start > t.genomic_end:
                continue  # no overlap
            putative_break2_transcripts.add(t)
    
    temp = itertools.product(
        putative_break1_transcripts, putative_break2_transcripts)

    # assume that the transcript partners cannot be two different transcripts from the same gene
    # if the transcripts have the same gene they must be the same
    # transcript
    # we want t
    combinations = []
    for t1, t2 in temp:
        if t1 is not None and t2 is not None:
            if t1.gene == t2.gene and t1 != t2:
                continue
            if bp.opposing_strands is not None \
                    and bp.opposing_strands != (t1.gene.strand != t2.gene.strand): \
                    # if the stand combination does not match then ignore
                continue
            if bp.stranded:
                if bp.break1.strand != STRAND.NS:
                    if t1.gene.strand != bp.break1.strand:
                        continue
                if bp.break2.strand != STRAND.NS:
                    if t2.gene.strand != bp.break2.strand:
                        continue
        combinations.append(t1, t2)

    for t1, t2 in combinations:
        a = Annotation(bp, t1, t2)

        break1 = bp.break1
        if t1 is not None:
            break1 = break1 & Interval(t1.genomic_start, t1.genomic_end)
        break2 = bp.break2
        if t2 is not None:
            break2 = break2 & Interval(t2.genomic_start, t2.genomic_end)

        encompass = None
        if not bp.interchromosomal:
            try:
                encompass = Interval(break1[1] + 1, break2[0] - 1)
            except AttributeError:
                pass
        b1_dist = []
        b2_dist = []

        for gene in ref[bp.break1.chr]:
            if Interval.overlaps(gene, break1):
                a.genes_at_break1.add(gene)
            if encompass is not None and gene in encompass:
                a.encompassed_genes.add(gene)
            d = Interval.dist(break1, gene)
            b1_dist.append((gene, abs(d)))

        for gene in ref[bp.break2.chr]:
            if Interval.overlaps(gene, break2):
                a.genes_at_break2.add(gene)
            if encompass is not None and gene in encompass:
                a.encompassed_genes.add(gene)
            d = Interval.dist(break2, gene)
            b2_dist.append((gene, abs(d)))

        b1_dist = [
            (g, d) for g, d in b1_dist if gene not in a.genes_at_break1
            and gene not in a.genes_at_break2
            and gene not in a.encompassed_genes
        ]
        b2_dist = [
            (g, d) for g, d in b2_dist if gene not in a.genes_at_break1
            and gene not in a.genes_at_break2
            and gene not in a.encompassed_genes
        ]
        if len(b1_dist) > 0:
            mind = min(b1_dist, key=lambda x: x[1])
            for gene, d in b1_dist:
                if d == mind:
                    a.nearest_gene_break1.add((gene, d))
        if len(b2_dist) > 0:
            mind = min(b2_dist, key=lambda x: x[1])
            for gene, d in b2_dist:
                if d == mind:
                    a.nearest_gene_break2.add((gene, d))
        annotations.append(a)
    return annotations




def main():
    pass
    # load the reference genes
    # grab putative transcripts for fusions


if __name__ == '__main__':
    main()
