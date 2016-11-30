
"""

"""
import itertools
# import TSV
from structural_variant.constants import STRAND
from structural_variant.interval import Interval
from structural_variant.breakpoint import BreakpointPair
from structural_variant import annotate


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

    def add_gene(self, gene):
        if gene.chr not in [self.breakpoint_pair.break1.chr, self.breakpoint_pair.break2.chr]:
            raise AttributeError('cannot add gene not on the same chromosome as either breakpoint')

        if not self.bpp.interchromosomal:
            try:
                encompassment = Interval(self.breakpoint_pair.break1.end + 1, self.breakpoint_pair.break2.start - 1)
                if gene in encompassment:
                    self.encompassed_genes.add(gene)
            except AttributeError:
                pass
        if Interval.overlaps(gene, self.breakpoint_pair.break1) and gene.chr == self.breakpoint_pair.break1.chr:
            self.genes_at_break1.add(gene)
        if Interval.overlaps(gene, self.breakpoint_pair.break2) and gene.chr == self.breakpoint_pair.break2.chr:
            self.genes_at_break2.add(gene)

        if gene in self.genes_at_break1 or gene in self.genes_at_break2 or gene in self.encompassed_genes:
            return

        d1 = Interval.dist(gene, self.breakpoint_pair.break1)
        d2 = Interval.dist(gene, self.breakpoint_pair.break2)

        if self.breakpoint_pair.interchromosomal:
            if gene.chr == self.breakpoint_pair.break1.chr:
                self.nearest_gene_break1.add((gene, d1))
            elif gene.chr == self.breakpoint_pair.break2.chr:
                self.nearest_gene_break2.add((gene, d2))
        else:
            if d1 < 0:
                self.nearest_gene_break1.add((gene, d1))
            if d2 > 0:
                self.nearest_gene_break2.add((gene, d2))

        temp = set()

        tmin = [d for g, d in self.nearest_gene_break1 if d < 0]
        tmax = [d for g, d in self.nearest_gene_break1 if d > 0]
        tmin = 0 if len(tmin) == 0 else max(tmin)
        tmax = 0 if len(tmax) == 0 else min(tmax)

        for gene, dist in self.nearest_gene_break1:
            if tmin != 0 and dist == tmin:
                temp.add((gene, dist))
            elif tmax != 0 and dist == tmax:
                temp.add((gene, dist))

        self.nearest_gene_break1 = temp

        temp = set()

        tmin = [d for g, d in self.nearest_gene_break2 if d < 0]
        tmax = [d for g, d in self.nearest_gene_break2 if d > 0]
        tmin = 0 if len(tmin) == 0 else max(tmin)
        tmax = 0 if len(tmax) == 0 else min(tmax)

        for gene, dist in self.nearest_gene_break2:
            if tmin != 0 and dist == tmin:
                temp.add((gene, dist))
            elif tmax != 0 and dist == tmax:
                temp.add((gene, dist))

        self.nearest_gene_break2 = temp


def read_validation_file(self):
    pass


def gather_annotations(ref, bp):  # TODO
    """
    Args:
        ref (Dict[str,List[Gene]]): the list of reference genes hashed by chromosomes
        breakpoint_pairs (List[BreakpointPair]): breakpoint pairs we wish to annotate as events

    each annotation is defined by the annotations selected at the breakpoints
    the other annotations are given relative to this
    the annotation at the breakpoint can be a transcript or an intergenic region

    """
    annotations = []

    break1_pos, break1_neg = annotate.gather_breakpoint_annotations(ref, bp.break1)
    break2_pos, break2_neg = annotate.gather_breakpoint_annotations(ref, bp.break2)

    combinations = []

    if bp.stranded:
        if bp.break1.strand == STRAND.POS:
            if bp.break1.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_pos, break2_pos))
            else:
                combinations.extend(itertools.product(break1_pos, break2_neg))
        else:
            if bp.break1.strand == STRAND.POS:
                combinations.extend(itertools.product(break1_neg, break2_pos))
            else:
                combinations.extend(itertools.product(break1_neg, break2_neg))
    else:
        if bp.opposing_strands:
            combinations.extend(itertools.product(break1_pos, break2_neg))
            combinations.extend(itertools.product(break1_neg, break2_pos))
        else:
            combinations.extend(itertools.product(break1_pos, break2_pos))
            combinations.extend(itertools.product(break1_neg, break2_neg))

    for a1, a2 in combinations:
        b1_itvl = bp.break1 & a1
        b2_itvl = bp.break2 & a2

        bpp = BreakpointPair.copy(bp)
        bp.break1.start = b1_itvl[0]
        bp.break1.end = b1_itvl[1]
        bp.break2.start = b2_itvl[0]
        bp.break2.end = b2_itvl[1]

        a = Annotation(bpp, a1, a2)

        for gene in ref[bp.break1.chr]:
            a.add_gene(gene)
        if bp.interchromosomal:
            for gene in ref[bp.break2.chr]:
                a.add_gene(gene)
        annotations.add(a)
    return annotations


def main():
    pass
    # load the reference genes
    # grab putative transcripts for fusions


if __name__ == '__main__':
    main()
