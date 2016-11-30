import TSV
from structural_variant.interval import Interval
from structural_variant.constants import *
import itertools
from structural_variant.error import *
from structural_variant.breakpoint import BreakpointPair
from Bio import SeqIO
import re


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

        if not self.breakpoint_pair.interchromosomal:
            try:
                encompassment = Interval(self.breakpoint_pair.break1.end + 1, self.breakpoint_pair.break2.start - 1)
                if gene in encompassment:
                    self.encompassed_genes.add(gene)
            except AttributeError:
                pass
        if Interval.overlaps(gene, self.breakpoint_pair.break1) and gene.chr == self.breakpoint_pair.break1.chr \
                and gene != self.transcript1.reference_object:
            self.genes_at_break1.add(gene)
        if Interval.overlaps(gene, self.breakpoint_pair.break2) and gene.chr == self.breakpoint_pair.break2.chr \
                and gene != self.transcript2.reference_object:
            self.genes_at_break2.add(gene)

        if gene in self.genes_at_break1 or gene in self.genes_at_break2 or gene in self.encompassed_genes \
                or gene == self.transcript1.reference_object or gene == self.transcript2.reference_object:
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


class Bio:
    """
    base class for biological type annotation objects
    """
    @property
    def key(self):
        raise NotImplementedError('abstract method must be overidden')

    def __init__(self, name=None, reference_object=None):
        """
        Args:
            name (string): name of the biological object
            reference_object: the object that is the 'parent' or reference for the current object/interval
        """
        self.name = name
        self.reference_object = reference_object


class BioInterval(Bio, Interval):
    """
    """
    def __init__(self, reference_object, start, end, name=None):
        Bio.__init__(self, reference_object=reference_object, name=name)
        Interval.__init__(self, start, end)


class IntergenicRegion(BioInterval):
    def __init__(self, chr, start, end, strand):
        BioInterval.__init__(self, chr, start, end)
        self.strand = strand


class Gene(BioInterval):
    """
    """
    def __init__(self, chr, start, end, name=None, strand=None, aliases=[]):
        """
        Args:
            chr (str): the chromosome
            name (str): the gene name/id i.e. ENSG0001
            strand (STRAND): the genomic strand '+' or '-'
            aliases (List[str]): a list of aliases. For example the hugo name could go here
        Example:
            >>> Gene('X', 'ENG0001', '+', ['KRAS'])
        """
        BioInterval.__init__(self, name=name, reference_object=chr, start=start, end=end)
        self.transcripts = set()
        self.strand = STRAND.enforce(strand)
        self.aliases = aliases

        if self.name is None or self.reference_object is None or self.strand == STRAND.NS:
            raise AttributeError('properties: name, reference_object/chr, and strand are required')

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    @property
    def key(self):
        return (self.name, self.strand, self.chr, self.start, self.end)


class Transcript(BioInterval):
    """
    """
    def __init__(
        self,
        cds_start=None,
        cds_end=None,
        genomic_start=None,
        genomic_end=None,
        gene=None,
        name=None,
        strand=None,
        exons=[],
        domains=[]
    ):
        """ creates a new transcript object

        Args:
            gene (Gene, optional): the gene this transcript belongs to
            name (str, optional): the name of the transcript or external db id. For example ENTS0001
            cds_start (int): the position (wrt the first exon) where translation would begin
            cds_end (int): the position (wrt the first exon) where translation would terminate
            strand (STRAND, optional): the strand the transcript occurs on

        """
        if genomic_start is None and len(exons) > 0:
            genomic_start = min([e[0] for e in exons])
        if genomic_end is None and len(exons) > 0:
            genomic_end = max([e[1] for e in exons])

        BioInterval.__init__(self, reference_object=gene, name=name, start=genomic_start, end=genomic_end)

        self.exons = set()
        self.domains = set()
        self._strand = strand

        if self._strand and self.gene and self.gene.strand != self._strand:
            raise AttributeError('strand does not match reference object')

        try:
            self.cds_start = int(cds_start) if cds_start is not None else None
            self.cds_end = int(cds_end) if cds_end is not None else None
            if cds_start is not None and cds_end is not None:
                if self.cds_start > self.cds_end:
                    raise AttributeError('cds_end must be >= cds_start')
        except ValueError:
            raise AttributeError('cds_end and/or cds_start must be integers')

        for e in exons:
            if isinstance(e, Exon):
                self.exons.add(e)
                e.reference_object = self
            else:
                Exon(e[0], e[1], transcript=self)

        if self.gene is not None:
            self.gene.transcripts.add(self)

        exons = sorted(self.exons, key=lambda x: (x.start, x.end))
        for i, exon in enumerate(exons):
            if i == 0:
                continue
            previous = exons[i - 1]
            if Interval.overlaps((previous.start, previous.end), (exon.start, exon.end)):
                raise AttributeError('exons cannot overlap within a transcript')

        for d in domains:
            d.reference_object = self
            self.domains.add(d)

    @property
    def gene(self):
        """(:class:`~structural_variant.annotate.Gene`): the gene this transcript belongs to"""
        return self.reference_object

    @property
    def strand(self):
        if self.gene is not None:
            return self.gene.strand
        return self._strand

    def get_exons(self):
        return sorted(self.exons, key=lambda x: x.start)

    def genomic_length(self):
        return self.genomic_end() - self.genomic_start() + 1

    @property
    def genomic_start(self):
        return self.start

    @property
    def genomic_end(self):
        return self.end

    def _exon_genomic_to_cdna_mapping(self):
        mapping = {}

        exons = self.get_exons()
        if self.strand == STRAND.POS:
            pass
        elif self.strand == STRAND.NEG:
            exons.reverse()
        else:
            raise StrandSpecificityError('cannot convert without strand information')

        l = 0
        for e in exons:
            mapping[(e.start, e.end)] = (l + 1, l + len(e))
            l += len(e)
        return mapping

    def _exon_cdna_to_genomic_mapping(self):
        mapping = {}
        for k, v in self._exon_genomic_to_cdna_mapping():
            mapping[v] = k
        return mapping

    def convert_genomic_to_cdna(self, pos):
        mapping = self._exon_genomic_to_cdna_mapping()
        return Interval.convert_pos(mapping, pos)

    @property
    def key(self):
        return (self.gene, self.name, self.start, self.end, self.cds_start, self.cds_end)

    def exon_number(self, exon):
        for i, e in enumerate(self.get_exons()):
            if e is exon:
                return i
        raise AttributeError('can only calculate phase on associated exons')


class Exon(BioInterval):
    """
    """
    def __init__(self, start, end, transcript=None, name=None):
        """
        Args:
            start (int): the genomic start position
            end (int): the genomic end position
            name (str, optional): the name of the exon
            transcript (Transcript, optional): the 'parent' transcript this exon belongs to
        Raises:
            AttributeError: if the exon start > the exon end
        Example:
            >>> Exon(15, 78)
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end)
        if self.transcript is not None:
            self.transcript.exons.add(self)

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.Transcript`): the transcript this exon belongs to"""
        return self.reference_object

    def __getitem__(self, index):
        try:
            index = int(index)
        except ValueError:
            raise IndexError('indices must be integers', index)
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        raise IndexError('index out of bounds', index)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key

    @property
    def key(self):
        return (self.transcript, self.start, self.end, self.name)


class Domain(Bio):
    """
    """
    def __init__(self, name, regions, transcript=None):
        """
        Args:
            name (str): the name of the domain i.e. PF00876
            regions (List[Interval]): the amino acid ranges that are part of the domain
            transcript (Transcript, optional): the 'parent' transcript this domain belongs to
        Raises:
            AttributeError: if the end of any region is less than the start
        Example:
            >>> Domain('DNA binding domain', [(1, 4), (10, 24)], transcript)
        """
        Bio.__init__(self, name=name, reference_object=transcript)
        self.name = name
        self.regions = sorted(list(set(regions)))  # remove duplicates

        for i, region in enumerate(self.regions):
            if region[0] > region[1]:
                raise AttributeError('domain region start must be <= end')
        self.regions = Interval.min_nonoverlapping(*self.regions)
        if self.transcript is not None:
            self.transcript.domains.add(self)

    @property
    def transcript(self):
        """(:class:`~structural_variant.annotate.Transcript`): the transcript this domain belongs to"""
        return self.reference_object

    @property
    def key(self):
        return tuple([self.name, self.transcript] + self.regions)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key


def load_masking_regions(filepath):
    """
    reads a file of regions. The expect input format for the file is tab-delimited and

    +---------------+---------------+-----------------------+
    | column name   | value type    | description           |
    +===============+===============+=======================+
    | chr           | string        | the chromosome name   |
    +---------------+---------------+-----------------------+
    | start         | int           | the start position    |
    +---------------+---------------+-----------------------+
    | end           | int           | the end position      |
    +---------------+---------------+-----------------------+
    | name          | string        | label for the region  |
    +---------------+---------------+-----------------------+

    Args:
        filepath (string): path to the input tab-delimited file
    Returns:
        Dict[str,List[BioInterval]]:
            a dictionary keyed by chromosome name with values of lists of regions on the chromosome

    """
    header, rows = TSV.read_file(
        filepath,
        require=['chr', 'start', 'end', 'name'],
        cast={'start': int, 'end': int, 'chr': lambda x: re.sub('^chr', '', x)}
    )
    regions = {}
    for row in rows:
        r = BioInterval(reference_object=row['chr'], start=row['start'], end=row['end'], name=row['name'])
        regions.setdefault(r.reference_object, []).append(r)
    return regions


def load_reference_genes(filepath):
    """
    given a file in the std input format (see below) reads and return a list of genes (and sub-objects)

    +-----------------------+---------------------------+-----------------------------------------------------------+
    | column name           | example                   | description                                               |
    +=======================+===========================+===========================================================+
    | ensembl_transcript_id | ENST000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | ensembl_gene_id       | ENSG000001                |                                                           |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | strand                | -1                        | positive or negative 1                                    |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_start     | 44                        | where translation begins relative to the start of the cdna|
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | cdna_coding_end       | 150                       | where translation terminates                              |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | genomic_exon_ranges   | 100-201;334-412;779-830   | semi-colon demitited exon start/ends                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | AA_domain_ranges      | DBD:220-251,260-271       | semi-colon delimited list of domains                      |
    +-----------------------+---------------------------+-----------------------------------------------------------+
    | hugo_names            | KRAS                      | hugo gene name                                            |
    +-----------------------+---------------------------+-----------------------------------------------------------+

    Args:
        filepath (string): path to the input tab-delimited file

    Returns:
        Dict[str,List[Gene]]: a dictionary keyed by chromosome name with values of list of genes on the chromosome
    """
    def parse_exon_list(row):
        exons = []
        for temp in row.split(';'):
            try:
                s, t = temp.split('-')
                exons.append(Exon(int(s), int(t)))
            except:
                print('exon error:', temp, row)
        return exons

    def parse_domain_list(row):
        domains = []
        for d in row.split(';'):
            try:
                name, temp = d.rsplit(':')
                temp = temp.split(',')
                temp = [x.split('-') for x in temp]
                temp = [(int(x), int(y)) for x, y in temp]
                d = Domain(name, temp)
            except:
                print('error in d:', d, row)
        return domains

    def nullable_int(row):
        try:
            row = int(row)
        except ValueError:
            row = TSV.null(row)
        return row

    def parse_strand(row):
        if row in ['-1', '-']:
            return STRAND.NEG
        elif row in ['1', '+', '+1']:
            return STRAND.POS
        raise ValueError('cast to strand failed')

    header, rows = TSV.read_file(
        filepath,
        require=[
            'ensembl_gene_id',
            'chr',
            'ensembl_transcript_id'
        ],
        add={
            'cdna_coding_start': '',
            'cdna_coding_end': '',
            'AA_domain_ranges': '',
            'genomic_exon_ranges': '',
            'hugo_names': '',
            'transcript_genomic_start': '',
            'transcript_genomic_end': ''
        },
        cast={
            'genomic_exon_ranges': parse_exon_list,
            'AA_domain_ranges': parse_domain_list,
            'cdna_coding_end': nullable_int,
            'cdna_coding_start': nullable_int,
            'transcript_genomic_end': nullable_int,
            'transcript_genomic_start': nullable_int,
            'strand': parse_strand,
            'gene_start': int,
            'gene_end': int
        }
    )
    genes = {}
    for row in rows:
        g = Gene(
            row['chr'],
            row['gene_start'],
            row['gene_end'],
            name=row['ensembl_gene_id'],
            strand=row['strand'],
            aliases=row['hugo_names'].split(';')
        )
        if g.name in genes:
            g = genes[g.name]
        else:
            genes[g.name] = g

        t = Transcript(
            name=row['ensembl_transcript_id'],
            gene=g,
            genomic_start=row['transcript_genomic_start'],
            genomic_end=row['transcript_genomic_end'],
            exons=row['genomic_exon_ranges'],
            domains=row['AA_domain_ranges'],
            cds_start=row['cdna_coding_start'],
            cds_end=row['cdna_coding_end']
        )

    ref = {}
    for g in genes.values():
        if g.chr not in ref:
            ref[g.chr] = []
        ref[g.chr].append(g)
    return ref


def overlapping_transcripts(ref_ann, breakpoint):
    """
    Args:
        ref_ann (Dict[str,List[Gene]]): the reference list of genes split by chromosome
        breakpoint (Breakpoint): the breakpoint in question
    Returns:
        List[Transcript]: a list of possible transcripts
    """
    putative_annotations = []
    for gene in ref_ann[breakpoint.chr]:
        for transcript in gene.transcripts:
            if breakpoint.strand != STRAND.NS and transcript.strand != STRAND.NS \
                    and transcript.strand != breakpoint.strand:
                continue
            if Interval.overlaps(breakpoint, transcript):
                putative_annotations.append(transcript)
    return putative_annotations


def gather_breakpoint_annotations(ref_ann, breakpoint):

    pos_overlapping_transcripts = []
    neg_overlapping_transcripts = []
    for gene in ref_ann[breakpoint.chr]:
        for t in gene.transcripts:
            if Interval.overlaps(t, breakpoint):
                if STRAND.compare(t.strand, STRAND.POS):
                    pos_overlapping_transcripts.append(t)
                if STRAND.compare(t.strand, STRAND.NEG):
                    neg_overlapping_transcripts.append(t)

    pos_intervals = Interval.min_nonoverlapping(*pos_overlapping_transcripts)
    neg_intervals = Interval.min_nonoverlapping(*neg_overlapping_transcripts)
    
    temp = []
    # before the first?
    if len(pos_intervals) > 0:
        first = pos_intervals[0]
        last = pos_intervals[-1]
        if breakpoint < first:
            temp.append(IntergenicRegion(breakpoint.chr, breakpoint[0], first[0] - 1, STRAND.POS))
        if breakpoint[1] > last[1]:
            temp.append(IntergenicRegion(breakpoint.chr, last[1] + 1, breakpoint[1], STRAND.POS))

        for i, curr in enumerate(pos_intervals):
            if i > 0:
                prev = pos_intervals[i - 1]
                try:
                    temp.append(IntergenicRegion(breakpoint.chr, prev[1] + 1, curr[0] - 1, STRAND.POS))
                except AttributeError:
                    pass
    else:
        temp.append(IntergenicRegion(breakpoint.chr, breakpoint.start, breakpoint.end, STRAND.POS))
    pos_overlapping_transcripts.extend(temp)

    temp = []
    # before the first?
    if len(neg_intervals) > 0:
        first = neg_intervals[0]
        last = neg_intervals[-1]
        if breakpoint < first:
            temp.append(IntergenicRegion(breakpoint.chr, breakpoint[0], first[0] - 1, STRAND.NEG))
        if breakpoint[1] > last[1]:
            temp.append(IntergenicRegion(breakpoint.chr, last[1] + 1, breakpoint[1], STRAND.NEG))

        for i, curr in enumerate(neg_intervals):
            if i > 0:
                prev = neg_intervals[i - 1]
                try:
                    temp.append(IntergenicRegion(breakpoint.chr, prev[1] + 1, curr[0] - 1, STRAND.NEG))
                except AttributeError:
                    pass
    else:
        temp.append(IntergenicRegion(breakpoint.chr, breakpoint.start, breakpoint.end, STRAND.NEG))
    neg_overlapping_transcripts.extend(temp)

    return sorted(pos_overlapping_transcripts), sorted(neg_overlapping_transcripts)


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

    break1_pos, break1_neg = gather_breakpoint_annotations(ref, bp.break1)
    break2_pos, break2_neg = gather_breakpoint_annotations(ref, bp.break2)
    
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
        bpp.break1.start = b1_itvl[0]
        bpp.break1.end = b1_itvl[1]
        bpp.break2.start = b2_itvl[0]
        bpp.break2.end = b2_itvl[1]

        a = Annotation(bpp, a1, a2)

        for gene in ref[bp.break1.chr]:
            a.add_gene(gene)
        if bp.interchromosomal:
            for gene in ref[bp.break2.chr]:
                a.add_gene(gene)
        annotations.append(a)
    return annotations


def load_reference_genome(filename):
    HUMAN_REFERENCE_GENOME = None
    with open(filename, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    return HUMAN_REFERENCE_GENOME

if __name__ == '__main__':
    import doctest
    doctest.testmod()
