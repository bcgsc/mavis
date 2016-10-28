import TSV
from structural_variant.interval import Interval
from structural_variant.constants import *
from structural_variant.error import *
from Bio import SeqIO
import re

TTRACK = 'transformed_from'


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


class Gene(Bio):
    """
    """
    def __init__(self, chr=None, name=None, strand=None, aliases=[]):
        """
        Args:
            chr (str): the chromosome
            name (str): the gene name/id i.e. ENSG0001
            strand (STRAND): the genomic strand '+' or '-'
            aliases (List[str]): a list of aliases. For example the hugo name could go here

        """
        Bio.__init__(self, name=name, reference_object=chr)
        self.transcripts = set()
        self.strand = strand
        self.aliases = aliases

        if self.name is None or self.reference_object is None or self.strand is None:
            raise AttributeError('properties: name, reference_object/chr, and strand are required')

    @property
    def chr(self):
        """returns the name of the chromosome that this gene resides on"""
        return self.reference_object

    @property
    def key(self):
        return (self.name, self.strand, self.chr)


class Transcript(Bio):
    """
    """
    def __init__(self, cds_start, cds_end, gene=None, name=None, strand=None, exons=[]):
        """ creates a new transcript object

        Args:
            gene (Gene, optional): the gene this transcript belongs to
            name (str, optional): the name of the transcript or external db id. For example ENTS0001
            cds_start (int): the position (wrt the first exon) where translation would begin
            cds_end (int): the position (wrt the first exon) where translation would terminate
            strand (STRAND, optional): the strand the transcript occurs on

        """
        Bio.__init__(self, reference_object=gene, name=name)

        self.exons = set()
        self.domains = set()
        self._strand = strand

        if self._strand and self.gene and self.gene.strand != self._strand:
            raise AttributeError('strand does not match reference object')

        try:
            self.cds_start = int(cds_start)
            self.cds_end = int(cds_end)
        except ValueError:
            raise AttributeError('cds_end and/or cds_start must be integers')

        if self.cds_start > self.cds_end:
            raise AttributeError('cds_end must be >= cds_start')

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

    @property
    def gene(self):
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

    def genomic_start(self):
        return min([e.start for e in self.exons])

    def genomic_end(self):
        return max([e.end for e in self.exons])

    def length(self):
        return sum([len(e) for e in self.exons])

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
        return (self.gene, self.name, self.cds_start, self.cds_end)

    def exon_number(self, exon):
        for i, e in enumerate(self.get_exons()):
            if e is exon:
                return i
        raise AttributeError('can only calculate phase on associated exons')

    def trim(self, breakpoint):  # TODO
        if breakpoint.orient == ORIENT.NS:
            raise AttributeError('cannot trim without specifying orientation of the breakpoint')
        # if breakpoint.start != breakpoint.end:
        #    raise AttributeError('cannot trim non-specific breakpoints')
        # t = Transcript(strand = self.gene.strand)
        setattr(t, TTRACK, self)
        # recall that breakpoint orientation is given wrt to the forward strand
        # of the genome
        exon_num = 0
        found_in_previous_intron = False
        found_in_exon = False
        exons = sorted(self.exons, key=lambda x: (x.start, x.end))

        exon_num, found_in_previous_intron = Transcript.position_in_range(
            [(e.start, e.end) for e in exons], breakpoint)

        if not found_in_previous_intron and breakpoint.start != breakpoint.end:
            raise AttributeError('cannot trim a transcript with a non-specific exonic breakpoint')
        if exon_num == 0 or exon_num == len(exons):
            raise UserWarning('cannot trim a  transcript if the breakpoint is outside the exons')
        print(exon_num, found_in_previous_intron)
        new_exons = []
        if breakpoint.orient == ORIENT.LEFT:
            # =====------> OR  <====-------
            for e in exons[:exon_num]:
                temp = Exon(e.start, e.end)
                setattr(temp, TTRACK, e)
                new_exons.append(temp)
            if found_in_exon:
                # create the final truncated exon
                e = exons[exon_num]
                temp = Exon(e.start, breakpoint.start)
                setattr(temp, TTRACK, e)
                new_exons.append(temp)
        else:
            # -----======> OR <----=======
            if found_in_exon:
                # create the first truncated exon
                e = exons[exon_num]
                temp = Exon(breakpoint.start, e.end, transcript=t)
                setattr(temp, TTRACK, e)
                new_exons.append(temp)
            for e in exons[exon_num:]:
                temp = Exon(e.start, e.end, transcript=t)
                setattr(temp, TTRACK, e)
                new_exons.append(temp)
        print(new_exons)

        # translation starts after the last exon
        if self.cds_start > new_exons[-1].end:
            pass
        # translation ends before the first exon
        elif self.cds_end < new_exons[0]:
            pass
        # now pull out all the domains
        new_domains = []
        for domain in self.domains:
            # figure out where the breakpoint lands in the set of ranges
            # find the breakpoint position in AA
            domain_num, in_previous = position_in_range(domain.regions, b)
            for start, end in domain.regions:
                pass
        return t


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
        """
        BioInterval.__init__(self, name=name, reference_object=transcript, start=start, end=end)
        if self.transcript is not None:
            self.transcript.exons.add(self)

    @property
    def transcript(self):
        """alias for the reference_object attribute"""
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
        """
        the transcript this exon is associated with
        """
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
        retain=['chr', 'start', 'end', 'name'],
        cast = {'start': 'int', 'end': 'int'},
        transform = {'chr': lambda x: re.sub('^chr', '', x) }
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
    header, rows = TSV.read_file(
        filepath,
        require=[
            'ensembl_gene_id',
            'strand',
            'hugo_names',
            'chr',
            'genomic_exon_ranges',
            'ensembl_transcript_id',
            'transcript_genomic_start',
            'transcript_genomic_end',
            'cdna_coding_start',
            'cdna_coding_end',
            'AA_domain_ranges'
        ]
    )
    genes = {}
    for row in rows:
        if row['strand'] == '-1':
            row['strand'] = '-'
        else:
            row['strand'] = '+'
        g = Gene(name=row['ensembl_gene_id'], strand=row['strand'], chr=row['chr'],
                 aliases=row['hugo_names'].split(';'))
        if g.name in genes:
            g = genes[g.name]
        else:
            genes[g.name] = g
        exons = [x.split('-') for x in row['genomic_exon_ranges'].split(';')]
        exons = [Exon(int(s), int(t)) for s, t in exons]

        t = Transcript(
            name=row['ensembl_transcript_id'],
            gene=g,
            exons=exons,
            cds_start=row['cdna_coding_start'],
            cds_end=row['cdna_coding_end']
        )
        if row['AA_domain_ranges']:
            domains = []
            for d in row['AA_domain_ranges'].split(';'):
                name, temp = d.rsplit(':')
                temp = temp.split(',')
                temp = [x.split('-') for x in temp]
                temp = [(int(x), int(y)) for x, y in temp]
                d = Domain(name, temp, transcript=t)

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
            if breakpoint.strand != STRAND.NS and transcript.strand != breakpoint.strand:
                continue
            if Interval.overlaps(
                    breakpoint,
                    (transcript.genomic_start(), transcript.genomic_end())):
                putative_annotations.append(transcript)
    return putative_annotations


def annotations(ref, breakpoint_pairs):  # TODO
    """
    Args:
        ref \a required (type: \b Dict<str,List<Gene>>) the list of reference genes hashed by chromosomes
        breakpoint_pairs \a required (type: \b List<BreakpointPair>) breakpoint pairs we wish to annotate as events
    """
    annotations = []

    for bp in breakpoint_pairs:
        putative_break1_transcripts = [None]
        putative_break2_transcripts = [None]

        for gene in ref[bp.break1.chr]:
            for t in gene.transcripts:
                if bp.break1.end < t.genomic_start or bp.break1.start > t.genomic_end:
                    continue  # no overlap
                putative_break1_transcripts.append(t)
        for gene in ref[bp.break2.chr]:
            for t in gene.transcripts:
                if bp.break2.end < t.genomic_start or bp.break2.start > t.genomic_end:
                    continue  # no overlap
                putative_break2_transcripts.append(t)

        temp = itertools.product(
            putative_break1_transcripts, putative_break2_transcripts)
        combinations = []

        # assume that the transcript partners cannot be two different transcripts from the same gene
        # if the transcripts have the same gene they must be the same
        # transcript
        for t1, t2 in temp:
            if t1 is None or t2 is None:
                combinations.append((t1, t2))
                continue
            elif t1.gene == t2.gene and t1 != t2:
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
            combinations.append((t1, t2))

        # now we have a list of putative combinations


def load_reference_genome(filename):
    HUMAN_REFERENCE_GENOME = None
    with open(filename, 'rU') as fh:
        HUMAN_REFERENCE_GENOME = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    return HUMAN_REFERENCE_GENOME

if __name__ == '__main__':
    import doctest
    doctest.testmod()
