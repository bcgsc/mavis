import TSV
from structural_variant.interval import Interval
from structural_variant.constants import *
from structural_variant.error import *
from Bio import SeqIO

TTRACK = 'transformed_from'


class Bio:

    @property
    def key(self):
        raise NotImplementedError('abstract method must be overidden')

    def __repr__(self):
        cls = self.__class__.__name__
        return cls + str(tuple([k for k in self.key if k is not None]))

    def __init__(self, name=None, reference_object=None):
        self.name = None
        self.reference_object = None


class Gene(Bio):

    def __init__(self, chr=None, name=None, strand=None, aliases=[]):
        Bio.__init__(self, name=name, reference_object=chr)
        self.transcripts = set()
        self.strand = strand
        self.aliases = aliases

        if self.name is None or self.reference_object is None or self.strand is None:
            raise AttributeError('properties: name, reference_object/chr, and strand are required')

    @property
    def chr(self):
        return self.reference_object

    @property
    def key(self):
        return (self.name, self.strand, self.chr)


class Transcript(Bio):

    def __init__(self, gene=None, name=None, cds_start=None, cds_end=None, strand=STRAND.NS, exons=[]):
        """ creates a new transcript object

        Args:
            gene (Gene, optional): the gene this transcript belongs to
            name (str, optional): the name of the transcript or external db id. For example ENTS0001
            cds_start (int): the position (wrt the first exon) where translation would begin
            cds_end (int): the position (wrt the first exon) where translation would terminate
            strand (STRAND, default=STRAND.NS): the strand the transcript occurs on

        """
        Bio.__init__(self, reference_object=gene, name=name)

        self.exons = set()
        self.domains = set()
        self._strand = strand

        try:
            self.cds_start = int(cds_start)
            self.cds_end = int(cds_end)
        except ValueError:
            raise AttributeError('cds_end and/or cds_start are required and must be integers')

        if self.cds_start > self.cds_end:
            raise AttributeError('cds_end must be >= cds_start')

        for e in exons:
            self.exons.add(e)
            e.transcript = self

        if self.gene is not None:
            self.gene.transcripts.add(self)

        exons = sorted(self.exons, key=lambda x: (x.start, x.end))
        for i, exon in enumerate(exons):
            if i == 0:
                continue
            previous = exons[i - 1]
            if Interval(exon.start, exon.end).overlap(Interval(previous.start, previous.end)):
                raise AttributeError(
                    'exons cannot overlap within a transcript')

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
        return sum([e.length() for e in self.exons])

    def _exon_genomic_to_cdna_mapping(self):
        mapping = {}

        exons = self.get_exons()
        if self.strand == STRAND.POS:
            pass
        elif self.strand == STRAND.NEG:
            exons.reverse()
        else:
            raise StrandSpecificityError(
                'cannot convert without strand information')

        l = 0
        for e in exons:
            mapping[(e.start, e.end)] = (l + 1, l + e.length())
            l += e.length()
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
            raise AttributeError(
                'cannot trim without specifying orientation of the breakpoint')
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
            raise AttributeError(
                'cannot trim a transcript with a non-specific exonic breakpoint')
        if exon_num == 0 or exon_num == len(exons):
            raise UserWarning(
                'cannot trim a  transcript if the breakpoint is outside the exons')
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


class Exon(Bio):

    def __init__(self, start, end, transcript=None, name=None):
        Bio.__init__(self, name=name, reference_object=transcript)
        self.start = int(start)
        self.end = int(end)
        self.name = name
        if self.start > self.end:
            raise AttributeError('exon start must be <= exon end')
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

    def length(self):
        return self.end - self.start + 1


class Domain(Bio):

    def __init__(self, name, regions, transcript=None):
        Bio.__init__(self, name=name, reference_object=transcript)
        self.name = name
        self.regions = sorted(list(set(regions)))  # remove duplicates
        for i, region in enumerate(self.regions):
            if region[0] > region[1]:
                raise AttributeError('domain region start must be <= end')
            if i > 0:
                previous_region = self.regions[i - 1]
                if previous_region[1] >= region[0]:
                    raise AttributeError(
                        'regions cannot overlap', previous_region, region)
        if self.transcript is not None:
            self.transcript.domains.add(self)

    @property
    def transcript(self):
        return self.reference_object

    @property
    def key(self):
        return tuple([self.name, self.transcript] + self.regions)

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key


def load_reference_genes(filepath):
    """
    given a file in the std input format (see below) reads and return a list of genes (and sub-objects)

    Expected column headers:
    - ensembl_transcript_id
    - ensembl_gene_id
    - strand
        - [-1] gene is on the reverse strand
        - else: gene is on the forward strand
    - transcript_genomic_start
    - transcript_genomic_end
    - cdna_coding_start
        - position where translation would start, given wrt the transcript start being 0
    - cdna_coding_end
        - position where translation would terminate (last base of last AA), given wrt the transcript start being 0
    - AA_domain_ranges
        - a semi-colon delimited list of domain names and their respective amino acid intervals
        - general pattern
            - <domain name>:<start>-<end>(,<start>-<end>)*(;<domain name>:<start>-<end>(,<start>-<end>)*)*
    - genomic_exon_ranges
        - list of genome start/end positions for exons. semi-colon delimited
        - general pattern
            - <start>-<end>(;<start>-<end>)*
    - hugo_names
        - semi-colon delimited list of hugo names mapped to the ensembl gene id
    """
    header, rows = TSV.read_file(filepath)
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
            genomic_start=row['transcript_genomic_start'],
            genomic_end=row['transcript_genomic_end'],
            cds_start=row['cdna_coding_start'],
            cds_end=row['cdna_coding_end']
        )
        if row['AA_domain_ranges']:
            domains = []
            for d in row['AA_domain_ranges'].split(';'):
                name, temp = d.rsplit(':')
                temp = [x.split('-') for x in temp.split(',')]
                temp = [(int(x), int(y)) for x, y in temp]
                d = Domain(name, temp, transcript=t)
        print(t)
        print(['Exon({0}, {1})'.format(e.start, e.end) for e in exons])

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
            if breakpoint.pos.overlaps((transcript.genomic_start(), transcript.genomic_end())):
                putative_annotations.append(transcript)
    return putative_annotations


def annotations(ref, breakpoint_pairs):  # TODO
    """
    @param ref \a required (type: \b Dict<str,List<Gene>>) the list of reference genes hashed by chromosomes
    @param breakpoint_pairs \a required (type: \b List<BreakpointPair>) breakpoint pairs we wish to annotate as events
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
