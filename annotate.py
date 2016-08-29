import TSV
import svgwrite
from sv import Interval

class Bio:
    @property
    def key(self):
        raise NotImplementedError('abstract method must be overidden')

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key

    def __repr__(self):
        cls = self.__class__.__name__
        return cls + str(self.key)

class Gene(Bio):
    def __init__(self, name, **kwargs):
        self.name = name
        self.transcripts = set()
        self.strand = kwargs.pop('strand')
        self.chr = kwargs.pop('chr')
        self.aliases = kwargs.pop('aliases', [])
    
    @property
    def key(self):
        return (self.name, self.strand, self.chr)

class Transcript(Bio):
    def __init__(self, name, **kwargs):
        self.gene = kwargs.pop('gene', None)
        self.name = name
        
        self.exons = set()
        self.domains = set()
        
        #self.genomic_start = int(kwargs.pop('genomic_start'))
        #self.genomic_end = int(kwargs.pop('genomic_end'))
        self.cds_start = int(kwargs.pop('cds_start')) # genomic position where coding begins
        self.cds_end = int(kwargs.pop('cds_end')) # genomic position where coding ends

        if self.cds_start > self.cds_end:
            raise AttributeError('cds_end must be >= cds_start')

        for e in kwargs.pop('exons', []):
            self.exons.add(e)
            e.transcript = self

        if self.gene is not None:
            self.gene.transcripts.add(self)
        else:
            self._strand = kwargs.pop('strand', STRAND.NS)

        exons = sorted(self.exons, key=lambda x: (x.start, x.end))
        for i, exon in enumerate(exons):
            if i == 0:
                continue
            previous = exons[i - 1]
            if Interval(exon.start, exon.end).overlap(Interval(previous.start, previous.end)):
                raise AttributeError('exons cannot overlap within a transcript')
    
    @property
    def get_exons(self):
        return sorted(self.exons, key=lambda x: x.start)

    @property
    def genomic_length(self):
        return self.genomic_end - self.genomic_start + 1
    
    @property
    def genomic_start(self):
        return min([e.start for e in self.exons])

    @property
    def genomic_end(self):
        return max([e.end for e in self.exons])
    
    @property
    def length(self):
        return sum([e.length for e in self.exons])

    @property
    def strand(self):
        if self.gene is not None:
            return self.gene.strand
        return self._strand

    @property
    def key(self):
        return (self.gene, self.name, self.genomic_start, self.genomic_end)
    
    def exon_number(self, exon):
        for i, e in enumerate(self.get_exons()):
            if e is exon:
                return i
        raise AttributeError('can only calculate phase on associated exons')

    def start_phase(self, exon):
        if exon not in self.exons:
            raise AttributeError('can only calculate phase on associated exons')
        
        if self.strand == STRAND.POS:
            if exon.start < self.cds_start or exon.start > self.cds_end:
                return PHASE.NA
            cds = 0
            for e in sorted(self.exons, key=lambda x: x.start):
                if e is exon:
                    return (cds + 1) % 3
                if self.cds_start <= e.start and self.cds_end >= e.end:
                    # covers the entire exon
                    cds += e.length
                elif self.cds_start >= e.start and self.cds_end <= e.end:
                    # cds is entirely within the current exon
                    cds += self.cds_end - self.cds_start + 1
                elif self.cds_start >= e.start and self.cds_start <= e.end:
                    # covers the end of the exon
                    cds += e.end - self.cds_start + 1
                elif self.cds_end <= e.end and self.cds_end >= e.start:
                    # covers the start of the exon
                    cds += self.cds_end - e.start + 1
        elif self.strand == STRAND.NEG:
            cds = 0
            for e in sorted(self.exons, key=lambda x: x.start, reverse=True):
                if self.cds_start <= e.start and self.cds_end >= e.end:
                    # covers the entire exon
                    cds += e.length
                elif self.cds_start >= e.start and self.cds_end <= e.end:
                    # cds is entirely within the current exon
                    cds += self.cds_end - self.cds_start + 1
                elif self.cds_start >= e.start and self.cds_start <= e.end:
                    # covers the end of the exon
                    cds += e.end - self.cds_start + 1
                elif self.cds_end <= e.end and self.cds_end >= e.start:
                    # covers the start of the exon
                    cds += self.cds_end - e.start + 1
        else:
            raise AttributeError('cannot calculate phase of exons when the strand of the transcript is not specified')
        raise AttributeError('can only calculate phase on associated exons')

    def trim(self, breakpoint):
        if breakpoint.orient == ORIENT.NS:
            raise AttributeError('cannot trim without specifying orientation of the breakpoint')
        #if breakpoint.start != breakpoint.end:
        #    raise AttributeError('cannot trim non-specific breakpoints')
        t = Transcript(strand = self.gene.strand)
        # recall that breakpoint orientation is given wrt to the forward strand of the genome
        exon_num = 0
        found_in_previous_intron = False
        found_in_exon = False
        exons = sorted(self.exons, key=lambda x: (x.start, x.end))
        
        while exon_num < len(self.exons):
            current_exon = exons[exon_num]

            if breakpoint.start >= current_exon.start \
                    and breakpoint.end <= current_exon.end \
                    and breakpoint.start == breakpoint.end:
                # breakpoint range is fully contained in the current exon and specific
                found_in_exon = True
                break
            elif Interval(breakpoint.start, breakpoint.end).overlap(Interval(current_exon.start, current_exon.end)):
                # the breakpoint range overlaps the current exon but is non-specific
                raise AttributeError('cannot trim a transcript with a non-specific exonic breakpoint')
            elif exon_num == 0: # first exon
                if breakpoint.end < current_exon.start:
                    # before the first exon
                    raise AttributeError('cannot trim a transcript where the breakpoint is outside the transcript')
            else:
                # check the previous exon
                previous_exon = exons[exon_num - 1]
                if breakpoint.start > previous_exon.end and breakpoint.end < current_exon.start:
                    found_in_previous_intron = True
                    break
            exon_num += 1
        if not found_in_exon and not found_in_previous_intron:
            # after the last exon
            raise AttributeError('cannot trim a transcript where the breakpoint is outside the transcript')
        if breakpoint.orient == ORIENT.LEFT:
            # =====------> OR  <====-------
            exons = []
            for e in exons[:exon_num]:
                temp = Exon(e.start, e.end)
                setattr(temp, 'transformed_from', e)
                exons.append(temp)
            if found_in_exon:
                # create the final truncated exon
                e = exons[exon_num]
                temp = Exon(e.start, breakpoint.start, transcript=t)
                setattr(temp, 'transformed_from', e)
        else:
            # -----======> OR <----=======
            if found_in_exon:
                # create the final truncated exon
                e = exons[exon_num]
                temp = Exon(breakpoint.start, e.end, transcript=t)
                setattr(temp, 'transformed_from', e)
            for e in exons[exon_num:]:
                temp = Exon(e.start, e.end, transcript=t)
                setattr(temp, 'transformed_from', e)
        return t


class Exon(Bio):
    def __init__(self, start, end, **kwargs):
        self.transcript = kwargs.pop('transcript', None)
        self.start = int(start)
        self.end = int(end)
        self.name = kwargs.pop('name', None)
        if self.start > self.end:
            raise AttributeError('exon start must be <= exon end')
        if self.transcript is not None:
            self.transcript.exons.add(self)
    
    @property
    def key(self):
        return (self.transcript, self.start, self.end, self.name)

    @property
    def length(self):
        return self.end - self.start + 1

class Domain(Bio):
    def __init__(self, name, regions, **kwargs):
        self.name = name
        self.transcript = kwargs.pop('transcript', None)
        self.regions = sorted(list(set(regions))) # remove duplicates
        for s, t in regions:
            if s > t:
                raise AttributeError('domain region start must be <= end')
        if self.transcript is not None:
            self.transcript.domains.add(self)
    
    @property
    def key(self):
        return tuple([self.name, self.transcript] + self.regions)

def load_reference(filepath):
    header, rows = TSV.read_file(filepath)
    genes = {}
    for row in rows:
        if row['strand'] == '-1':
            row['strand'] = '-'
        else:
            row['strand'] = '+'
        g = Gene(row['ensembl_gene_id'], strand = row['strand'], chr = row['chr'],
                aliases = row['hugo_names'].split(';'))
        if g.name in genes:
            g = genes[g.name]
        else:
            genes[g.name] = g
        exons = [ x.split('-') for x in row['genomic_exon_ranges'].split(';')]
        exons = [Exon(int(s), int(t)) for s, t in exons]

        t = Transcript(
                row['ensembl_transcript_id'], 
                gene = g, 
                exons = exons,
                genomic_start = row['transcript_genomic_start'],
                genomic_end = row['transcript_genomic_end'],
                cds_start = row['cdna_coding_start'],
                cds_end = row['cdna_coding_end']
                )
        if row['AA_domain_ranges']:
            domains = []
            for d in row['AA_domain_ranges'].split(';'):
                name, temp = d.rsplit(':')
                temp = [ x.split('-') for x in temp.split(',') ]
                temp = [ (int(x), int(y)) for x, y in temp ]
                d = Domain(name, temp, transcript = t)
    
    ref = {}
    for g in genes.values():
        if g.chr not in ref:
            ref[g.chr] = []
        ref[g.chr].append(g)
    return ref

def annotations(ref, breakpoint_pairs):
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
                    continue # no overlap
                putative_break1_transcripts.append(t)
        for gene in ref[bp.break2.chr]:
            for t in gene.transcripts:
                if bp.break2.end < t.genomic_start or bp.break2.start > t.genomic_end:
                    continue # no overlap
                putative_break2_transcripts.append(t)
        
        temp = itertools.product(putative_break1_transcripts, putative_break2_transcripts)
        combinations = []
        
        # assume that the transcript partners cannot be two different transcripts from the same gene
        # if the transcripts have the same gene they must be the same transcript
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
ref = load_reference('/home/creisle/svn/ensembl_flatfiles/ensembl69_transcript_exons_and_domains_20160808.tsv')

# '10:61665878' '10:43612031'
pos = 61665878
#pos = 43612031
for gene in ref['10']:
    for t in gene.transcripts:
        if pos < t.genomic_start or pos > t.genomic_end:
            continue
        print('t', t, t.gene.aliases)
