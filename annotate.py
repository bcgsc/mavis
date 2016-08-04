import TSV
import svgwrite

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
        
        self.genomic_start = int(kwargs.pop('genomic_start'))
        self.genomic_end = int(kwargs.pop('genomic_end'))
        self.cds_start = int(kwargs.pop('cds_start'))
        self.cds_end = int(kwargs.pop('cds_end'))

        for e in kwargs.pop('exons', []):
            s, t = e
            if s > t:
                raise AttributeError('exon start/end tuples must have start <= end')
            e = Exon(s, t, transcript = self)
        if self.gene is not None:
            self.gene.transcripts.add(self)
    
    @property
    def key(self):
        return (self.gene, self.name, self.genomic_start, self.genomic_end, self.cds_start, self.cds_end)

class Exon(Bio):
    def __init__(self, start, end, **kwargs):
        self.transcript = kwargs.pop('transcript', None)
        self.start = int(start)
        self.end = int(end)
        self.name = kwargs.pop('name', None)
        if self.transcript is not None:
            self.transcript.exons.add(self)
    
    @property
    def key(self):
        return (self.transcript, self.start, self.end, self.name)

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
        exons = [(int(s), int(t)) for s, t in exons]
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

ref = load_reference('temp')

# '10:61665878' '10:43612031'
pos = 61665878
#pos = 43612031
for gene in ref['10']:
    for t in gene.transcripts:
        if pos < t.genomic_start or pos > t.genomic_end:
            continue
        print('t', t, t.gene.aliases)
