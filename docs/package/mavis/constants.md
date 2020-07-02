# mavis.constants

module responsible for small utility functions and constants used throughout the structural_variant package

## PROGNAME

```python
PROGNAME = 'mavis'
```

## EXIT_OK

```python
EXIT_OK = 0
```

## EXIT_ERROR

```python
EXIT_ERROR = 1
```

## EXIT_INCOMPLETE

```python
EXIT_INCOMPLETE = 2
```

## COMPLETE_STAMP

```python
COMPLETE_STAMP = 'MAVIS.COMPLETE'
```

## SUBCOMMAND

```python
SUBCOMMAND = MavisNamespace(
    ANNOTATE='annotate',
    VALIDATE='validate',
    SETUP='setup',
    SCHEDULE='schedule',
    CLUSTER='cluster',
    PAIR='pairing',
    SUMMARY='summary',
    CONFIG='config',
    CONVERT='convert',
    OVERLAY='overlay',
)
```

## CODON_SIZE

```python
CODON_SIZE = 3
```

## GAP

```python
GAP = '-'
```

## ORIENT

```python
ORIENT = MavisNamespace(LEFT='L', RIGHT='R', NS='?')
```

## PROTOCOL

```python
PROTOCOL = MavisNamespace(GENOME='genome', TRANS='transcriptome')
```

## DISEASE_STATUS

```python
DISEASE_STATUS = MavisNamespace(DISEASED='diseased', NORMAL='normal')
```

## STRAND

```python
STRAND = MavisNamespace(POS='+', NEG='-', NS='?')
```

## SVTYPE

```python
SVTYPE = MavisNamespace(
    DEL='deletion',
    TRANS='translocation',
    ITRANS='inverted translocation',
    INV='inversion',
    INS='insertion',
    DUP='duplication',
)
```

## CIGAR

```python
CIGAR = MavisNamespace(M=0, I=1, D=2, N=3, S=4, H=5, P=6, X=8, EQ=7)  # noqa
""":class:`MavisNamespace`: Enum-like. For readable cigar values
```

## NA_MAPPING_QUALITY

```python
NA_MAPPING_QUALITY = 255
```

## PYSAM_READ_FLAGS

```python
PYSAM_READ_FLAGS = MavisNamespace(
    REVERSE=16,
    MATE_REVERSE=32,
    UNMAPPED=4,
    MATE_UNMAPPED=8,
    FIRST_IN_PAIR=64,
    LAST_IN_PAIR=128,
    SECONDARY=256,
    MULTIMAP=1,
    SUPPLEMENTARY=2048,
    TARGETED_ALIGNMENT='ta',
    RECOMPUTED_CIGAR='rc',
    BLAT_RANK='br',
    BLAT_SCORE='bs',
    BLAT_ALIGNMENTS='ba',
    BLAT_PERCENT_IDENTITY='bi',
    BLAT_PMS='bp',
)
```

## DNA_ALPHABET

```python
DNA_ALPHABET = alphabet = Gapped(ambiguous_dna, '-')
```

## alphabet

```python
DNA_ALPHABET = alphabet = Gapped(ambiguous_dna, '-')
```

## DNA_ALPHABET.match

```python
DNA_ALPHABET.match = lambda x, y: _match_ambiguous_dna(x, y)
```

## FLAGS

```python
FLAGS = MavisNamespace(LQ='LOWQUAL')
```

## READ_PAIR_TYPE

```python
READ_PAIR_TYPE = MavisNamespace(RR='RR', LL='LL', RL='RL', LR='LR')
```

## CALL_METHOD

```python
CALL_METHOD = MavisNamespace(
    CONTIG='contig',
    SPLIT='split reads',
    FLANK='flanking reads',
    SPAN='spanning reads',
    INPUT='input',
)
```

## GENE_PRODUCT_TYPE

```python
GENE_PRODUCT_TYPE = MavisNamespace(SENSE='sense', ANTI_SENSE='anti-sense')
```

## PRIME

```python
PRIME = MavisNamespace(FIVE=5, THREE=3)
```

## START_AA

```python
START_AA = 'M'
```

## STOP_AA

```python
STOP_AA = '*'
```

## GIEMSA_STAIN

```python
GIEMSA_STAIN = MavisNamespace(
    GNEG='gneg',
    GPOS33='gpos33',
    GPOS50='gpos50',
    GPOS66='gpos66',
    GPOS75='gpos75',
    GPOS25='gpos25',
    GPOS100='gpos100',
    ACEN='acen',
    GVAR='gvar',
    STALK='stalk',
)
```

## COLUMNS

```python
COLUMNS = MavisNamespace(
    tracking_id='tracking_id',
    library='library',
    cluster_id='cluster_id',
    cluster_size='cluster_size',
    validation_id='validation_id',
    annotation_id='annotation_id',
    product_id='product_id',
    event_type='event_type',
    pairing='pairing',
    inferred_pairing='inferred_pairing',
    gene1='gene1',
    gene1_direction='gene1_direction',
    gene2='gene2',
    gene2_direction='gene2_direction',
    gene1_aliases='gene1_aliases',
    gene2_aliases='gene2_aliases',
    gene_product_type='gene_product_type',
    transcript1='transcript1',
    transcript2='transcript2',
    fusion_splicing_pattern='fusion_splicing_pattern',
    fusion_cdna_coding_start='fusion_cdna_coding_start',
    fusion_cdna_coding_end='fusion_cdna_coding_end',
    fusion_mapped_domains='fusion_mapped_domains',
    fusion_sequence_fasta_id='fusion_sequence_fasta_id',
    fusion_sequence_fasta_file='fusion_sequence_fasta_file',
    fusion_protein_hgvs='fusion_protein_hgvs',
    annotation_figure='annotation_figure',
    annotation_figure_legend='annotation_figure_legend',
    genes_encompassed='genes_encompassed',
    genes_overlapping_break1='genes_overlapping_break1',
    genes_overlapping_break2='genes_overlapping_break2',
    genes_proximal_to_break1='genes_proximal_to_break1',
    genes_proximal_to_break2='genes_proximal_to_break2',
    break1_chromosome='break1_chromosome',
    break1_position_start='break1_position_start',
    break1_position_end='break1_position_end',
    break1_orientation='break1_orientation',
    exon_last_5prime='exon_last_5prime',
    exon_first_3prime='exon_first_3prime',
    break1_strand='break1_strand',
    break1_seq='break1_seq',
    break2_chromosome='break2_chromosome',
    break2_position_start='break2_position_start',
    break2_position_end='break2_position_end',
    break2_orientation='break2_orientation',
    break2_strand='break2_strand',
    break2_seq='break2_seq',
    opposing_strands='opposing_strands',
    stranded='stranded',
    protocol='protocol',
    disease_status='disease_status',
    tools='tools',
    call_method='call_method',
    break1_ewindow='break1_ewindow',
    break1_ewindow_count='break1_ewindow_count',
    break1_ewindow_practical_coverage='break1_ewindow_practical_coverage',
    break1_homologous_seq='break1_homologous_seq',
    break1_split_read_names='break1_split_read_names',
    break1_split_reads='break1_split_reads',
    break1_split_reads_forced='break1_split_reads_forced',
    break2_ewindow='break2_ewindow',
    break2_ewindow_count='break2_ewindow_count',
    break2_ewindow_practical_coverage='break2_ewindow_practical_coverage',
    break2_homologous_seq='break2_homologous_seq',
    break2_split_read_names='break2_split_read_names',
    break2_split_reads='break2_split_reads',
    break2_split_reads_forced='break2_split_reads_forced',
    contig_alignment_query_consumption='contig_alignment_query_consumption',
    contig_alignment_score='contig_alignment_score',
    contig_alignment_query_name='contig_alignment_query_name',
    contig_read_depth='contig_read_depth',
    contig_break1_read_depth='contig_break1_read_depth',
    contig_break2_read_depth='contig_break2_read_depth',
    contig_alignment_rank='contig_alignment_rank',
    contig_build_score='contig_build_score',
    contig_remap_score='contig_remap_score',
    contig_remap_coverage='contig_remap_coverage',
    contig_remapped_read_names='contig_remapped_read_names',
    contig_remapped_reads='contig_remapped_reads',
    contig_seq='contig_seq',
    contig_strand_specific='contig_strand_specific',
    contigs_assembled='contigs_assembled',
    call_sequence_complexity='call_sequence_complexity',
    spanning_reads='spanning_reads',
    spanning_read_names='spanning_read_names',
    flanking_median_fragment_size='flanking_median_fragment_size',
    flanking_pairs='flanking_pairs',
    flanking_pairs_compatible='flanking_pairs_compatible',
    flanking_pairs_read_names='flanking_pairs_read_names',
    flanking_pairs_compatible_read_names='flanking_pairs_compatible_read_names',
    flanking_stdev_fragment_size='flanking_stdev_fragment_size',
    linking_split_read_names='linking_split_read_names',
    linking_split_reads='linking_split_reads',
    raw_break1_half_mapped_reads='raw_break1_half_mapped_reads',
    raw_break1_split_reads='raw_break1_split_reads',
    raw_break2_half_mapped_reads='raw_break2_half_mapped_reads',
    raw_break2_split_reads='raw_break2_split_reads',
    raw_flanking_pairs='raw_flanking_pairs',
    raw_spanning_reads='raw_spanning_reads',
    untemplated_seq='untemplated_seq',
    filter_comment='filter_comment',
    cdna_synon='cdna_synon',
    protein_synon='protein_synon',
    supplementary_call='supplementary_call',
    net_size='net_size',
    repeat_count='repeat_count',
    assumed_untemplated='assumed_untemplated',
)
```

## class MavisNamespace

Namespace to hold module constants

### MavisNamespace.\_\_init\_\_()

```python
def __init__(self, *pos, **kwargs):
```


### MavisNamespace.discard()

Remove a variable if it exists

```python
def discard(self, attr):
```

**Args**

- attr

### MavisNamespace.get\_env\_name()

Get the name of the corresponding environment variable

```python
def get_env_name(self, attr):
```

**Args**

- attr

**Examples**

```python
nspace = MavisNamespace(a=1)
nspace.get_env_name('a')
'MAVIS_A'
```


### MavisNamespace.get\_env\_var()

retrieve the environment variable definition of a given attribute

```python
def get_env_var(self, attr):
```

**Args**

- attr

### MavisNamespace.parse\_listable\_string()

Given some string, parse it into a list

```python
@classmethod
def parse_listable_string(cls, string, cast_type=str, nullable=False):
```

**Args**

- string
- cast_type
- nullable

**Examples**

```python
MavisNamespace.parse_listable_string('1,2,3', int)
[1, 2, 3]
MavisNamespace.parse_listable_string('1;2,None', int, True)
[1, 2, None]
```











### MavisNamespace.copy\_from()

Copy variables from one namespace onto the current namespace

```python
def copy_from(self, source, attrs=None):
```

**Args**

- source
- attrs

### MavisNamespace.get()

get an attribute, return a default (if given) if the attribute does not exist

```python
def get(self, key, *pos):
```

**Args**

- key

**Examples**

```python
nspace = MavisNamespace(thing=1, otherthing=2)
nspace.get('thing', 2)
1
nspace.get('nonexistant_thing', 2)
2
nspace.get('nonexistant_thing')
Traceback (most recent call last):
....
```


### MavisNamespace.keys()

get the attribute keys as a list

```python
def keys(self):
```

**Examples**

```python
MavisNamespace(thing=1, otherthing=2).keys()
['thing', 'otherthing']
```


### MavisNamespace.values()

get the attribute values as a list

```python
def values(self):
```

**Examples**

```python
MavisNamespace(thing=1, otherthing=2).values()
[1, 2]
```


### MavisNamespace.enforce()

checks that the current namespace has a given value

```python
def enforce(self, value):
```

**Args**

- value

**Returns**

: the input value

**Raises**

- `KeyError`: the value did not exist

**Examples**

```python
nspace = MavisNamespace(thing=1, otherthing=2)
nspace.enforce(1)
1
nspace.enforce(3)
Traceback (most recent call last):
....
```


### MavisNamespace.reverse()

for a given value, return the associated key

```python
def reverse(self, value):
```

**Args**

- value: the value to get the key/attribute name for

**Raises**

- `KeyError`: the value is not unique
- `KeyError`: the value is not assigned

**Examples**

```python
nspace = MavisNamespace(thing=1, otherthing=2)
nspace.reverse(1)
'thing'
```




### MavisNamespace.type()

returns the type

```python
def type(self, attr, *pos):
```

**Args**

- attr

**Examples**

```python
nspace = MavisNamespace(thing=1, otherthing=2)
nspace.type('thing')
<class 'int'>
```


### MavisNamespace.define()

Get the definition of a given attribute or return a default (when given) if the attribute does not exist

```python
def define(self, attr, *pos):
```

**Args**

- attr

**Returns**

- `str`: definition for the attribute

**Raises**

- `KeyError`: the attribute does not exist and a default was not given

**Examples**

```python
nspace = MavisNamespace()
nspace.add('thing', 1, defn='I am a thing')
nspace.add('otherthing', 2)
nspace.define('thing')
'I am a thing'
nspace.define('otherthing')
Traceback (most recent call last):
....
nspace.define('otherthing', 'I am some other thing')
'I am some other thing'
```


### MavisNamespace.add()

Add an attribute to the name space

```python
def add(
    self,
    attr,
    value,
    defn=None,
    cast_type=None,
    nullable=False,
    env_overwritable=False,
    listable=False,
):
```

**Args**

- attr (`str`): name of the attribute being added
- value: the value of the attribute
- defn (`str`): the definition, will be used in generating documentation and help menus
- cast_type (`callable`): the function to use in casting the value
- nullable (`bool`): True if this attribute can have a None value
- env_overwritable (`bool`): True if this attribute will be overriden by its environment variable equivalent
- listable (`bool`): True if this attribute can have multiple values

**Examples**

```python
nspace = MavisNamespace()
nspace.add('thing', 1, int, 'I am a thing')
nspace = MavisNamespace()
nspace.add('thing', 1, int)
nspace = MavisNamespace()
nspace.add('thing', 1)
nspace = MavisNamespace()
nspace.add('thing', value=1, cast_type=int, defn='I am a thing')
```




## float\_fraction()

cast input to a float

```python
def float_fraction(num):
```

**Args**

- num: input to cast

**Returns**

: float

**Raises**

- `TypeError`: if the input cannot be cast to a float or the number is not between 0 and 1

## reverse\_complement()

wrapper for the Bio.Seq reverse_complement method

```python
def reverse_complement(s):
```

**Args**

- s (`str`): the input DNA sequence

**Returns**

: :class:`str`: the reverse complement of the input sequence Warning: assumes the input is a DNA sequence

**Examples**

```python
reverse_complement('ATCCGGT')
'ACCGGAT'
```


## translate()

given a DNA sequence, translates it and returns the protein amino acid sequence

```python
def translate(s, reading_frame=0):
```

**Args**

- s (`str`): the input DNA sequence
- reading_frame (`int`): where to start translating the sequence

**Returns**

- `str`: the amino acid sequence
