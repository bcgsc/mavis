# mavis.annotate.file_io

module which holds all functions relating to loading reference files

## REFERENCE_DEFAULTS

```python
REFERENCE_DEFAULTS = WeakMavisNamespace()
```

## class ReferenceFile







### ReferenceFile.load()

load (or return) the contents of a reference file and add it to the cache if enabled

```python
def load(self, ignore_cache=False, verbose=True):
```

**Args**

- ignore_cache
- verbose


## load\_masking\_regions()

reads a file of regions. The expect input format for the file is tab-delimited and
the header should contain the following columns

- chr: the chromosome
- start: start of the region, 1-based inclusive
- end: end of the region, 1-based inclusive
- name: the name/label of the region

For example:

.. code-block:: text

#chr    start   end     name
chr20   25600000        27500000        centromere

```python
def load_masking_regions(*filepaths):
```

**Returns**

: :class:`dict` of :class:`list` of :class:`BioInterval` by :class:`str`: a dictionary keyed by chromosome name with values of lists of regions on the chromosome

**Examples**

```python
m = load_masking_regions('filename')
m['1']
[BioInterval(), BioInterval(), ...]
```


## load\_reference\_genes()

*Deprecated* Use :func:`load_annotations` instead

```python
def load_reference_genes(*pos, **kwargs):
```

## load\_annotations()

loads gene models from an input file. Expects a tabbed or json file.

```python
def load_annotations(*filepaths, warn=DEVNULL, reference_genome=None, best_transcripts_only=False):
```

**Returns**

: :class:`dict` of :class:`list` of :class:`~mavis.annotate.genomic.Gene` by :class:`str`: lists of genes keyed by chromosome name

## parse\_annotations\_json()

parses a json of annotation information into annotation objects

```python
def parse_annotations_json(data, reference_genome=None, best_transcripts_only=False, warn=DEVNULL):
```

**Args**

- data
- reference_genome
- best_transcripts_only
- warn

## convert\_tab\_to\_json()

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

```python
def convert_tab_to_json(filepath, warn=DEVNULL):
```

**Args**

- filepath (`str`): path to the input tab-delimited file
- warn

**Returns**

: :class:`dict` of :class:`list` of :any:`Gene` by :class:`str`: a dictionary keyed by chromosome name with values of list of genes on the chromosome

**Examples**

```python
ref = load_reference_genes('filename')
ref['1']
[Gene(), Gene(), ....]
Warning:
does not load translations unless then start with 'M', end with '*' and have a length of multiple 3
```



## load\_templates()

primarily useful if template drawings are required and is not necessary otherwise
assumes the input file is 0-indexed with [start,end) style. Columns are expected in
the following order, tab-delimited. A header should not be given

1. name
2. start
3. end
4. band_name
5. giemsa_stain

for example

.. code-block:: text

chr1    0       2300000 p36.33  gneg
chr1    2300000 5400000 p36.32  gpos25

```python
def load_templates(*filepaths):
```

**Returns**

: :class:`list` of :class:`Template`: list of the templates loaded
