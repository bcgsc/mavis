# mavis.assemble

## class Contig







### Contig.remap\_depth()

the average depth of remapped reads over a give range of the contig sequence

```python
def remap_depth(self, query_range=None):
```

**Args**

- query_range (`Interval`): 1-based inclusive range


## class DeBruijnGraph

**inherits** `nx.DiGraph`

wrapper for a basic digraph
enforces edge weights

### DeBruijnGraph.get\_edge\_freq()

returns the freq from the data attribute for a specified edge

```python
def get_edge_freq(self, n1, n2):
```

**Args**

- n1
- n2

### DeBruijnGraph.add\_edge()

add a given edge to the graph, if it exists add the frequency to the existing frequency count

```python
def add_edge(self, n1, n2, freq=1):
```

**Args**

- n1
- n2
- freq


### DeBruijnGraph.trim\_tails\_by\_freq()

for any paths where all edges are lower than the minimum weight trim

```python
def trim_tails_by_freq(self, min_weight):
```

**Args**

- min_weight (`int`): the minimum weight for an edge to be retained

### DeBruijnGraph.trim\_forks\_by\_freq()

for all nodes in the graph, if the node has an out-degree > 1 and one of the outgoing
edges has freq < min_weight. then that outgoing edge is deleted

```python
def trim_forks_by_freq(self, min_weight):
```

**Args**

- min_weight

### DeBruijnGraph.trim\_noncutting\_paths\_by\_freq()

trim any low weight edges where another path exists between the source and target
of higher weight

```python
def trim_noncutting_paths_by_freq(self, min_weight):
```

**Args**

- min_weight

### DeBruijnGraph.get\_sinks()

returns all nodes with an outgoing degree of zero

```python
def get_sinks(self, subgraph=None):
```

**Args**

- subgraph

### DeBruijnGraph.get\_sources()

returns all nodes with an incoming degree of zero

```python
def get_sources(self, subgraph=None):
```

**Args**

- subgraph


## digraph\_connected\_components()

the networkx module does not support deriving connected
components from digraphs (only simple graphs)
this function assumes that connection != reachable
this means there is no difference between connected components
in a simple graph and a digraph

```python
def digraph_connected_components(graph, subgraph=None):
```

**Args**

- graph (`networkx.DiGraph`): the input graph to gather components from
- subgraph

**Returns**

: :class:`list` of :class:`list`: returns a list of compnents which are lists of node names

## pull\_contigs\_from\_component()

builds contigs from the a connected component of the assembly DeBruijn graph

```python
def pull_contigs_from_component(
    assembly, component, min_edge_trim_weight, assembly_max_paths, log=DEVNULL
):
```

**Args**

- assembly (`DeBruijnGraph`): the assembly graph
- component (`list`): list of nodes which make up the connected component
- min_edge_trim_weight (`int`): the minimum weight to not remove a non cutting edge/path
- assembly_max_paths (`int`): the maximum number of paths allowed before the graph is further simplified
- log (`function`): the log function

**Returns**

: :class:`Dict` of :class:`int` by :class:`str`: the paths/contigs and their scores

## filter\_contigs()

given a list of contigs, removes similar contigs to leave the highest (of the similar) scoring contig only

```python
def filter_contigs(contigs, assembly_min_uniq=0.01):
```

**Args**

- contigs
- assembly_min_uniq

## assemble()

for a set of sequences creates a DeBruijnGraph
simplifies trailing and leading paths where edges fall
below a weight threshold and the return all possible unitigs/contigs

drops any sequences too small to fit the kmer size

```python
def assemble(
    sequences,
    kmer_size,
    min_edge_trim_weight=3,
    assembly_max_paths=20,
    assembly_min_uniq=0.01,
    min_complexity=0,
    log=lambda *pos, **kwargs: None,
    **kwargs
):
```

**Args**

- sequences (`:class:`list` of :class:`str``): a list of strings/sequences to assemble
- kmer_size: see :term:`assembly_kmer_size` the size of the kmer to use
- min_edge_trim_weight: see :term:`assembly_min_edge_trim_weight`
- assembly_max_paths: see :term:`assembly_max_paths`
- assembly_min_uniq
- min_complexity
- log (`function`): the log function

**Returns**

: :class:`list` of :class:`Contig`: a list of putative contigs

## kmers()

for a sequence, compute and return a list of all kmers of a specified size

```python
def kmers(s, size):
```

**Args**

- s (`str`): the input sequence
- size (`int`): the size of the kmers

**Returns**

: :class:`list` of :class:`str`: the list of kmers

**Examples**

```python
kmers('abcdef', 2)
['ab', 'bc', 'cd', 'de', 'ef']
```
