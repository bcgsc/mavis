# mavis.util

## ENV_VAR_PREFIX

```python
ENV_VAR_PREFIX = 'MAVIS_'
```

## LOG

```python
LOG = Log()
```

## DEVNULL

```python
DEVNULL = Log(level=None)
```

## class Log

wrapper aroung the builtin logging to make it more readable

### Log.\_\_init\_\_()

```python
def __init__(self, indent_str='  ', indent_level=0, level=logging.INFO):
```

**Args**

- indent_str
- indent_level
- level







## class NullableType




## class WeakMavisNamespace

**inherits** `MavisNamespace`




## cast()

cast a value to a given type

```python
def cast(value, cast_func):
```

**Args**

- value
- cast_func

**Examples**

```python
cast('1', int)
1
```


## soft\_cast()

cast a value to a given type, if the cast fails, cast to null

```python
def soft_cast(value, cast_type):
```

**Args**

- value
- cast_type

**Examples**

```python
cast(None, int)
None
cast('', int)
None
```



## bash\_expands()

expand a file glob expression, allowing bash-style brackets.

```python
def bash_expands(*expressions):
```

**Returns**

- `list`: a list of files

**Examples**

```python
bash_expands('./{test,doc}/*py')
[...]
```


## log\_arguments()

output the arguments to the console

```python
def log_arguments(args):
```

**Args**

- args (`Namespace`): the namespace to print arguments for

## mkdirp()

Make a directory or path of directories. Suppresses the error that is normally raised when the directory already exists

```python
def mkdirp(dirname):
```

**Args**

- dirname

## filter\_on\_overlap()

filter a set of breakpoint pairs based on overlap with a set of genomic regions

```python
def filter_on_overlap(bpps, regions_by_reference_name):
```

**Args**

- bpps (`:class:`list` of :class:`~mavis.breakpoint.BreakpointPair``): list of breakpoint pairs to be filtered
- regions_by_reference_name (`:class:`dict` of :class:`list` of :class:`~mavis.annotate.base.BioInterval` by :class:`str``): regions to filter against




## get\_connected\_components()

for a dictionary representing an adjacency matrix of undirected edges returns the connected components

```python
def get_connected_components(adj_matrix):
```

**Args**

- adj_matrix

## generate\_complete\_stamp()

writes a complete stamp, optionally including the run time if start_time is given

```python
def generate_complete_stamp(output_dir, log=DEVNULL, prefix='MAVIS.', start_time=None):
```

**Args**

- output_dir (`str`): path to the output dir the stamp should be written in
- log (`function`): function to print logging messages to
- prefix (`str`): prefix for the stamp name
- start_time (`int`): the start time

**Examples**

```python
generate_complete_stamp('some_output_dir')
'some_output_dir/MAVIS.COMPLETE'
```




## read\_bpp\_from\_input\_file()

reads a file using the tab module. Each row is converted to a breakpoint pair and
other column data is stored in the data attribute

```python
def read_bpp_from_input_file(
    filename, expand_orient=False, expand_strand=False, expand_svtype=False, **kwargs
):
```

**Args**

- filename (`str`): path to the input file
- expand_orient
- expand_strand
- expand_svtype

**Returns**

: :class:`list` of :any:`BreakpointPair`: a list of pairs

**Examples**

```python
read_bpp_from_input_file('filename')
[BreakpointPair(), BreakpointPair(), ...]
One can also validate other expected columns that will go in the data attribute using the usual arguments
to the tab.read_file function
.. code-block:: python
read_bpp_from_input_file('filename', cast={'index': int})
[BreakpointPair(), BreakpointPair(), ...]
```
