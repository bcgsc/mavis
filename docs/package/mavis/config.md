# mavis.config

## CONVERT_OPTIONS

```python
CONVERT_OPTIONS = WeakMavisNamespace()
```

## class CustomHelpFormatter

**inherits** `argparse.ArgumentDefaultsHelpFormatter`

subclass the default help formatter to stop default printing for required arguments





## class RangeAppendAction

**inherits** `argparse.Action`

allows an argument to accept a range of arguments

### RangeAppendAction.\_\_init\_\_()

```python
def __init__(self, nmin=1, nmax=None, **kwargs):
```

**Args**

- nmin
- nmax



## class LibraryConfig

**inherits** `MavisNamespace`

holds library specific configuration information

### LibraryConfig.\_\_init\_\_()

```python
def __init__(
    self,
    library,
    protocol,
    disease_status,
    bam_file=None,
    inputs=None,
    read_length=None,
    median_fragment_size=None,
    stdev_fragment_size=None,
    strand_specific=False,
    strand_determining_read=2,
    **kwargs
):
```

**Args**

- library
- protocol
- disease_status
- bam_file
- inputs
- read_length
- median_fragment_size
- stdev_fragment_size
- strand_specific
- strand_determining_read



### LibraryConfig.build()

Builds a library config section and gathers the bam stats

```python
@staticmethod
def build(
    library,
    protocol,
    bam_file,
    inputs,
    annotations=None,
    log=DEVNULL,
    distribution_fraction=0.98,
    sample_cap=3000,
    sample_bin_size=1000,
    sample_size=500,
    **kwargs
):
```

**Args**

- library
- protocol
- bam_file
- inputs
- annotations
- log
- distribution_fraction
- sample_cap
- sample_bin_size
- sample_size



## class MavisConfig

**inherits** `MavisNamespace`



### MavisConfig.read()

reads the configuration settings from the configuration file

```python
@staticmethod
def read(filepath):
```

**Args**

- filepath (`str`): path to the input configuration file

**Returns**

- `class`: `list` of :class:`Namespace`: namespace arguments for each library



## validate\_section()

given a dictionary of values, returns a new dict with the values casted to their appropriate type or set
to a default if the value was not given

```python
def validate_section(section, namespace, use_defaults=False):
```

**Args**

- section
- namespace
- use_defaults

## get\_metavar()

For a given argument type, returns the string to be used for the metavar argument in add_argument

```python
def get_metavar(arg_type):
```

**Args**

- arg_type

**Examples**

```python
get_metavar(bool)
'{True,False}'
```


## nameable\_string()

A string that can be used for library and/or filenames

```python
def nameable_string(input_string):
```

**Args**

- input_string

## augment\_parser()

Adds options to the argument parser. Separate function to facilitate the pipeline steps
all having a similar look/feel

```python
def augment_parser(arguments, parser, required=None):
```

**Args**

- arguments
- parser
- required
