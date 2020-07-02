# mavis.main

## check\_overlay\_args()

parse the overlay options and check the formatting

```python
def check_overlay_args(args, parser):
```

**Args**

- args
- parser

## overlay\_main()

generates an overlay diagram

```python
def overlay_main(
    gene_name,
    output,
    buffer_length,
    read_depth_plots,
    markers,
    annotations,
    drawing_width_iter_increase,
    max_drawing_retries,
    min_mapping_quality,
    ymax_color='#FF0000',
    **kwargs
):
```

**Args**

- gene_name
- output
- buffer_length
- read_depth_plots
- markers
- annotations
- drawing_width_iter_increase
- max_drawing_retries
- min_mapping_quality
- ymax_color


## main()

sets up the parser and checks the validity of command line args
loads reference files and redirects into subcommand main functions

```python
def main(argv=None):
```

**Args**

- argv (`list`): List of arguments, defaults to command line arguments
