# Migrating

## Migrating from v2 to v3

There are major changes from v2 to v3 of MAVIS.

### Tab File Headers

Tab file headers no longer start with `#`. Any lines starting with a pound will be treated
as comments. This will apply to mavis-style inputs as well as any tab delimited
reference files

### Configuration

MAVIS no longer uses command line arguments, config files, and environment variables for
configuration. Instead all configurable settings are controlled via a single input JSON
config file

### Scheduling

MAVIS is now integrated with snakemake instead of handling its own scheduling

## Reference Annotation Files

MAVIS no longer supports the previously deprecated tab-delimited format of the annotations file. If you are still using these files in your project we have provided a script to automatically convert them to the newer format in the tools directory

```bash
python src/tools/migrate_mavis_annotations_to_jsonl.py \
    /path/to/tab/file.tab \
    /path/to/new/json/file.json
```
