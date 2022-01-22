# Annotation Only

Sometimes you have a set of variants and would simply like to run the annotate step of MAVIS to visualize and annotate them.

First you need to create your basic config to tell MAVIS where the reference files you want to use are and some minimal information about the library/sample you want to process.

Here is an example config where the user has created a minimal input file in the MAVIS standard input file format. We convert it to expand any unknowns (ex. SV type if left blank)

```json
{
    "libraries": {
        "my_library": {
            "assign": ["my_converted_file"],
            "disease_status": "normal",
            "protocol": "genome"
        }
    },
    "convert": {
        "my_converted_file": {
            "inputs": ["/path/to/file/structural_variants.txt"],
            "file_type": "mavis"
         }
    },
    "cluster.split_only": true,
    "skip_stage.validate": true,
    "output_dir": "my_output_dir",
    "reference.annotations": "/path/to/mavis/reference_files/ensembl79_hg38_annotations.json",
    "reference.template_metadata": "/path/to/mavis/reference_files/hg38_cytoBand.txt",
    "reference.reference_genome": "/path/to/hg38_no_alt/genome/hg38_no_alt.fa",
    "reference.masking": "/path/to/mavis/reference_files/masking_hg38.adjusted.tab",
    "reference.dgv_annotation": "/path/to/mavis/reference_files/dgv_hg38_annotations.tab"
}
```

Another example is given in the MAVIS tests folder under `tests/mini-tutorial.annotate_only.config.json` which looks like this

```json
{
    "annotate.draw_fusions_only": false,
    "convert": {
        "mock_converted": {
            "inputs": [
                "tests/data/mock_sv_events.tsv"
            ],
            "file_type": "mavis",
            "assume_no_untemplated": true
        }
    },
    "skip_stage.validate": true,
    "cluster.uninformative_filter": true,
    "cluster.limit_to_chr": null,
    "cluster.min_clusters_per_file": 5,
    "libraries": {
        "mock-A47933": {
            "assign": [
                "tests/data/mock_trans_sv_events.tsv"
            ],
            "bam_file": "tests/data/mock_trans_reads_for_events.sorted.bam",
            "disease_status": "diseased",
            "protocol": "transcriptome",
            "strand_specific": true
        },
        "mock-A36971": {
            "assign": [
                "mock_converted"
            ],
            "bam_file": "tests/data/mock_reads_for_events.sorted.bam",
            "disease_status": "diseased",
            "protocol": "genome",
            "strand_specific": false
        }
    },
    "output_dir": "output_dir",
    "reference.annotations": [
        "tests/data/mock_annotations.json"
    ],
    "reference.dgv_annotation": [
        "tests/data/mock_dgv_annotation.txt"
    ],
    "reference.masking": [
        "tests/data/mock_masking.tab"
    ],
    "reference.reference_genome": [
        "tests/data/mock_reference_genome.fa"
    ],
    "reference.template_metadata": [
        "tests/data/cytoBand.txt"
    ]
}
```

Either of these configurations can be run with the following command simply by changing the configfile argument

```bash
snakemake -j 1 \
    --configfile tests/mini-tutorial.annotate_only.config.json \
    -s Snakefile
```
