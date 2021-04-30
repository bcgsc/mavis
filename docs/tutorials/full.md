# MAVIS (Full) Tutorial

The following tutorial is an introduction to running MAVIS. You will
need to download the tutorial data. Additionally the instructions
pertain to running MAVIS on a [SLURM](../../glossary/#slurm)
cluster. This tutorial will require more resources than the
[mini-tutorial](../../tutorials/mini/) above.

## Getting the Tutorial Data

The tutorial data can be downloaded from the link below. Note that it
may take a while as the download is \~29GB

```text
wget http://www.bcgsc.ca/downloads/mavis/tutorial_data.tar.gz
tar -xvzf tutorial_data.tar.gz
```

The expected contents are

| Path                               | Description                                                                                                              |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| README                             | Information regarding the other files in the directory                                                                   |
| L1522785992\_expected\_events.tab  | The events that we expect to find, either experimentally validated or 'spiked' in                                        |
| L1522785992\_normal.sorted.bam     | Paired normal library BAM file                                                                                           |
| L1522785992\_normal.sorted.bam.bai | BAM index                                                                                                                |
| L1522785992\_trans.sorted.bam      | Tumour transcriptome BAM file                                                                                            |
| L1522785992\_trans.sorted.bam.bai  | BAM index file                                                                                                           |
| L1522785992\_tumour.sorted.bam     | Tumour genome BAM file                                                                                                   |
| L1522785992\_tumour.sorted.bam.bai | BAM index file                                                                                                           |
| breakdancer-1.4.5/                 | Contains the [BreakDancer](../../glossary/#breakdancer) output which was run on the tumour genome BAM file               |
| breakseq-2.2/                      | Contains the [BreakSeq](../../glossary/#breakseq) output which was run on the tumour genome BAM file                     |
| chimerascan-0.4.5/                 | Contains the [ChimeraScan](../../glossary/#chimerascan) output which was run on the tumour transcriptome BAM file        |
| defuse-0.6.2/                      | Contains the [deFuse](../../glossary/#defuse) output which was run on the tumour transcriptome BAM file                  |
| manta-1.0.0/                       | Contains the [Manta](../../glossary/#manta) output which was run on the tumour genome and paired normal genome BAM files |

## Downloading the Reference Inputs

Run the following to download the hg19 reference files

```bash
wget https://raw.githubusercontent.com/bcgsc/mavis/master/tools/get_hg19_reference_files.sh
bash get_hg19_reference_files.sh
```

## Creating the Config File

Most settings can be left as defaults, however you will need to fill out the `libraries` and
`convert` sections to tell MAVIS how to convert your inputs and what libraries to expect.

### Libraries Settings

For this example, because we want to determine which events are
germline/somatic we are going to pass all genome calls to both genomes.
We can use either full file paths (if the input is already in the
standard format) or the alias from a conversion (the first argument
given to the convert option)

```json
{
    "libraries": {
        "L1522785992-normal": { // keyed by library name
            "assign": [ // these are the names of the input files (or conversion aliases) to check for this library
                "breakdancer",
                "breakseq",
                "manta"
            ],
            "bam_file": "tutorial_data/L1522785992_normal.sorted.bam",
            "disease_status": "normal",
            "protocol": "genome"
        },
        "L1522785992-trans": {
            "assign": [
                "chimerascan",
                "defuse"
            ],
            "bam_file": "tutorial_data/L1522785992_trans.sorted.bam",
            "disease_status": "diseased",
            "protocol": "transcriptome",
            "strand_specific": true
        },
        "L1522785992-tumour": {
            "assign": [
                "breakdancer",
                "breakseq",
                "manta"
            ],
            "bam_file": "tutorial_data/L1522785992_tumour.sorted.bam",
            "disease_status": "diseased",
            "protocol": "genome"
        }
    }
}
```

### Convert Settings

If they are raw tool output as in the current example you will need to
use the convert argument to tell MAVIS the file type

```json
{
    "convert": {
        "breakdancer": {  // conversion alias/key
            "assume_no_untemplated": true,
            "file_type": "breakdancer",  // input/file type
            "inputs": [
                "tutorial_data/breakdancer-1.4.5/*txt"
            ]
        },
        "breakseq": {
            "assume_no_untemplated": true,
            "file_type": "breakseq",
            "inputs": [
                "tutorial_data/breakseq-2.2/breakseq.vcf.gz"
            ]
        },
        "chimerascan": {
            "assume_no_untemplated": true,
            "file_type": "chimerascan",
            "inputs": [
                "tutorial_data/chimerascan-0.4.5/chimeras.bedpe"
            ]
        },
        "defuse": {
            "assume_no_untemplated": true,
            "file_type": "defuse",
            "inputs": [
                "tutorial_data/defuse-0.6.2/results.classify.tsv"
            ]
        },
        "manta": {
            "assume_no_untemplated": true,
            "file_type": "manta",
            "inputs": [
                "tutorial_data/manta-1.0.0/diploidSV.vcf.gz",
                "tutorial_data/manta-1.0.0/somaticSV.vcf"
            ]
        }
    }
}
```

### Top-level Settings

Finally you will need to set output directory and the reference files

```json
{
  "output_dir": "output_dir_full",  // where to output files
  "reference.aligner_reference": [
      "reference_inputs/hg19.2bit"
  ],
  "reference.annotations": [
      "reference_inputs/ensembl69_hg19_annotations.json"
  ],
  "reference.dgv_annotation": [
      "reference_inputs/dgv_hg19_variants.tab"
  ],
  "reference.masking": [
      "reference_inputs/hg19_masking.tab"
  ],
  "reference.reference_genome": [
      "reference_inputs/hg19.fa"
  ],
  "reference.template_metadata": [
      "reference_inputs/cytoBand.txt"
  ]
}
```

## Running the Workflow

In order to run the snakemake file you will need to have the config validation module
`mavis_config` installed which has minimal dependencies.

```bash
pip install mavis_config
```

You are now ready to run the workflow

```bash
snakemake --jobs 100 --configfile=tests/full-tutorial.config.json
```

## Analyzing the Output

The best place to start with looking at the MAVIS output is the summary
folder which contains the final results. For column name definitions see
the [glossary](../../outputs/columns).

```text
output_dir/summary/mavis_summary_all_L1522785992-normal_L1522785992-trans_L1522785992-tumour.tab
```
