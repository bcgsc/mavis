# Non-python Dependencies

MAVIS integrates with
[SV callers](./sv_callers.md),
[job schedulers](#job-schedulers), and
[aligners](#aligners). While some of
these dependencies are optional, all currently supported options are
detailed below. The versions column in the tables below list all the
versions which were tested for each tool. Each version listed is known
to be compatible with MAVIS.

## Job Schedulers

MAVIS v3 uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to handle job scheduling

## Aligners

Two aligners are supported [bwa](../../glossary/#bwa) and
[blat](../../glossary/#blat) (default). These are both included in the docker image by default.

| Name                                           | Version(s)              | Environment Setting       |
| ---------------------------------------------- | ----------------------- | ------------------------- |
| [blat](../../glossary/#blat)                   | `36x2` `36`             | `MAVIS_ALIGNER=blat`      |
| [bwa mem <bwa>](../../glossary/#bwa mem <bwa>) | `0.7.15-r1140` `0.7.12` | `MAVIS_ALIGNER='bwa mem'` |

!!! note
    When setting the aligner you will also need to set the
    [aligner_reference](../../configuration/settings/#aligner_reference) to match
