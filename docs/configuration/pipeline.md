# Running the Pipeline

## Running MAVIS using a Job Scheduler

MAVIS v3 uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to handle job scheduling
and setup

The MAVIS pipeline is highly configurable. Some pipeline steps
(cluster, validate) are optional and can be automatically skipped. The
standard pipeline is
far-left.

The most common use case is running the pipeline through snakemake

```bash
snakemake -j <MAX JOBS> --configfile <YOUR CONFIG>
```

If you are submitting to a cluster, use the [snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)

```bash
snakemake -j <MAX JOBS> --configfile <YOUR CONFIG> --profile <YOUR PROFILE NAME>
```

This will submit a series of jobs with dependencies.

To use the mavis docker container through singularity, instead of installing mavis via pip, add the
[`--use-singularity`](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#running-jobs-in-containers)
flag.

```bash
snakemake -j <MAX JOBS> --configfile <YOUR CONFIG> --profile <YOUR PROFILE NAME> --use-singularity`
```
