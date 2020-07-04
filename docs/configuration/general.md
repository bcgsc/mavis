# Getting Started

An exhaustive list of the various configurable settings can be found [here](../settings)

## Pipeline Configuration File

The pipeline can be run in steps or it can be configured using a
configuration file and setup in a single step. Scripts will be generated
to run all steps following clustering. The configuration file can be
built from scratch or a template can be output as shown below

```bash
mavis config --write template.cfg
```

This will create a template config file called template.cfg which can
then be edited by the user. However this will be a simple config with no
library information. To generate a configuration file with the library
information as well as estimates for the fragment size parameters more
inputs are required (see
[generating the config file](../../tutorials/full/#generating-the-config-file) for more information).

## Environment Variables

Most of the default settings can be changed by using environment
variables. The value given by the environment variables will be used as
the new default. Config or command-line parameters will still override
these settings.

All environment variables are prefixed with MAVIS and an underscore.
Otherwise the variable name is the same as that used for the command
line parameter or config setting (uppercased). For example to change the
default minimum mapping quality used during the validate stage

```bash
export MAVIS_MIN_MAPPING_QUALITY=10
```

## Adjusting the Resource Requirements

### Choosing the Number of Validation/Annotation Jobs

MAVIS chooses the number of jobs to split validate/annotate stages into
based on two settings: [max_files](../../configuration/settings/#max_files) and
[min_clusters_per_file](../../configuration/settings/#min-clusters-per-file).

For example, in the following situation say you have: 1000 clusters,
`max_files=10`, and `min_clusters_per_file=10`. Then MAVIS will set up
10 validation jobs each with 100 events.

However, if `min_clusters_per_file=500`, then MAVIS would only set up 2
jobs each with 500 events. This is because
[min_clusters_per_file](../../configuration/settings/#min-clusters-per-file) takes precedence
over [max_files](../../configuration/settings/#max_files).

Splitting into more jobs will lower the resource requirements per job
(see [resource requirements](../performance/)). The memory and time requirements for validation are linear
with respect to the number of events to be validated.

### Uninformative Filter

For example, if the user is only interested in events in genes, then the
[uninformative_filter](../../configuration/settings/#uninformative_filter) can be used. This
will drop all events that are not within a certain distance
([max_proximity](../../configuration/settings/#max_proximity)) to any annotation in
the annotations reference file. These events will be dropped prior to
the validation stage which results in significant speed up.

This can be set using the environment variable

```bash
export MAVIS_UNINFORMATIVE_FILTER=True
```

or in the pipeline config file

```text
[cluster]
uninformative_filter = True
```

or as a command line argument to the cluster stage

```bash
mavis cluster --uninformative_filter True ....
```
