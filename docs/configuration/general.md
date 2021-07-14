# Getting Started

An exhaustive list of the various configurable settings can be found [here](../settings)

## Pipeline Configuration File

The pipeline can be run in steps or it can be configured using a JSON
configuration file and setup in a single step. Scripts will be generated
to run all steps following clustering.

The config schema is found in the mavis package under `mavis/schemas/config.json`

Top level settings follow the pattern `<section>.<setting>`. The convert and library
sections are nested objects.

## Adjusting the Resource Requirements

### Choosing the Number of Validation/Annotation Jobs

MAVIS chooses the number of jobs to split validate/annotate stages into
based on two settings: [cluster.max_files](../../configuration/settings/#clustermax_files) and
[cluster.min_clusters_per_file](../../configuration/settings/#clustermin-clusters-per-file).

For example, in the following situation say you have: 1000 clusters,
`cluster.max_files=10`, and `cluster.min_clusters_per_file=10`. Then MAVIS will set up
10 validation jobs each with 100 events.

However, if `cluster.min_clusters_per_file=500`, then MAVIS would only set up 2
jobs each with 500 events. This is because
[cluster.min_clusters_per_file](../../configuration/settings/#clustermin-clusters-per-file) takes precedence
over [custer.max_files](../../configuration/settings/#clustermax_files).

Splitting into more jobs will lower the resource requirements per job
(see [resource requirements](../performance/)). The memory and time requirements for validation are linear
with respect to the number of events to be validated.

### Uninformative Filter

For example, if the user is only interested in events in genes, then the
[cluster.uninformative_filter](../../configuration/settings/#clusteruninformative_filter) can be used. This
will drop all events that are not within a certain distance
([cluster.max_proximity](../../configuration/settings/#clustermax_proximity)) to any annotation in
the annotations reference file. These events will be dropped prior to
the validation stage which results in significant speed up.
