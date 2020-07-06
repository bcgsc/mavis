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

| Path                               | Description                                                                                                         |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| README                             | Information regarding the other files in the directory                                                              |
| L1522785992\_expected\_events.tab  | The events that we expect to find, either experimentally validated or 'spiked' in                                 |
| L1522785992\_normal.sorted.bam     | Paired normal library BAM file                                                                                      |
| L1522785992\_normal.sorted.bam.bai | BAM index                                                                                                           |
| L1522785992\_trans.sorted.bam      | Tumour transcriptome BAM file                                                                                       |
| L1522785992\_trans.sorted.bam.bai  | BAM index file                                                                                                      |
| L1522785992\_tumour.sorted.bam     | Tumour genome BAM file                                                                                              |
| L1522785992\_tumour.sorted.bam.bai | BAM index file                                                                                                      |
| breakdancer-1.4.5/                 | Contains the [BreakDancer](../../glossary/#breakdancer) output which was run on the tumour genome BAM file               |
| breakseq-2.2/                      | Contains the [BreakSeq](../../glossary/#breakseq) output which was run on the tumour genome BAM file                     |
| chimerascan-0.4.5/                 | Contains the [ChimeraScan](../../glossary/#chimerascan) output which was run on the tumour transcriptome BAM file        |
| defuse-0.6.2/                      | Contains the [deFuse](../../glossary/#defuse) output which was run on the tumour transcriptome BAM file                  |
| manta-1.0.0/                       | Contains the [Manta](../../glossary/#manta) output which was run on the tumour genome and paired normal genome BAM files |

## Downloading the Reference Inputs

Run the following to download the hg19 reference files and set up the
environment variables for configuring MAVIS

```bash
wget https://raw.githubusercontent.com/bcgsc/mavis/master/tools/get_hg19_reference_files.sh
bash get_hg19_reference_files.sh
source reference_inputs/hg19_env.sh
```

## Generating the Config File

The [config](../../background/citations/#pipeline-config) command
does most of the work of creating the config for you but there are a few
things you need to tell it

1.  **Where your bams are and what library they belong to**

```text
--library L1522785992-normal genome normal False tutorial_data/L1522785992_normal.sorted.bam
--library L1522785992-tumour genome diseased False tutorial_data/L1522785992_tumour.sorted.bam
--library L1522785992-trans transcriptome diseased True tutorial_data/L1522785992_trans.sorted.bam
```

1.  **Where your SV caller output files (events) are**

If they are raw tool output as in the current example you will need to
use the convert argument to tell MAVIS the file type

```text
--convert breakdancer tutorial_data/breakdancer-1.4.5/*txt breakdancer
--convert breakseq tutorial_data/breakseq-2.2/breakseq.vcf.gz breakseq
--convert chimerascan tutorial_data/chimerascan-0.4.5/chimeras.bedpe chimerascan
--convert defuse tutorial_data/defuse-0.6.2/results.classify.tsv defuse
--convert manta tutorial_data/manta-1.0.0/diploidSV.vcf.gz tutorial_data/manta-1.0.0/somaticSV.vcf manta
```

!!! note
    For older versions of MAVIS the convert command may require the path to
    the file(s) be quoted and the strandedness be specified (default is
    False)


3.  **Which events you should validate in which libraries**

For this example, because we want to determine which events are
germline/somatic we are going to pass all genome calls to both genomes.
We can use either full file paths (if the input is already in the
standard format) or the alias from a conversion (the first argument
given to the convert option)

```text
--assign L1522785992-trans chimerascan defuse
--assign L1522785992-tumour breakdancer breakseq manta
--assign L1522785992-normal breakdancer breakseq manta
```

Putting this altogether with a name to call the config, we have the
command to generate the pipeline config. You should expect this step
with these inputs to take about \~5GB memory.

```bash
mavis config \
    --library L1522785992-normal genome normal False tutorial_data/L1522785992_normal.sorted.bam \
    --library L1522785992-tumour genome diseased False tutorial_data/L1522785992_tumour.sorted.bam \
    --library L1522785992-trans transcriptome diseased True tutorial_data/L1522785992_trans.sorted.bam \
    --convert breakdancer tutorial_data/breakdancer-1.4.5/*txt breakdancer \
    --convert breakseq tutorial_data/breakseq-2.2/breakseq.vcf.gz breakseq \
    --convert chimerascan tutorial_data/chimerascan-0.4.5/chimeras.bedpe chimerascan \
    --convert defuse tutorial_data/defuse-0.6.2/results.classify.tsv defuse \
    --convert manta tutorial_data/manta-1.0.0/diploidSV.vcf.gz tutorial_data/manta-1.0.0/somaticSV.vcf manta \
    --assign L1522785992-trans chimerascan defuse \
    --assign L1522785992-tumour breakdancer breakseq manta  \
    --assign L1522785992-normal breakdancer breakseq manta \
    -w mavis.cfg
```

## Setting Up the Pipeline

The next step is running the setup stage. This will perform conversion, clustering, and creating the
submission scripts for the other stages.

```bash
mavis setup mavis.cfg -o output_dir/
```

At this stage you should have something that looks like this. For
simplicity not all files/directories have been shown.

    output_dir/
    |-- build.cfg
    |-- converted_inputs
    |   |-- breakdancer.tab
    |   |-- breakseq.tab
    |   |-- chimerascan.tab
    |   |-- defuse.tab
    |   `-- manta.tab
    |-- L1522785992-normal_normal_genome
    |   |-- annotate
    |   |   |-- batch-aUmErftiY7eEWvENfSeJwc-1/
    |   |   `-- submit.sh
    |   |-- cluster
    |   |   |-- batch-aUmErftiY7eEWvENfSeJwc-1.tab
    |   |   |-- cluster_assignment.tab
    |   |   |-- clusters.bed
    |   |   |-- filtered_pairs.tab
    |   |   `-- MAVIS-batch-aUmErftiY7eEWvENfSeJwc.COMPLETE
    |   `-- validate
    |       |-- batch-aUmErftiY7eEWvENfSeJwc-1/
    |       `-- submit.sh
    |-- pairing
    |   `-- submit.sh
    `-- summary
        `-- submit.sh

## Submitting Jobs to the Cluster

The last step is simple, ssh to your head node of your
[SLURM](../../glossary/#slurm) cluster (or run locally if you
have configured [remote_head_ssh](../../configuration/settings/#remote_head_ssh) and
run the schedule step. This will submit the jobs and create the
dependency chain

```bash
ssh head_node
mavis schedule -o output_dir --submit
```

The schedule step also acts as a built-in checker and can be run to
check for errors or if the pipeline has completed.

```bash
mavis schedule -o output_dir
```

This should give you output something like below (times may vary) after
your run completed correctly.

```text
                      MAVIS: 2.0.0
                      hostname: gphost08.bcgsc.ca
[2018-06-02 19:47:56] arguments
                        command = 'schedule'
                        log = None
                        log_level = 'INFO'
                        output = 'output_dir/'
                        resubmit = False
                        submit = False
[2018-06-02 19:48:01] validate
                        MV_L1522785992-normal_batch-aUmErftiY7eEWvENfSeJwc (1701000) is COMPLETED
                          200 tasks are COMPLETED
                          run time: 609
                        MV_L1522785992-tumour_batch-aUmErftiY7eEWvENfSeJwc (1701001) is COMPLETED
                          200 tasks are COMPLETED
                          run time: 669
                        MV_L1522785992-trans_batch-aUmErftiY7eEWvENfSeJwc (1701002) is COMPLETED
                          23 tasks are COMPLETED
                          run time: 1307
[2018-06-02 19:48:02] annotate
                        MA_L1522785992-normal_batch-aUmErftiY7eEWvENfSeJwc (1701003) is COMPLETED
                          200 tasks are COMPLETED
                          run time: 622
                        MA_L1522785992-tumour_batch-aUmErftiY7eEWvENfSeJwc (1701004) is COMPLETED
                          200 tasks are COMPLETED
                          run time: 573
                        MA_L1522785992-trans_batch-aUmErftiY7eEWvENfSeJwc (1701005) is COMPLETED
                          23 tasks are COMPLETED
                          run time: 537
[2018-06-02 19:48:07] pairing
                        MP_batch-aUmErftiY7eEWvENfSeJwc (1701006) is COMPLETED
                          run time: 466
[2018-06-02 19:48:07] summary
                        MS_batch-aUmErftiY7eEWvENfSeJwc (1701007) is COMPLETED
                          run time: 465
                      parallel run time: 3545
                      rewriting: output_dir/build.cfg
                      run time (hh/mm/ss): 0:00:11
                      run time (s): 11
```

The parallel run time reported corresponds to the sum of the slowest job
for each stage and does not include any queue time etc.

## Analyzing the Output

The best place to start with looking at the MAVIS output is the summary
folder which contains the final results. For column name definitions see
the [glossary](../../outputs/columns).

    output_dir/summary/mavis_summary_all_L1522785992-normal_L1522785992-trans_L1522785992-tumour.tab
