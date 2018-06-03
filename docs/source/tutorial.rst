
....

.. _mini-tutorial:

MAVIS (Mini) Tutorial
-----------------------

This tutorial is based on the data included in the tests folder of MAVIS. The data files are very small and
this tutorial is really only intended for testing a MAVIS install. The data here is simulated and results are
not representitive of the typical events you would see reported from MAVIS. For a more complete tutorial with
actual fusion gene examples, please see the :ref:`full-tutorial` below.


The first step is to clone or download a zip of the MAVIS repository (https://github.com/bcgsc/mavis). You will need the tests directory.
The tag you check out should correspond to the MAVIS version you have installed

.. code:: bash

    git clone https://github.com/bcgsc/mavis.git
    git checkout v2.0.0
    mv mavis/tests .
    rm -r mavis

Now you should have a folder called ``tests`` in your current directory. You will need to specify the scheduler
if you want to test one that is not the default. For example

.. code:: bash

    export MAVIS_SCHEDULER=LOCAL

Since this is a trivial example, it can easily be run locally. By default MAVIS in local mode will run a maximum of 1 less than
the current cpu count processes. If you are running other things on the same machine you may find it useful to set this directly.

.. code:: bash

    export MAVIS_CONCURRENCY_LIMIT=2

The above will limit mavis to running 2 processes concurrently.

Now you are ready to run MAVIS itself. This can be done in two commands (since the config file we are going to use is already built).
First set up the pipeline

.. code:: bash

    mavis setup tests/data/pipeline_config.cfg -o output_dir

Now if you run the schedule step (without the submit flag, schedule acts as a checker) you should see something like


.. code:: bash

    mavis schedule -o output_dir/


::

                          MAVIS: 1.8.4
                          hostname: gphost08.bcgsc.ca
    [2018-06-01 12:19:31] arguments
                            command = 'schedule'
                            log = None
                            log_level = 'INFO'
                            output = 'output_dir/'
                            resubmit = False
                            submit = False
    [2018-06-01 12:19:31] validate
                            MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                            MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 is NOT SUBMITTED
    [2018-06-01 12:19:31] annotate
                            MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                            MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 is NOT SUBMITTED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 is NOT SUBMITTED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 is NOT SUBMITTED
    [2018-06-01 12:19:31] pairing
                            MP_batch-s4W2Go4tinn49nkhSuusrE is NOT SUBMITTED
    [2018-06-01 12:19:31] summary
                            MS_batch-s4W2Go4tinn49nkhSuusrE is NOT SUBMITTED
                          rewriting: output_dir/build.cfg

Adding the submit argument will start the pipeline

.. code:: bash

    mavis schedule -o output_dir/ --submit

After this completes, run schedule without the submit flag again and you should see something like

::

                          MAVIS: 1.8.4
                          hostname: gphost08.bcgsc.ca
    [2018-06-01 13:15:28] arguments
                            command = 'schedule'
                            log = None
                            log_level = 'INFO'
                            output = 'output_dir/'
                            resubmit = False
                            submit = False
    [2018-06-01 13:15:28] validate
                            MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 (zQJYndSMimaoALwcSSiYwi) is COMPLETED
                            MV_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 (BHFVf3BmXVrDUA5X4GGSki) is COMPLETED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 (tUpx3iabCrpR9iKu9rJtES) is COMPLETED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 (hgmH7nqPXZ49a8yTsxSUWZ) is COMPLETED
                            MV_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 (cEoRN582An3eAGALaSKmpJ) is COMPLETED
    [2018-06-01 13:15:28] annotate
                            MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-1 (tMHiVR8ueNokhBDnghXYo6) is COMPLETED
                            MA_mock-A36971_batch-s4W2Go4tinn49nkhSuusrE-2 (AsNpNdvUyhNtKmRZqRSPpR) is COMPLETED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-1 (k7qQiAzxfC2dnZwsGH7BzD) is COMPLETED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-2 (dqAuhhcVKejDvHGBXn22xb) is COMPLETED
                            MA_mock-A47933_batch-s4W2Go4tinn49nkhSuusrE-3 (eB69Ghed2xAdp2VRdaCJBf) is COMPLETED
    [2018-06-01 13:15:28] pairing
                            MP_batch-s4W2Go4tinn49nkhSuusrE (6LfEgBtBsmGhQpLQp9rXmi) is COMPLETED
    [2018-06-01 13:15:28] summary
                            MS_batch-s4W2Go4tinn49nkhSuusrE (HDJhXgKjRmseahcQ7mgNoD) is COMPLETED
                          rewriting: output_dir/build.cfg
                          run time (hh/mm/ss): 0:00:00
                          run time (s): 0

If you see the above, then MAVIS has completed correctly!


....


.. _full-tutorial:

MAVIS (Full) Tutorial
-----------------------

The following tutorial is an introduction to running MAVIS. You will need to download the tutorial data. Additionally
the instructions pertain to running MAVIS on a :term:`SLURM` cluster. This tutorial will require more resources than
the :ref:`mini-tutorial` above.

Getting the Tutorial Data
.............................

The tutorial data can be downloaded from the link below. Note that it may take a while as the download is ~29GB

.. code::

    wget http://www.bcgsc.ca/downloads/mavis/tutorial_data.tar.gz
    tar -xvzf tutorial_data.tar.gz


The expected contents are

.. list-table::
    :header-rows: 1

    *   - Path
        - Description
    *   - README
        - Information regarding the other files in the directory
    *   - L1522785992_expected_events.tab
        - The events that we expect to find, either experimentally validated or 'spiked' in
    *   - L1522785992_normal.sorted.bam
        - Paired normal library BAM file
    *   - L1522785992_normal.sorted.bam.bai
        - BAM index
    *   - L1522785992_trans.sorted.bam
        - Tumour transcriptome BAM file
    *   - L1522785992_trans.sorted.bam.bai
        - BAM index file
    *   - L1522785992_tumour.sorted.bam
        - Tumour genome BAM file
    *   - L1522785992_tumour.sorted.bam.bai
        - BAM index file
    *   - breakdancer-1.4.5/
        - Contains the :term:`BreakDancer` output which was run on the tumour genome BAM file
    *   - breakseq-2.2/
        - Contains the :term:`BreakSeq` output which was run on the tumour genome BAM file
    *   - chimerascan-0.4.5/
        - Contains the :term:`ChimeraScan` output which was run on the tumour transcriptome BAM file
    *   - defuse-0.6.2/
        - Contains the :term:`deFuse` output which was run on the tumour transcriptome BAM file
    *   - manta-1.0.0/
        - Contains the :term:`Manta` output which was run on the tumour genome and paired normal genome BAM files


Downloading the Reference Inputs
.................................

Run the following to download the hg19 reference files and set up the environment variables for configuring MAVIS

.. code:: bash

    wget https://raw.githubusercontent.com/bcgsc/mavis/master/tools/get_hg19_reference_files.sh
    bash get_hg19_reference_files.sh
    source reference_inputs/hg19_env.sh


.. _example-generating-the-conf:

Generating the Config File
.............................

The :ref:`config <pipeline-config>` command does most of the work of creating the config for you but there are a few things you need to tell it

1. **Where your bams are and what library they belong to**

.. code:: text

    --library L1522785992-normal genome normal False tutorial_data/L1522785992_normal.sorted.bam
    --library L1522785992-tumour genome diseased False tutorial_data/L1522785992_tumour.sorted.bam
    --library L1522785992-trans transcriptome diseased True tutorial_data/L1522785992_trans.sorted.bam

2. **Where your SV caller output files (events) are**

If they are raw tool output as in the current example you will need to use the convert argument to tell MAVIS the file type

.. code:: text

    --convert breakdancer tutorial_data/breakdancer-1.4.5/*txt breakdancer
    --convert breakseq tutorial_data/breakseq-2.2/breakseq.vcf.gz breakseq
    --convert chimerascan tutorial_data/chimerascan-0.4.5/chimeras.bedpe chimerascan
    --convert defuse tutorial_data/defuse-0.6.2/results.classify.tsv defuse
    --convert manta tutorial_data/manta-1.0.0/diploidSV.vcf.gz tutorial_data/manta-1.0.0/somaticSV.vcf manta

.. note::

    For older versions of MAVIS the convert command may require the path to the file(s) be quoted and the strandedness be specified (default is False)


3. **Which events you should validate in which libraries**

For this example, because we want to determine which events are germline/somatic we are going to pass all genome
calls to both genomes. We can use either full file paths (if the input is already in the standard format)
or the alias from a conversion (the first argument given to the convert option)

.. code:: text

    --assign L1522785992-trans chimerascan defuse
    --assign L1522785992-tumour breakdancer breakseq manta
    --assign L1522785992-normal breakdancer breakseq manta

Putting this altogether with a name to call the config, we have the command to generate the pipeline config. You should
expect this step with these inputs to take about ~5GB memory.

.. code:: bash

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


Setting Up the Pipeline
.........................

The next step is :ref:`running the setup stage <pipeline-standard>`. This will perform conversion, clustering, and creating the
submission scripts for the other stages.

.. code:: bash

    mavis setup mavis.cfg -o output_dir/

At this stage you should have something that looks like this.
For simplicity not all files/directories have been shown.

::

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


Submitting Jobs to the Cluster
..................................

The last step is simple, ssh to your head node of your :term:`SLURM` cluster (or run locally if you have configured
:term:`remote_head_ssh`) and run the schedule step. This will submit the jobs and create the dependency chain

.. code:: bash

    ssh head_node
    mavis schedule -o output_dir --submit

The schedule step also acts as a built-in checker and can be run to check for errors or if the pipeline has completed.

.. code:: bash

    mavis schedule -o output_dir

This should give you output something like below (times may vary) after your run completed correctly.

::

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

The parallel run time reported corresponds to the sum of the slowest job for each stage and does not include any queue time etc.


Analyzing the Output
.....................

The best place to start with looking at the MAVIS output is the summary folder which contains the
final results. For column name definitions see the :ref:`glossary <glossary-column-names>`.

::

    output_dir/summary/mavis_summary_all_L1522785992-normal_L1522785992-trans_L1522785992-tumour.tab
