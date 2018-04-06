
....


MAVIS Tutorial
----------------

The following tutorial is an introduction to running MAVIS. You will need to download the tutorial data. Additionally
the instructions pertain to running MAVIS on a :term:`SLURM` cluster.

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


It is also assumed that the user has set the MAVIS environment variables for the required hg19 reference files. For example, the file
your source might look like

.. code:: bash

    export MAVIS_REFERENCE_GENOME=/reference_inputs/hg19.fa
    export MAVIS_ANNOTATIONS=/reference_inputs/ensembl69_hg19_annotations.json
    export MAVIS_MASKING=/reference_inputs/hg19_masking.tab
    export MAVIS_DGV_ANNOTATION=/reference_inputs/dgv_hg19_variants.tab
    export MAVIS_ALIGNER_REFERENCE=/reference_inputs/hg19.2bit
    export MAVIS_TEMPLATE_METADATA=/reference_inputs/cytoBand.txt


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

.. note::

    If you haven't set defaults for the reference input files (using environment variables) you may need to pass those arguments as well.


Setting Up the Pipeline
.........................

The next step is :ref:`running the pipeline stage <pipeline-standard>`. This will perform conversion, clustering, and creating the :term:`SLURM`/:term:`SGE`
submission scripts for the other stages.

.. code:: bash

    mavis pipeline mavis.config -o output_dir/

At this stage you should have something that looks like this. 
For simplicity not all files/directories have been shown.

.. code:: text

    output_dir/
    |-- converted_inputs
    |   |-- breakdancer.tab
    |   |-- breakseq.tab
    |   |-- chimerascan.tab
    |   |-- defuse.tab
    |   `-- manta.tab
    |-- L1522785992-normal_normal_genome
    |   |-- annotate
    |   |   `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |   |       `-- submit.sh
    |   |-- cluster
    |   |   |-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1.tab
    |   |   |-- cluster_assignment.tab
    |   |   |-- clusters.bed
    |   |   |-- filtered_pairs.tab
    |   |   `-- MAVIS.COMPLETE
    |   `-- validate
    |       `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |           `-- submit.sh
    |-- L1522785992-trans_diseased_transcriptome
    |   |-- annotate
    |   |   `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |   |       `-- submit.sh
    |   |-- cluster
    |   |   |-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1.tab
    |   |   |-- cluster_assignment.tab
    |   |   |-- clusters.bed
    |   |   |-- filtered_pairs.tab
    |   |   `-- MAVIS.COMPLETE
    |   `-- validate
    |       `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |           `-- submit.sh
    |-- L1522785992-tumour_diseased_genome
    |   |-- annotate
    |   |   `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |   |       `-- submit.sh
    |   |-- cluster
    |   |   |-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1.tab
    |   |   |-- cluster_assignment.tab
    |   |   |-- clusters.bed
    |   |   |-- filtered_pairs.tab
    |   |   `-- MAVIS.COMPLETE
    |   `-- validate
    |       `-- batch-EK5Nx7xmfrbX9Vhuz2S7LR-1
    |           `-- submit.sh
    |-- pairing
    |   `-- submit.sh
    |-- submit_pipeline_batch-EK5Nx7xmfrbX9Vhuz2S7LR.sh
    `-- summary
        `-- submit.sh


Submitting Jobs to the Cluster
..................................

The last step is simple, ssh to your head node of your :term:`SLURM` cluster and run
the submit pipeline script. This contains the sbatch commands to submit all the individual 
jobs and chain their dependencies

.. code:: bash

    ssh head_node
    cd output_dir/
    bash submit_pipeline_batch-d8zWDrpkbBuxj7eTahuH3K.sh

You should see the output (with a different job number)

.. code:: bash

    Submitted batch job 1120999

This is the job number of the mavis summary job. When this job is complete your MAVIS run is complete.
To check that everything ran correctly MAVIS has a built-in checker.

.. code:: bash

    mavis checker -o output_dir

This should give you output something like below (times may vary) if your run completed correctly.

Analyzing the Output
.....................

The best place to start with looking at the MAVIS output is the summary folder which contains the 
final results. For column name definitions see the :ref:`glossary <glossary-column-names>`.

.. code:: text

    output_dir/summary/mavis_summary_all_L1522785992-normal_L1522785992-trans_L1522785992-tumour.tab
