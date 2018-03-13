Examples
=========

Tumour Pair Full Pipeline
----------------------------

The following example is a tumour normal genome pair with a tumour transcriptome, all of which was set up to run on :term:`SLURM`.

Before you Start
....................

Before you can run MAVIS, you will need the output from SV callers and :term:`BAM` files for any libraries that you wish to validate. In this example
we used output from :term:`DELLY`, :term:`Manta`, :term:`Trans-ABySS`, :term:`Chimerascan`, and :term:`DeFuse`. The results 
looked something like below (note that for simplicity only relevant files are shown).

For the purposes of this example we have 3 libraries

- NG1: A normal genome (blood)
- TG1: A diseased genome (tumour biopsy)
- TT1: A diseased transcriptome (tumour transcriptome)

.. code:: text

    trans-abyss_TT1/
    |-- fusions/
    |   |-- fusions.tsv
    |   |-- ITD.tsv
    |   |-- local.tsv
    |   |-- LSR.tsv
    |   |-- PTD.tsv
    |   `-- sense_fusion.tsv
    `-- indels/
        `-- events_exons_novel.tsv
    
    trans-abyss_TG1/
    |-- fusions/
    |   |-- fusions.tsv
    |   |-- ITD.tsv
    |   |-- local.tsv
    |   |-- LSR.tsv
    |   |-- PTD.tsv
    |   `-- sense_fusion.tsv
    `-- indels/
        `-- events_exons_novel.tsv

    defuse_TT1/
    `-- results/
        `-- results.classify.tsv

    chimeracan_TT1/
    `-- chimeras.bedpe
    
    manta_TG1_NG1/
    `-- results/
        `-- variants/
            |-- diploidSV.vcf
            `-- somaticSV.vcf

    delly_TG1_NG1/
    |-- *.bcf
    `-- combined.vcf

The appropriate :term:`BAM` files might look like

.. code:: text

    TT1.bam
    TT1.bam.bai
    TG1.bam
    TG1.bam.bai
    NG1.bam
    NG1.bam.bai

If you have all of these inputs then you are ready to generate the MAVIS pipeline config file.

Generating the Config File
.............................

The :ref:`config <pipeline-config>` command does most of the work of creating the config for you but there are a few things you need to tell it

1. **Where your bams are and what library they belong to**

.. code:: text

    --library NG1 genome normal False NG1.bam
    --library TG1 genome diseased False TG1.bam
    --library TT1 transcriptome diseased True TT1.bam

2. **Where your SV caller output files (events) are**
   
If they are raw tool output as in the current example you will need to use the convert argument to tell MAVIS the file type

.. code:: text

    --convert ta_trans trans-abyss_TT1/fusions/*.tsv transabyss
    --convert defuse defuse_TT1/results/results.classify.tsv defuse
    --convert chimera chimeracan_TT1/chimeras.bedpe chimerascan
    --convert delly delly_TG1_NG1/combined.vcf delly
    --convert manta manta_TG1_NG1/results/variant/{Diploid,Somatic}.vcf
    --convert ta_genome trans-abyss_TG1/fusions/*.tsv transabyss
    --convert ta_indels trans-abyss_TG1/indels/events_exons_novel.tsv transabyss

.. note::

    For older versions of MAVIS the convert command may require the path to the file(s) be quoted and the strandedness be specified (default is False)


3. **Which events you should validate in which libraries**

For this example, because we want to determine which events are germline/somatic we are going to pass all genome 
calls to both genomes. We can use either full file paths (if the input is already in the standard format) 
or the alias from a conversion (the first argument given to the convert option)

.. code:: text
    
    --assign TT1 ta_indels ta_trans chimera defuse
    --assign TG1 ta_indels ta_genome manta delly
    --assign NG1 ta_indels ta_genome manta delly

Putting this altogether with a name to call the config, we have the command to generate the pipeline config

.. code:: bash

    mavis config \
    --library NG1 genome normal False NG1.bam \
    --library TG1 genome diseased False TG1.bam \
    --library TT1 transcriptome diseased True TT1.bam \
    --convert ta_trans trans-abyss_TT1/fusions/*.tsv transabyss \
    --convert defuse defuse_TT1/results/results.classify.tsv defuse \
    --convert chimera chimeracan_TT1/chimeras.bedpe chimerascan \
    --convert delly delly_TG1_NG1/combined.vcf delly \
    --convert manta manta_TG1_NG1/results/variant/{Diploid,Somatic}.vcf \
    --convert ta_genome trans-abyss_TG1/fusions/*.tsv transabyss \
    --convert ta_indels trans-abyss_TG1/indels/events_exons_novel.tsv trans \
    --assign TT1 ta_indels ta_trans chimera defuse \
    --assign TG1 ta_indels ta_genome manta delly \
    --assign NG1 ta_indels ta_genome manta delly \
    -w mavis.cfg

.. note::

    If you haven't set defaults for the reference input files (using environment variables) you may need to pass those arguments as well.


Setting Up the Pipeline
.........................

The next step is :ref:`running the pipeline stage <pipeline-standard>`. This will perform conversion, clustering, and creating the :term:`SLURM`/:term:`SGE`
submission scripts for the other stages.

.. code:: bash

    mavis pipeline mavis.config -o output_dir/

Once complete this will look something like this. For simplicity only two of the job directories/files are shown (the number of jobs in configurable)
and only the contents of one of the library directories

.. code:: text

    output_dir/
    |-- converted_inputs/
    |   |-- delly.tab
    |   |-- chimera.tab
    |   |-- defuse.tab
    |   |-- manta.tab
    |   |-- ta_genome.tab
    |   |-- ta_indels.tab
    |   `-- ta_trans.tab
    |-- NG1_normal_genome/
    |   |-- cluster/
    |   |   |-- batch-d8zWDrpkbBuxj7eTahuH3K-1.tab
    |   |   |-- batch-d8zWDrpkbBuxj7eTahuH3K-2.tab
    |   |   |-- cluster_assignment.tab
    |   |   |-- clusters.bed
    |   |   |-- filtered_pairs.tab
    |   |   `-- MAVIS.COMPLETE
    |   |-- validate/
    |   |   |-- batch-d8zWDrpkbBuxj7eTahuH3K-1/
    |   |   |   `-- submit.sh
    |   |   `-- batch-d8zWDrpkbBuxj7eTahuH3K-2/
    |   |       `-- submit.sh
    |   `-- annotate/
    |   |   |-- batch-d8zWDrpkbBuxj7eTahuH3K-1/
    |   |   |   `-- submit.sh
    |   |   `-- batch-d8zWDrpkbBuxj7eTahuH3K-2/
    |   |       `-- submit.sh
    |-- TG1_diseased_genome/
    |-- TT1_diseased_genome/
    |-- pairing/
    |   `-- submit.sh
    |-- summary/
    |   `-- submit.sh
    `-- submit_pipeline_batch-d8zWDrpkbBuxj7eTahuH3K.sh

Submitting the Pipeline to SLURM
..................................

The last step is simple, ssh to your head node of your :term:`SLURM` cluster and run
the submit pipeline script. This contains the sbatch commands to submit all the individual 
jobs and chain their dependencies

.. code:: bash

    ssh head_node
    cd output_dir/
    bash submit_pipeline_batch-d8zWDrpkbBuxj7eTahuH3K.sh


