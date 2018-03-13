.. _pipeline:

Running the Pipeline
-----------------------

Getting Help
................

All steps in the MAVIS pipeline are called following the main mavis entry point. The usage menu can be viewed
by running without any arguments, or by giving the -h/--help option

**Example:**

.. code:: bash

    mavis -h


Help sub-menus can be found by giving the pipeline step followed by no arguments or the -h options

.. code:: bash

    mavis cluster -h


Running MAVIS using a Job Scheduler
.........................................

The default setup and main 'pipeline' step of MAVIS is set up to use a job scheduler on a compute cluster. Two schedulers are currently
supported: :term:`SLURM` and :term:`SGE`. Using the pipeline step
will generate submission scripts and a wrapper bash script for the user to execute on their cluster head node.

.. figure:: _static/pipeline_options.svg
    :width: 100%

    The MAVIS pipeline is highly configurable. Some pipeline steps (cluster, validate) are optional and can be automatically skipped.
    The standard pipeline is far-left.

.. _pipeline-standard:

Standard
+++++++++++

The most common use case is :ref:`auto-generating a configuration file <pipeline-config>` and then running the pipeline setup step.
The pipeline setup step will run clustering and create scripts for running the other steps.

.. code:: bash

    mavis config .... -w config.cfg
    mavis pipeline config.cfg -o /path/to/top/output_dir

This will create submission scripts as follows

.. code:: text

    output_dir/
    |-- library1/
    |   |-- validate/<jobdir>/submit.sh
    |   `-- annotate/<jobdir>/submit.sh
    |-- library2/
    |   |-- validate/<jobdir>/submit.sh
    |   `-- annotate/<jobdir>/submit.sh
    |-- pairing/submit.sh
    |-- summary/submit.sh
    `-- submit_pipeline_<batchid>.sh

The submit_pipeline_<batchid>.sh is the wrapper script which can be executed on the head node

.. code:: bash

    ssh cluster_head_node
    cd /path/to/output_dir
    bash submit_pipeline_<batchid>.sh

This will submit a series of jobs with dependencies. 

.. _pipeline-dependency-graph:


.. figure:: _static/pipeline_dependency_graph.svg
    :width: 100%

    Dependency graph of MAVIS jobs for the standard pipeline setup. The notation on the arrows indicates the 
    SLURM setting on the job to add the dependency on the previous job.


Non-Standard
+++++++++++++++

To set up a non-standard pipeline and skip steps use the skip stage option.

.. code:: bash

    mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage cluster

.. code:: bash

    mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage validate

Or to skip both clustering and validation, simply call the option twice.

.. code:: bash

    mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage cluster --skip_stage validate

.. note::

    skipping clustering will still produce and output directory and files, but no merging will be done

Configuring Scheduler Settings
+++++++++++++++++++++++++++++++

There are multiple ways to configure the scheduler settings. Some of the configurable options are listed below

- :term:`queue` ``MAVIS_QUEUE``
- :term:`memory_limit` ``MAVIS_MEMORY_LIMIT``
- :term:`time_limit` ``MAVIS_TIME_LIMIT``
- :term:`import_env` ``MAVIS_IMPORT_ENV``
- :term:`scheduler` ``MAVIS_SCHEDULER``

For example to set the job queue default using an :ref:`environment variable <config-environment>`

.. code:: bash

    export MAVIS_QUEUE=QUEUENAME

Or to give it as an argument during :ref:`config generation <pipeline-config>`

.. code:: bash

    mavis config -w /path/to/config --queue QUEUENAME

Finally it can also be added to the config file manually

.. code:: text

    [schedule]
    queue = QUEUENAME


Troubleshooting Dependency Failures
++++++++++++++++++++++++++++++++++++++

The most common error to occur when running MAVIS on the cluster is a memory or time limit exception. These can be detected by running the checker or looking for dependency failures reported on the cluster. The suffix of the job name will be a number and will correspond to the suffix of the job directory. Note that the following example commands are :term:`SLURM`-specific and do not apply to :term:`SGE`.

.. code::

    mavis checker -o /path/to/output/dir

This will report any failed jobs. For example if this were a memory issue for one of the validation jobs we might expect to see something like below in the checker output

.. code:: text

    validate FAIL
        1 jobs CRASHED (jobs: <job number>)
            slurmstepd: error: exceeded job memory limit (jobs: <job number>)

Each job has its own submission script. The values for memory/time limits or other parameters can be edited by editing the script header and resubmitting the failed job.

.. code:: bash
    
    cd /path/to/output/dir/<library>*/validate/*-<job number>
    vim submit.sh  # edit the header to change memory/time/etc.
    sbatch submit.sh

Now you will need to edit the annotation job which was dependent on it. For convenience it is easiest to find the job by name. The annotation job name will be the same as the validation job name except that instead of the ``MV_`` prefix it will begin with the ``MA_`` prefix.

.. code:: bash
    
    squeue -n MA_*-<job number>

And then change the dependency to be the new validation job

.. code:: bash
        
    scontrol update job=<annotation job id> Dependency=afterok:<new validation job id>

If memory errors are frequent then it would be better to adjust the default values (:term:`trans_validation_memory`, :term:`validation_memory`, :term:`time_limit`)


