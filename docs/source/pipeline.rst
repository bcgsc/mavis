
Running the Pipeline
-----------------------

Getting Help
................

All steps in the MAVIS pipeline are called following the main mavis entry point. The usage menu can be viewed
by running without any arguments, or by giving the -h/--help option

**Example:**

.. code-block:: bash

    >>> mavis -h


Help sub-menus can be found by giving the pipeline step followed by no arguments or the -h options

.. code-block:: bash

    >>> mavis cluster -h


Pipeline Use Cases
.......................

.. figure:: _static/pipeline_options.svg
    :width: 100%

    The MAVIS pipeline is highly configurable. Some pipeline steps (cluster, validate) are optional and can be automatically skipped. 
    The standard pipeline is far-left.


Standard
+++++++++++

The most common use case is auto-generating a configuration file and then running the pipeline setup step.
The pipeline setup step will run clustering and create scripts for running the other steps.

.. code-block:: bash

    >>> mavis config .... -w config.cfg
    >>> mavis pipeline config.cfg -o /path/to/top/output_dir

This will create submission scripts as follows

.. code:: text

    output_dir/
    |-- library1/
    |   |-- validation/<jobdir>/submit.sh
    |   `-- annotation/<jobdir>/submit.sh
    |-- library2/
    |   |-- validation/<jobdir>/submit.sh
    |   `-- annotation/<jobdir>/submit.sh
    |-- pairing/submit.sh
    |-- summary/submit.sh
    `-- submit_pipeline_<batchid>.sh

The scripts are meant for submission to an `SGE <http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html>`_
or `SLURM <https://slurm.schedmd.com>`_ cluster. The summary job is held on the pairing job, the pairing job is held on all the annotation jobs,
and the annotation jobs are held on their validation jobs. To submit all jobs simply bash the submit pipeline script in the top-level
directory

.. code-block:: bash

    >>> ssh cluster_head_node
    >>> cd /path/to/output_dir
    >>> bash submit_pipeline_<batchid>.sh


Non-Standard
+++++++++++++++

To set up a non-standard pipeline and skip steps use the skip stage option.

.. code:: bash

    >>> mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage cluster

.. code:: bash

    >>> mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage validate

Or to skip both clustering and validation, simply call the option twice.

.. code:: bash

    >>> mavis pipeline /path/to/config -o /path/to/output/dir --skip_stage cluster --skip_stage validate

.. note::

    skipping clustering will still produce and output directory and files, but no merging will be done


