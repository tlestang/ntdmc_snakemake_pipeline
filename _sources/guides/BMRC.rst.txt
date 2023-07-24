Oxford BMRC
===========

This document will guide you across the set up and execution fo the
NTDMC inference pipeline on the BMRC cluster.

Set up
------

Start by downloading the workflow and supporting files.  You can do
this with Git

.. code-block:: shell

   git clone https://github.com/tlestang/ntdmc_snakemake_pipeline.git


You can also download the reposiroty as an archive from the GitHub repository page.

Most of the workflow rules are is executed as cluster jobs. Snakemake
automatically handles job status monitoring and submission. In order
to ensure the job submission script submitted by Snakemake loads the
right modules, we must specify a job submission template script.

Copy the job submission template script into the pipeline directory:

.. code-block:: shell

   cp /well/hollingsworth/shared/inference_pipeline/BMRC_jobscript.sh .

Input data
----------

The pipeline repository comes ready with a small dataset
`data/quickstart.csv` covering a few locations in Ethiopia, Uganda,
Senegal, Burkina Faso, Nigeria and Tanzania.

A larger dataset is available in the shared directory:

.. code-block::

   /well/hollingsworth/shared/inference_pipeline/data/diggle.csv

Configuration
-------------

There are two types of configuration we need to take care of:

1. Parameters relative to where the workflow is running, in this case
   the BMRC (e.g. job submission command). These parameters are
   unlikely to change across runs.
2. Parameters of the workflow itself (e.g. number of disease parameter
   values to forward simulate). These parameters are likely to change
   across runs.

Cluster parameters
~~~~~~~~~~~~~~~~~~

Cluster parameters could be specified manually as command line options
when invoking `snakemake`.  Because their value is not expected to
change much, let's group them in a configuration file instead.  This
file is referred to as a cluster *profile* is snakemake jargon.

The pipeline repository already contains a default profile
`profiles/default/config.yml`.  Let's make a copy of it, creating a
new BMRC profile:

.. code-block:: shell

   mkdir profiles/BMRC
   cp profiles/default/config.yml profiles/BMRC

The default profile is almost ready to go.  The only line to change is
the submission command (`cluster` entry).

.. code-block:: diff

   --- profiles/default/config.yaml
   +++ profiles/BMRC/config.yaml
   @@ -1,4 +1,4 @@
   -cluster: "sbatch"
   +cluster: "sbatch -A hollingsworth.prj -p short"
   jobs: 10
   latency-wait: 60
   max-jobs-per-second: 1

This ensure rule jobs will be attached to the `hollingsworth` project
allocation, and submitted to the `short` queue.

Workflow parameters
~~~~~~~~~~~~~~~~~~~

Similarly to cluster parameters, values for workflow parameters can be
specified in a YAML configuration file.  The pipeline repository comes
with a default workflow configuration file `config/default.yml` that
should be good enough to start with.

.. note::

   It is also possible to specify value for workflow parameters
   as command line options to `snakemake`:

   .. code-block:: shell

      snakemake --config <keyword>=<value>

   Values specified at the command line take precedence over values
   specified in a configuration file.


Execution
---------

The pipeline is executed by submitting a main snakemake job to the
queuing system.  You can think of this job as a puppeteer: it is
responsible for submitting jobs for rules, monitoring their status and
checking their output.

Before we submit the main snakemake job to the queuing system, we must
specify a template jobscript for rule job submission.  You can copy
the `BMRC_jobscript.sh` from the shared directory:

.. code-block:: shell

   cp /well/hollingsworth/shared/inference_pipeline/BMRC_jobscript.sh .

We are now ready to submit the main snakemake job.  Let's start by
writing the submission script:

.. code-block:: shell

   # snakemake-job.sh
   #!/bin/bash

   #SBATCH -J mainsnake
   #SBATCH -o mainsnake-%j.out
   #SBATCH -e mainsnake-%j.err

   echo "------------------------------------------------"
   echo "Run on host: "`hostname`
   echo "Operating system: "`uname -s`
   echo "Username: "`whoami`
   echo "Started at: "`date`
   echo "------------------------------------------------"

   module load snakemake/7.22.0-foss-2022a
   snakemake --profile profiles/BMRC \
             --configfile configs/default.yml \
             --jobscript=BMRC_jobscript.sh \
             --group-components resimulate=4

To save some typing, you can also find this script in the shared
directory:

.. code-block:: shell

   cp /well/hollingsworth/shared/inference_pipeline/snakemake-job.sh .

Lastly, you can submit the main smakemake job with:

.. code-block:: shell

   sbatch -A hollingsworth.prj -p short snakemake-job.sh

.. code-block::

   user@rescomp1$ squeue -u user
             JOBID PARTITION     NAME
          22886203     short snakejob
          22886204     short snakejob
          22886205     short snakejob
          22886206     short snakejob
          22886182     short mainsnak
