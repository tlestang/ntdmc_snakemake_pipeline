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


Execution
---------

The pipeline is executed by submitting a snakemake job to the queuing
system.  You can think of this job as the "puppeteer": it is
responsible for monitoring jobs, output files, and submiting rule
jobs.

To submit the snakemake job, start by writing the submission script:

.. code-block:: shell

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

   snakemake --profile profiles/BMRC --jobscript=jobscript.sh \
       --group-components resimulate=4

       


You can save this script as e.g. `snakemake-job.sh`. Lastly, you can
submit the smakemake job with:

.. code-block:: shell

   sbatch -A hollingsworth.prj -p short snakemake-job.sh



   


