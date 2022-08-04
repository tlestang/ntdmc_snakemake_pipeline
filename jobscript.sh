#!/bin/bash

module load Python/3.7.4-GCCcore-8.3.0
module load snakemake/5.26.1-foss-2019b-Python-3.7.4
module load R/3.6.2-foss-2019b

# Print some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Job ID: $SGE_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# properties = {properties}

{exec_job}
