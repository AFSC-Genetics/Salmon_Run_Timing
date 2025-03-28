#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/anvb-raw_multiQC_%j.out

source /home/nhowe/venv/bin/activate
multiqc /home/nhowe/runtiming/anvil/fastqc/raw/