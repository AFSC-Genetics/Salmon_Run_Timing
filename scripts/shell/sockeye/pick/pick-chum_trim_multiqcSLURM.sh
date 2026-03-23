#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum-trim_multiQC.out

source /home/nhowe/venv/bin/activate
multiqc /scratch2/nhowe/runtiming/sock-pick/fastqc/trimmed/
