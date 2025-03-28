#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chum-trim_multiQC.out

source /home/nhowe/venv/bin/activate
multiqc /home/nhowe/runtiming/pink/fastqc/trimmed/