#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/pinkOdd/job_outfiles/pink-odd-trim_multiQC.out

source /home/nhowe/venv/bin/activate
multiqc /home/nhowe/runtiming/pinkOdd/fastqc/trimmed/