#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patrick.barry@noaa.gov
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink-trim_multiQC.out

source /home/pbarry/bin/hydraQC/bin/activate
multiqc /home/pbarry/AukePink_lcWGS/fastqc/trimmed/
