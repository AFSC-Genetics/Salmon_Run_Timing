#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mem=20G
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink-raw_multiQC.out

source /home/pbarry/bin/hydraQC/bin/activate
multiqc /home/pbarry/AukePink_lcWGS/fastqc/raw/
