#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/maturity/euclide/job_outfiles/euclide-trim_multiQC.out
#SBATCH --error=/home/nhowe/maturity/euclide/job_outfiles/euclide-trim_multiQC.err
#SBATCH --mem=10G

source /home/nhowe/venv/bin/activate
multiqc /home/nhowe/maturity/euclide/fastqc/trimmed/
