#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=coho-chum_wgp-pca
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/coho-chum_wholegenome_polymorphic_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py  -threads 10  -beagle /home/nhowe/runtiming/coho/gls/coho-chum_wholegenome_polymorphic.beagle.gz -o /home/nhowe/runtiming/coho/pca/coho-chum_wholegenome-polymorphic -pcadapt