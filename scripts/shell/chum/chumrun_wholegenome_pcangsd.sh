#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-05:00:00
#SBATCH --job-name=chumrun_wgp-pca
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=medmem
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_wholegenome_polymorphic_%A.out
#SBATCH --error=/home/nhowe/runtiming/chum/job_outfiles/chumrun_wholegenome_polymorphic_%A.err

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py  -threads 10  -beagle /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_polymorphic.beagle.gz -o /home/nhowe/runtiming/chum/pca/chumrun_wholegenome-polymorphic #-pcadapt
