#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=AukePink_wgp-pca-admix
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_wholegenome_polymorphic_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.beagle.gz -o /home/pbarry/AukePink_lcWGS/pca/AukePink_wholegenome-polymorphic -sites_save -snp_weights -pcadapt

