#!/bin/bash

#SBATCH -p medmem
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=AukePink_wgp-pca-admix
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_wholegenome_polymorphicEVEN_%A.out

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

pcangsd.py -threads 10 -beagle /home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.beagle.gz -o /home/pbarry/AukePink_lcWGS/pca/Even/AukePink_wholegenome-polymorphic -sites_save -filter ./EvenSamps_Bam.txt -snp_weights -pcadapt

