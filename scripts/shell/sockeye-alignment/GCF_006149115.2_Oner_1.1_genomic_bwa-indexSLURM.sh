#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_GCF_006149115.2_Oner_1.1_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/maturity/euclide/job_outfiles/bwa-index_GCF_006149115.2_Oner_1.1_genomic.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/nhowe/maturity/euclide/bwa/GCF_006149115.2_Oner_1.1_genomic /home/nhowe/reference_genomes/sockeye/GCF_006149115.2_Oner_1.1_genomic.fna
