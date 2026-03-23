#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_GCF_006149115.2_Oner_1.1_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/maturity/euclide/job_outfiles/fai_GCF_006149115.2_Oner_1.1_genomic.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/nhowe/reference_genomes/sockeye/GCF_006149115.2_Oner_1.1_genomic.fna
