#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_GCF_002021735.2_Okis_V2_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/fai_GCF_002021735.2_Okis_V2_genomic.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna
