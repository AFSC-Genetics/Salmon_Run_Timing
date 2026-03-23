#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_GCF_002021735.2_Okis_V2_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/bwa-index_GCF_002021735.2_Okis_V2_genomic.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/nhowe/runtiming/tocoho/bwa/GCF_002021735.2_Okis_V2_genomic /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna
