#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_GCF_023373465.1_Oket_V2_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/bwa-index_GCF_023373465.1_Oket_V2_genomic.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/nhowe/runtiming/coho/bwa/GCF_023373465.1_Oket_V2_genomic /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna
