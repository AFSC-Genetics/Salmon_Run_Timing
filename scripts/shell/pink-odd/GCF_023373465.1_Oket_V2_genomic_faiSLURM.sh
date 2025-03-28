#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_GCF_023373465.1_Oket_V2_genomic
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/pinkOdd/job_outfiles/fai_GCF_023373465.1_Oket_V2_genomic.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna
