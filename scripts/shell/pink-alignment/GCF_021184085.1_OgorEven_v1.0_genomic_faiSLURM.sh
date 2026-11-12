#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_GCF_021184085.1_OgorEven_v1.0_genomic
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/fai_GCF_021184085.1_OgorEven_v1.0_genomic.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/pbarry/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna
