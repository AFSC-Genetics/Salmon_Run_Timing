#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_GCF_021184085.1_OgorEven_v1.0_genomic
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/bwa-index_GCF_021184085.1_OgorEven_v1.0_genomic.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/pbarry/AukePink_lcWGS/bwa/GCF_021184085.1_OgorEven_v1.0_genomic /home/pbarry/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna
