#!/bin/bash
#SBATCH --job-name=samtools_PinkRef
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load bio/samtools

samtools faidx $HOME/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna NC_060182.1:72673514-72737012 > $HOME/AukePink_lcWGS/PinkLG10ref.fa

