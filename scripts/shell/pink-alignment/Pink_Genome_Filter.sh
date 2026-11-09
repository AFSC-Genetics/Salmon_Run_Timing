#!/bin/bash
#SBATCH --job-name=filter_genome
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

grep "^>" $HOME/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna | head -n 27 | awk '{print $1}' | sed 's/>//' > $HOME/Mappability/Pink/Pink_ChrNames.txt

module load bio/samtools/1.19

samtools faidx $HOME/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna
samtools faidx $HOME/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna -r $HOME/Mappability/Pink/Pink_ChrNames.txt > $HOME/ref_genomes/OgorEven_v1.0_primary_only.fna
