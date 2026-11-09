#!/bin/bash
#SBATCH --job-name=filter_genome
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

grep "^>" $HOME/ref_genomes/GCF_002021735.2_Okis_V2_genomic.fna | head -n 30 | awk '{print $1}' | sed 's/>//' > $HOME/Mappability/Coho/Coho_ChrNames.txt

module load bio/samtools/1.19

samtools faidx $HOME/ref_genomes/GCF_002021735.2_Okis_V2_genomic.fna
samtools faidx $HOME/ref_genomes/GCF_002021735.2_Okis_V2_genomic.fna -r $HOME/Mappability/Coho/Coho_ChrNames.txt > $HOME/ref_genomes/Okis_V2_primary_only.fna
