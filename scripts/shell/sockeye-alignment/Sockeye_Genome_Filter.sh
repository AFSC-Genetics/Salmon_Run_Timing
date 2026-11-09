#!/bin/bash
#SBATCH --job-name=filter_genome
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

grep "^>" $HOME/ref_genomes/GCF_006149115.2_Oner_1.1_genomic.fna | head -n 29 | awk '{print $1}' | sed 's/>//' > $HOME/Mappability/Sockeye/Sockeye_ChrNames.txt

module load bio/samtools/1.19

samtools faidx $HOME/ref_genomes/GCF_006149115.2_Oner_1.1_genomic.fna
samtools faidx $HOME/ref_genomes/GCF_006149115.2_Oner_1.1_genomic.fna -r $HOME/Mappability/Sockeye/Sockeye_ChrNames.txt > $HOME/ref_genomes/Oner_1.1_primary_only.fna
