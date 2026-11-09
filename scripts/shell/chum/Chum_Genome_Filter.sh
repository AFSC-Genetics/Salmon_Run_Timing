#!/bin/bash
#SBATCH --job-name=filter_genome
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

grep "^>" $HOME/ref_genomes/GCA_023373465.1_Oket_V2_genomic.fna | head -n 37 | awk '{print $1}' | sed 's/>//' > ChumChrNames.txt

module load bio/samtools/1.19

samtools faidx $HOME/ref_genomes/GCA_023373465.1_Oket_V2_genomic.fna
samtools faidx $HOME/ref_genomes/GCA_023373465.1_Oket_V2_genomic.fna -r ChumChrNames.txt > $HOME/ref_genomes/Oket_V2_primary_only.fna
