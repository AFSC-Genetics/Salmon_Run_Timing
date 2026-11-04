#!/bin/bash
#SBATCH --job-name=genmap_index
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -p standard,medmem

eval "$(mamba shell hook --shell bash)"

mamba activate genmap-1.3.0

genmap index -F $HOME/ref_genomes/Oket_V2_primary_only.fna -I $HOME/GenMap
