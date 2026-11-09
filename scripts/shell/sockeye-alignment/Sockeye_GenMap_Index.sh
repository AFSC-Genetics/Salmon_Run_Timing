#!/bin/bash
#SBATCH --job-name=genmap_index
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH -p standard,medmem

eval "$(mamba shell hook --shell bash)"

mamba activate genmap-1.3.0

genmap index -F $HOME/ref_genomes/Oner_1.1_primary_only.fna -I $HOME/Mappability/Sockeye/GenMap_Index
