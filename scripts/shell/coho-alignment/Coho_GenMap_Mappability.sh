#!/bin/bash
#SBATCH --job-name=genmap_map
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -p standard,medmem

eval "$(mamba shell hook --shell bash)"

mamba activate genmap-1.3.0

genmap map -K 150 -E 2 -T 16 -I $HOME/Mappability/Coho/GenMap_Index -O $HOME/Mappability/Coho/GenMap_Res -bg
