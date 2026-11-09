#!/bin/bash
#SBATCH --job-name=bedtools_map
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load bio/bedtools/2.31.1

bedtools intersect -a $HOME/Mappability/Sockeye/GenMap_Res/Oner_1.1_primary_only.genmap.bedgraph -b $HOME/Mappability/Sockeye/SockeyeGenome_ROI.bed -sorted > $HOME/Mappability/Sockeye/GenMap_Res/Sockeye_ROI_Mappability.bedgraph
