#!/bin/bash
#SBATCH --job-name=bedtools_map
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load bio/bedtools/2.31.1

bedtools intersect -a $HOME/Mappability/Chum/GenMap_Res/Oket_V2_primary_only.genmap.bedgraph -b $HOME/Mappability/Chum/ChumGenome_ROI.bed -sorted > $HOME/Mappability/Chum/GenMap_Res/Chum_ROI_Mappability.bedgraph
