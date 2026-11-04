#!/bin/bash
#SBATCH --job-name=bedtools_map
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load bio/bedtools/2.31.1

bedtools intersect -a $HOME/ChumMap/Oket_V2_primary_only.genmap.bedgraph -b $HOME/ChumGenome_ROI.bed -sorted > $HOME/ChumMap/Chum_ROI_Mappability.bedgraph
