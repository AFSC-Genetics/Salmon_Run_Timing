#!/bin/bash
#SBATCH --job-name=bedtools_map
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load bio/bedtools/2.31.1

bedtools intersect -a $HOME/Mappability/Coho/GenMap_Res/Okis_v2_primary_only.genmap.bedgraph -b $HOME/Mappability/Coho/CohoGenome_ROI.bed -sorted > $HOME/Mappability/Coho/GenMap_Res/Coho_ROI_Mappability.bedgraph
