#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_countdepth_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard,himem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

chrom=NC_034180.2

angsd -b /home/nhowe/runtiming/tocoho/tocoho_early_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/depth/tocoho_${chrom}_early -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -dumpCounts 3 -setminDepth 45 -setmaxDepth 900 -GL 1 -doMaf 1 -doMajorMinor 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/tocoho/tocoho_late_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/depth/tocoho_${chrom}_late -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -dumpCounts 3 -setminDepth 38 -setmaxDepth 760 -GL 1 -doMaf 1 -doMajorMinor 1 -only_proper_pairs 1
