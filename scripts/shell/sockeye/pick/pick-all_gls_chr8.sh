#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=0-20:00:00
#SBATCH --job-name=gls_pick-all
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_gls_minInd0.3_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

chrom=NC_068428.1

angsd -nThreads 8 -b /home/nhowe/runtiming/sock-pick/sock-chum_plusPick_bamslist.txt \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-r ${chrom}: \
	-out /home/nhowe/runtiming/sock-pick/gls/pick-all_${chrom}_minInd0.3 \
	-only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 15 -minQ 20 \
	-doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
	-minMaf 0.05 -SNP_pval 1e-10 \
	-setMinDepth 97 -setMaxDepth 1940 -minInd 58
