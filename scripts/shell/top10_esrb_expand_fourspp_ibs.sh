#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=esrb_ibs_fourspp
#SBATCH --output=/home/nhowe/runtiming/fourspecies/job_outfiles/fourspp_top10_esrb_IBS_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

gene=esrb_expand
chrom=NC_068449.1
firstpos=25414060 #expanded region based on boundary SNPs within peak that had fst > 0.9
lastpos=25501622 #expanded region based on boundary SNPs within peak that had fst > 0.9

angsd -b /home/nhowe/runtiming/fourspecies/scripts/twospp_esrb_top10_ibs_input.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -r ${chrom}:${firstpos}-${lastpos} -out /home/nhowe/runtiming/fourspecies/gls/twospp_top10_${gene} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -GL 2 -doMaf 1 -doMajorMinor 1 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -minMaf 0.05 -SNP_pval 1e-10 -only_proper_pairs 1
