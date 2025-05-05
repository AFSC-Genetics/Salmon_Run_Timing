#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=lrrc9_ibs_fourspp
#SBATCH --output=/home/nhowe/runtiming/fourspecies/job_outfiles/fourspp_top5_lrrc9_IBS_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

gene=lrrc9
chrom=NC_068455.1
firstpos=28128954
lastpos=28169980

angsd -b /home/nhowe/runtiming/fourspecies/scripts/fourspp_lrrc9_top5_ibs_input.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -r ${chrom}:${firstpos}-${lastpos} -out /home/nhowe/runtiming/fourspecies/gls/fourspp_top5_GL1_MM4_${gene} -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -GL 1 -doMaf 1 -doMajorMinor 4 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -minMaf 0.05 -SNP_pval 1e-10 -only_proper_pairs 1
