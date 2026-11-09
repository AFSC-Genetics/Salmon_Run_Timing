#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=esr1_ibs_fourspp
#SBATCH --output=/home/nhowe/runtiming/fourspecies/job_outfiles/fourspp_top5_esr1_IBS_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

gene=esr1
chrom=NC_068428.1
firstpos=23264160
lastpos=23645630

angsd -nThreads 10 -b /home/nhowe/runtiming/fourspecies/scripts/fourspp_${gene}_top5_ibs_input.txt \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-r ${chrom}:${firstpos}-${lastpos} \
	-out /home/nhowe/runtiming/fourspecies/gls/fourspp_GL1_MM4_top5_${gene} \
	-uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 -minQ 20 \
	-GL 1 -doMaf 1 -doMajorMinor 4 -doCounts 1 \
	-minMaf 0.05 -SNP_pval 1e-10 \
	-doCov 1 -makeMatrix 1 -doIBS 1
