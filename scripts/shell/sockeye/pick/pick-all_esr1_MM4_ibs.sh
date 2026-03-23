#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=er1_ibs_sockeye
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-all_esr1_IBS_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

gene=esr1
chrom=NC_068428.1
firstpos=23264160
lastpos=23645630

angsd -b /home/nhowe/runtiming/sock-pick/sock-chum_plusPick_bamslist.txt \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-r ${chrom}:${firstpos}-${lastpos} \
	-out /home/nhowe/runtiming/sock-pick/gls/pick-all_${gene}_MM4 \
	-nThreads 10 -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 15 -minQ 20 \
	-GL 1 -doMaf 1 -doMajorMinor 4 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 \
	-minMaf 0.05 -SNP_pval 1e-10
