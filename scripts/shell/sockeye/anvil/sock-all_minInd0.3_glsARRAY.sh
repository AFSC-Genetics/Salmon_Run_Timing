#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=plm_sock-all
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/sock-all_minInd0.3_chr35_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

contig=NC_068455.1

angsd -b /home/nhowe/runtiming/anvil/all_sockeye_bamslist.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -r ${contig}: -out /home/nhowe/runtiming/anvil/gls/sock-all_${contig}_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -minQ 20 -doCounts 1 -minInd 30 -setminDepth 100 -setmaxDepth 2000 -doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-10 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1
