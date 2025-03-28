#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=depth_check
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_depthcheck_angsd_%A.out
#SBATCH --error=/home/nhowe/runtiming/chum/job_outfiles/chumrun_depthcheck_%A.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

angsd -b /home/nhowe/runtiming/chum/chumrun_bamslist.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/gls/chumrun_depthcheck -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -doCounts 1 -GL 1 -doDepth 1 -only_proper_pairs 1
