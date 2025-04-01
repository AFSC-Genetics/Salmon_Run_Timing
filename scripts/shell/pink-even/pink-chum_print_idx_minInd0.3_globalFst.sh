#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --job-name=idx_pink-chum
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chum_chrom_idxPrint_minInd0.3_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=standard,medmem

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

#minInd 0.3*N
#minMapQ = 15

# print idx FST file for each chromosome

for idx_name in /home/nhowe/runtiming/pink/fst/pink-chum_NC_*_early-late_minInd0.3.sfs.pbs.fst.idx 
do
	echo ${idx_name}
	realSFS fst print ${idx_name} > ${idx_name}.txt
done

