#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --job-name=idx_chumrun
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_chrom_idxPrint_minInd0.3_minDepthHalf_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=standard,medmem

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

#minInd 0.3_minDepthHalf*N
#minMapQ = 20

# print idx FST file for each chromosome

for idx_name in /home/nhowe/runtiming/chum/fst/chumrun_NC_*_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.idx 
do
	echo ${idx_name}
	realSFS fst print ${idx_name} > ${idx_name}.txt
done

