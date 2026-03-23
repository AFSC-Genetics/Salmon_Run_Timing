#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fstidx_pick-chum
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_fst_analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=standard,medmem
#SBATCH --exclude=node[30]
#SBATCH --array=1-37%18

module purge
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/sock-pick/scripts/pick-chum_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

#realSFS /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_early_minInd0.3.saf.idx -fold 1 \
#	/home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_late_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs

#realSFS fst index /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_early_minInd0.3.saf.idx -fold 1 \
#	/home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_late_minInd0.3.saf.idx -fold 1 \
#	-sfs /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs \
#	-fstout /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs -whichFst 1

#realSFS fst stats2 /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.txt

realSFS fst print /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.idx > /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.idx.txt
