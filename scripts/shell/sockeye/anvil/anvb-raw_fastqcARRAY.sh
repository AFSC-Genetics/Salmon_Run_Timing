#!/bin/bash

#SBATCH --job-name=fqc_array_anvb
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/anvb-raw_fastqc_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=0-03:00:00
#SBATCH --array=1-54%34

module unload bio/fastqc/0.11.9
module load bio/fastqc/0.11.9

JOBS_FILE=/home/nhowe/runtiming/anvil/scripts/anvb-raw_fqcARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

fastqc ${fq} -o /home/nhowe/runtiming/anvil/fastqc/raw/
