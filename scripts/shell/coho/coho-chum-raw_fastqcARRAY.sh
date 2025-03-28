#!/bin/bash

#SBATCH --job-name=fqc_array_coho-chum
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/coho-chum-raw_fastqc_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=standard,medmem
#SBATCH --time=0-03:00:00
#SBATCH --array=1-192%48

module unload bio/fastqc/0.11.9
module load bio/fastqc/0.11.9

JOBS_FILE=/home/nhowe/runtiming/coho/scripts/coho-chum-raw_fqcARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

fastqc ${fq} -o /home/nhowe/runtiming/coho/fastqc/raw/
