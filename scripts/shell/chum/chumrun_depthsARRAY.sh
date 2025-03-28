#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --cpus-per-task=5
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrundepths_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=medmem
#SBATCH --time=0-12:00:00
#SBATCH --array=1-96%48

JOBS_FILE=/home/nhowe/runtiming/chum/scripts/chumrun_depthsARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	depth_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

touch /home/nhowe/runtiming/chum/bamtools/chumrun_depths.csv
python3 mean_cov_ind.py -i ${depth_file} -o /home/nhowe/runtiming/chum/bamtools/chumrun_depths.csv
