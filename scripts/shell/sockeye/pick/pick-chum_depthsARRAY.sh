#!/bin/bash

#SBATCH --job-name=depth_pick
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_depths_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=0-12:00:00
#SBATCH --partition=standard,medmem
#SBATCH --array=1-96%96

JOBS_FILE=/home/nhowe/runtiming/sock-pick/scripts/pick-chum_depthsARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	depth_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

touch /home/nhowe/runtiming/sock-pick/bamtools/pick-chum_depths.csv
python3 mean_cov_ind.py -i ${depth_file} -o /home/nhowe/runtiming/sock-pick/bamtools/pick-chum_depths.csv
