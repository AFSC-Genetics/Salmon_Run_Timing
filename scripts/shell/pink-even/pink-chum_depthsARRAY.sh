#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --cpus-per-task=5
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chumdepths_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --exclude=node[29,30]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=0-12:00:00
#SBATCH --array=1-154%24

JOBS_FILE=/home/nhowe/runtiming/pink/scripts/pink-chum_depthsARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	depth_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

touch /home/nhowe/runtiming/pink/bamtools/pink-chum_depths.csv
python3 mean_cov_ind.py -i ${depth_file} -o /home/nhowe/runtiming/pink/bamtools/pink-chum_depths.csv
