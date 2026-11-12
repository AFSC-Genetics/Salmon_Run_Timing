#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --job-name=pca_AukePink
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_pcangsd_%A-%a.out
#SBATCH --array=1-27%24

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

JOBS_FILE=/home/pbarry/AukePink_lcWGS/scripts/AukePink_pcangsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	beagle_file=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

chrom=$(echo $beagle_file | sed 's!^.*/!!')
chrom=${chrom%.beagle.gz}

pcangsd.py -threads 10 -beagle ${beagle_file} -o /home/pbarry/AukePink_lcWGS/pca/${chrom} -sites_save -pcadapt