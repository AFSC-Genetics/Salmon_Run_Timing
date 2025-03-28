#!/bin/bash

#SBATCH --cpus-per-task=5
#SBATCH --time=0-12:00:00
#SBATCH --job-name=ds_pink-chum
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chum_downsample-bams_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --exclude=node[29,30]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-13%13

module unload bio/samtools/1.11
module load bio/samtools/1.11

JOBS_FILE=/home/nhowe/runtiming/pink/scripts/pink-chum_13indivs_downsampleARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	sample_id=$(echo ${sample_line} | awk -F ":" '{print $2}')
	downsample_value=$(echo ${sample_line} | awk -F ":" '{print $3}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

samtools view -bo /home/nhowe/runtiming/pink/bamtools/pink-chum_${sample_id}_sorted_dedup_clipped_downsampled.bam -s ${downsample_value} /home/nhowe/runtiming/pink/bamtools/pink-chum_${sample_id}_sorted_dedup_clipped.bam
samtools depth -aa /home/nhowe/runtiming/pink/bamtools/pink-chum_${sample_id}_sorted_dedup_clipped_downsampled.bam | cut -f 3 | gzip > /home/nhowe/runtiming/pink/bamtools/pink-chum_${sample_id}_sorted_dedup_clipped_downsampled.depth.gz
samtools index /home/nhowe/runtiming/pink/bamtools/pink-chum_${sample_id}_sorted_dedup_clipped_downsampled.bam
