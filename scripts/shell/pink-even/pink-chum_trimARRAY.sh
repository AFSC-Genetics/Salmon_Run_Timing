#!/bin/bash

#SBATCH --job-name=trim
#SBATCH --cpus-per-task=4
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chum_trimming_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --exclude=node[30,31]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=0-12:00:00
#SBATCH --array=1-154%48

module unload bio/trimmomatic/0.39 bio/fastp/0.23.2
module load bio/trimmomatic/0.39 bio/fastp/0.23.2
JOBS_FILE=/home/nhowe/runtiming/pink/scripts/pink-chum_trimARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	fq_r1=$(echo ${sample_line} | awk -F ":" '{print $2}')
	fq_r2=$(echo ${sample_line} | awk -F ":" '{print $3}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

sample_id=$(echo $fq_r1 | sed 's!^.*/!!')
sample_id=${sample_id%%_*}

java -jar ${TRIMMOMATIC} PE -threads 4 -phred33 ${fq_r1} ${fq_r2} /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R1_paired.fq.gz /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R1_unpaired.fq.gz /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R2_paired.fq.gz /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R2_unpaired.fq.gz ILLUMINACLIP:/home/nhowe/reference_genomes/NexteraPE-PE.fa:2:30:10:1:true MINLEN:40
fastp --trim_poly_g -L -A --cut_right -i /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R1_paired.fq.gz -o /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_clipped_R1_paired.fq.gz -I /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_R2_paired.fq.gz -O /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_clipped_R2_paired.fq.gz -h /home/nhowe/runtiming/pink/trimmed/${sample_id}_trimmed_clipped_paired_report.html
