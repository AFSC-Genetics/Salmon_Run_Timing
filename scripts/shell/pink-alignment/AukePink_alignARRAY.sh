#!/bin/bash

#SBATCH --job-name=align
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_alignment_%A-%a.out
#SBATCH --time=7-00:00:00
#SBATCH --array=1-288%12

module unload aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5
module load aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5

JOBS_FILE=/home/pbarry/AukePink_lcWGS/scripts/AukePink_alignARRAY_input.txt
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

bwa mem -M -t 10 /home/pbarry/AukePink_lcWGS/bwa/GCF_021184085.1_OgorEven_v1.0_genomic ${fq_r1} ${fq_r2} 2> /home/pbarry/AukePink_lcWGS/bwa/AukePink_${sample_id}_bwa-mem.out > /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.sam

samtools view -bS -F 4 /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.sam > /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.bam
rm /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.sam

samtools view -h /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.bam | samtools view -buS - | samtools sort -o /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted.bam
rm /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.bam

java -jar $PICARD MarkDuplicates I=/home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted.bam O=/home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup.bam M=/home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
rm /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted.bam

bam clipOverlap --in /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup.bam --out /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup_clipped.bam --stats
rm /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup.bam

samtools depth -aa /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup_clipped.bam | cut -f 3 | gzip > /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}.depth.gz

samtools index /home/pbarry/AukePink_lcWGS/bamtools/AukePink_${sample_id}_sorted_dedup_clipped.bam