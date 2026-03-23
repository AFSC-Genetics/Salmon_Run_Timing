#!/bin/bash

#SBATCH --job-name=align
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_alignment_%A-%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=3-00:00:00
#SBATCH --partition=medmem,standard
#SBATCH --array=1-96%12

module purge
module load aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5

#mkdir -p /scratch2/nhowe/runtiming/sock-pick/bamtools

JOBS_FILE=/home/nhowe/runtiming/sock-pick/scripts/pick-chum_alignARRAY_input.txt
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

echo "1) bwa mem -------------------------------------------------------------"
bwa mem -M -t 10 /home/nhowe/runtiming/sock-pick/bwa/GCF_023373465.1_Oket_V2_genomic ${fq_r1} ${fq_r2} 2> /home/nhowe/runtiming/sock-pick/bwa/pick-chum_${sample_id}_bwa-mem.out > /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.sam

echo "2) samtools view -bS -F 4 ----------------------------------------------"
samtools view -bS -F 4 /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.sam > /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.bam
rm /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.sam

echo "3) sort: samtools view -h ----------------------------------------------"
samtools view -h /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.bam | samtools view -buS - | samtools sort -o /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted.bam
rm /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.bam

echo "4) java mark duplicates ----------------------------------------------"
java -jar $PICARD MarkDuplicates I=/scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted.bam O=/scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup.bam M=/scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
rm /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted.bam

echo "5) bam clipoverlap ----------------------------------------------------"
bam clipOverlap --in /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup.bam --out /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup_clipped.bam --stats
rm /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup.bam

echo "2) samtools depth -----------------------------------------------------"
samtools depth -aa /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup_clipped.bam | cut -f 3 | gzip > /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}.depth.gz

echo "2) samtools index -----------------------------------------------------"
samtools index /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_${sample_id}_sorted_dedup_clipped.bam

echo " -------- DONE -----------------------------------------------------------"
