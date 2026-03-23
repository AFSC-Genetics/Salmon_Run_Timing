#!/bin/bash

#SBATCH --job-name=align
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_alignment_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --exclude=node[01-08]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --time=3-00:00:00
#SBATCH --array=1-96%48

module unload aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5
module load aligners/bwa/0.7.17 bio/samtools/1.11 bio/bamtools/2.5.1 bio/picard/2.23.9 bio/bamutil/1.0.5

JOBS_FILE=/home/nhowe/runtiming/tocoho/scripts/tocoho_alignARRAY_input.txt
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

bwa mem -M -t 10 /home/nhowe/runtiming/tocoho/bwa/GCF_002021735.2_Okis_V2_genomic ${fq_r1} ${fq_r2} 2> /home/nhowe/runtiming/tocoho/bwa/tocoho_${sample_id}_bwa-mem.out > /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.sam

samtools view -bS -F 4 /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.sam > /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.bam
rm /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.sam

samtools view -h /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.bam | samtools view -buS - | samtools sort -o /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted.bam
rm /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.bam

java -jar $PICARD MarkDuplicates I=/home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted.bam O=/home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup.bam M=/home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_dups.log VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
rm /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted.bam

bam clipOverlap --in /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup.bam --out /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup_clipped.bam --stats
rm /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup.bam

samtools depth -aa /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup_clipped.bam | cut -f 3 | gzip > /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}.depth.gz

samtools index /home/nhowe/runtiming/tocoho/bamtools/tocoho_${sample_id}_sorted_dedup_clipped.bam
