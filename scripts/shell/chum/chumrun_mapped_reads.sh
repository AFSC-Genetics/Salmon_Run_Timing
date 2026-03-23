#!/bin/bash

#SBATCH --job-name=chumReads
#SBATCH --cpus-per-task=10
#SBATCH --time=2-24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/atka/job_outfiles/atka_mapped_reads.%j.out # update your out file directory

module unload bio/samtools/1.19
module load bio/samtools/1.19

cd /home/nhowe/atka/bamtools

touch atka_mapped_reads_concatenate.tsv

for file in atka_*_mapped_reads.tsv;
do
	sample=result=$(echo "$file" | awk -F '_' '{print $2}')
	reads=$(cat $file | cut -f 1 | sort | uniq | wc -l) 

	echo "${sample} ${reads}" | awk >> atka_mapped_reads_concatenate.tsv
done

