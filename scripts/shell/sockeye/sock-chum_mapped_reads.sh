#!/bin/bash

#SBATCH --job-name=sockReads
#SBATCH --cpus-per-task=10
#SBATCH --time=2-24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/sockeye/job_outfiles/sock-chum_mapped_reads.%j.out # update your out file directory

module unload bio/samtools/1.19
module load bio/samtools/1.19

cd /home/nhowe/runtiming/sockeye/bamtools

touch sock-chum_mapped_reads.csv

for i in *.bam;
do
	samtools view -F 0x4 ${i} | cut -f 1 | sort | uniq | wc -l >> sock-chum_mapped_reads.csv
done

