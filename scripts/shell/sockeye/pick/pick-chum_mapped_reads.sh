#!/bin/bash

#SBATCH --job-name=pickReads
#SBATCH --cpus-per-task=1
#SBATCH --time=2-24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_mapped_reads.%j.out # update your out file directory

module unload bio/samtools/1.19
module load bio/samtools/1.19

cd /scratch2/nhowe/runtiming/sock-pick/bamtools

> /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum_mapped_reads.csv

for i in /scratch2/nhowe/runtiming/sock-pick/bamtools/pick-chum*.bam;
do
	samtools view -F 0x4 ${i} | cut -f 1 | sort | uniq | wc -l >> /scratch2/nhowe/runtiming/sock-pick/pick-chum_mapped_reads.csv
done

