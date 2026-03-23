#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=tocoho_concat-mafs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=medmem
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minDepthHalf_concatenate-mafs_%A.out

module unload bio/angsd/0.933
module load bio/angsd/0.933

for i in /home/nhowe/runtiming/tocoho/gls/tocoho_*_minDepthHalf.mafs.gz 
do 
	zcat $i | tail -n +2 -q >> /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.mafs 
done

cut -f 1,2,3,4 /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.mafs > /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites
gzip /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.mafs

angsd sites index /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites
