#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --job-name=tocoho_concat-beagles
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minDepthHalf_concatenate-beagles_%A.out

zcat /home/nhowe/runtiming/tocoho/gls/tocoho_NC_068421.1_minDepthHalf.beagle.gz | head -n 1 > /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.beagle 

for i in /home/nhowe/runtiming/tocoho/gls/tocoho_*_minDepthHalf.beagle.gz 
	
do zcat $i | tail -n +2 -q >> /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.beagle
done

gzip /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.beagle
