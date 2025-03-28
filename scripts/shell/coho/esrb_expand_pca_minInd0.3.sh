#!/bin/bash

#SBATCH --job-name=esrbPCA
#SBATCH --time=0-1:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/esrb_pca.%j.out # update your out file directory
#SBATCH --error=/home/nhowe/runtiming/coho/job_outfiles/esrb_pca.%j.err # update your error readout directory

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

path=/home/nhowe/runtiming/coho
prefix=coho-chum
chrom=NC_068449.1
firstpos=25414060 #expanded region based on boundary SNPs within peak that had fst > 0.9
lastpos=25501622 #expanded region based on boundary SNPs within peak that had fst > 0.9

# call in beagle file for pop and chromosome
BEAGLE=${path}/gls/${prefix}_${chrom}_minInd0.3_minDepthHalf.beagle.gz

# cut 9a and 9b for the regions of interest
zcat ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_minDepthHalf_esrb.beagle.gz

pcangsd.py -threads 10 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_minDepthHalf_esrb.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}_esrb_minInd0.3_minDepthHalf
