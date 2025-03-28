#!/bin/bash

#SBATCH --job-name=lrrc9PCA
#SBATCH --time=0-1:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/sockeye/job_outfiles/lrrc9_pca.%j.out # update your out file directory
#SBATCH --error=/home/nhowe/runtiming/sockeye/job_outfiles/lrrc9_pca.%j.err # update your error readout directory

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

path=/home/nhowe/runtiming/sockeye
prefix=sock-chum
chrom=NC_068455.1
firstpos=28128954
lastpos=28169980

# call in beagle file for pop and chromosome
BEAGLE=${path}/gls/${prefix}_${chrom}_minInd0.3.beagle.gz

# cut 9a and 9b for the regions of interest
zcat ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_lrrc9_minInd0.3.beagle.gz

pcangsd.py -threads 10 -beagle ${path}/gls/${prefix}_${chrom}_lrrc9_minInd0.3.beagle.gz -o ${path}/pca/${prefix}_${chrom}_lrrc9_minInd0.3 
